-- | From the syntax tree to gates. Functions are inlined and loops unrolled
-- at compile time; integers known at compile time stay exact integers
-- (indices, bounds) and everything else is a value over circuit variables:
--
--   * a linear combination  c_1 v_1 + ... + c_k v_k + c_0      (no gate yet)
--   * a quadratic term       k a b + (linear part)               (no gate yet)
--
-- Gates appear only when a value has to become a variable: as the operand
-- of a multiplication, in an assertion, or when it does not fit one gate.
-- A PLONK gate holds q_M a b + q_L a + q_R b + q_O c + q_C, so a product of
-- two affine values, or an assertion over at most three variables, costs
-- one gate. With optimization off (-O0) every product is a gate of its own
-- and nothing is shared, which shows what the two rules save.
module Compile (compile, Options (..), defaultOptions) where

import Circuit
import Control.Monad (foldM, forM, unless, when)
import Control.Monad.Except (catchError, throwError)
import Control.Monad.State.Strict
import qualified Data.Map.Strict as M
import Syntax
import Text.Parsec.Pos (sourceColumn, sourceLine)

newtype Options = Options {optimize :: Bool}

defaultOptions :: Options
defaultOptions = Options True

-- | coefficients (mod r) and a constant
data Lin = Lin (M.Map Var Integer) Integer
  deriving (Show)

data Val
  = VC Integer -- known at compile time (an exact integer)
  | VL Lin
  | VQ Integer Var Var Lin -- k a b + lin
  | VArr [Val]
  deriving (Show)

data St = St
  { sNext :: Int,
    sHints :: [(Var, Hint)],
    sGates :: [Gate],
    sCse :: M.Map (Integer, Integer, Integer, Integer, Var, Var) Var,
    sScopes :: [M.Map String Val],
    sFuncs :: M.Map String Func,
    sDepth :: Int,
    sOpt :: Bool
  }

type C = StateT St (Either String)

compile :: Options -> Program -> Either String Circuit
compile opts prog = evalStateT go st0
  where
    st0 = St 1 [] [] M.empty [M.empty] (M.fromList [(fName f, f) | f <- progFuncs prog]) 0 (optimize opts)
    go = do
      -- the inputs, in declaration order
      (decls, binds, pubs) <- unzip3 <$> declare 0 (progParams prog)
      modify $ \s -> s {sScopes = [M.fromList binds]}
      r <- execBlock (progBody prog)
      case r of
        Just _ -> throwError "a circuit does not return a value"
        Nothing -> pure ()
      st <- get
      pure
        Circuit
          { cName = progName prog,
            cVars = sNext st,
            cInputs = decls,
            cPublic = concat pubs,
            cHints = reverse (sHints st),
            cGates = reverse (sGates st)
          }
    declare _ [] = pure []
    declare k (p : ps) = do
      let n = maybe 1 fromIntegral (pSize p)
      vs <- forM [0 .. n - 1] $ \i -> newVar (HInput (k + i))
      let val = case (pSize p, vs) of
            (Nothing, [v]) -> var v
            _ -> VArr (map var vs)
      rest <- declare (k + n) ps
      pure ((InputDecl (pName p) (pPublic p) n, (pName p, val), if pPublic p then vs else []) : rest)

-- ---- values ----

var :: Var -> Val
var v = VL (Lin (M.singleton v 1) 0)

norm :: Lin -> Val
norm (Lin m k) = let m' = M.filter (/= 0) (M.map fr m) in if M.null m' then VC (fr k) else VL (Lin m' (fr k))

addLin :: Lin -> Lin -> Lin
addLin (Lin a k) (Lin b j) = Lin (M.unionWith (+) a b) (k + j)

scaleLin :: Integer -> Lin -> Lin
scaleLin c (Lin m k) = Lin (M.map (* c) m) (k * c)

asLin :: Val -> Maybe Lin
asLin (VC c) = Just (Lin M.empty c)
asLin (VL l) = Just l
asLin _ = Nothing

addV :: Val -> Val -> C Val
addV (VC a) (VC b) = pure (VC (a + b))
addV (VArr _) _ = throwError "cannot add arrays"
addV _ (VArr _) = throwError "cannot add arrays"
addV (VQ k a b l) o | Just l2 <- asLin o = pure (VQ k a b (addLin l l2))
addV o q@(VQ {}) | Just _ <- asLin o = addV q o
addV q@(VQ {}) o@(VQ {}) = do
  v <- materialize o
  addV q (var v)
addV x y = case (asLin x, asLin y) of
  (Just a, Just b) -> pure (norm (addLin a b))
  _ -> throwError "internal: addV"

scaleV :: Integer -> Val -> C Val
scaleV c (VC x) = pure (VC (c * x))
scaleV c (VL l) = pure (norm (scaleLin c l))
scaleV c (VQ k a b l)
  | fr c == 0 = pure (VC 0)
  | otherwise = pure (VQ (fr (c * k)) a b (scaleLin c l))
scaleV _ (VArr _) = throwError "cannot scale an array"

mulV :: Val -> Val -> C Val
mulV (VC a) (VC b) = pure (VC (a * b))
mulV (VC k) v = scaleV k v
mulV v (VC k) = scaleV k v
mulV (VArr _) _ = throwError "cannot multiply arrays"
mulV _ (VArr _) = throwError "cannot multiply arrays"
mulV x y = do
  opt <- gets sOpt
  if opt
    then do
      (c1, a, d1) <- affine x
      (c2, b, d2) <- affine y
      pure (VQ (fr (c1 * c2)) a b (addLin (Lin (M.singleton a (c1 * d2)) 0) (Lin (M.singleton b (c2 * d1)) (d1 * d2))))
    else do
      a <- materialize x
      b <- materialize y
      var <$> emitEval 1 0 0 0 a b

-- | c v + d, making a variable first if the value has more than one
affine :: Val -> C (Integer, Var, Integer)
affine (VL (Lin m k)) | [(v, c)] <- M.toList m = pure (c, v, k)
affine v = do
  x <- materialize v
  pure (1, x, 0)

-- | a variable holding the value
materialize :: Val -> C Var
materialize (VC c)
  | fr c == 0 = pure 0
  | otherwise = emitEval 0 0 0 c 0 0
materialize (VArr _) = throwError "an array is not a single value"
materialize (VL (Lin m k)) = case M.toList (M.filter ((/= 0) . fr) m) of
  [] -> materialize (VC k)
  [(v, c)] | fr c == 1 && fr k == 0 -> pure v
  [(v, c)] -> emitEval 0 c 0 k v 0
  [(v, c), (u, d)] -> emitEval 0 c d k v u
  ((v, c) : (u, d) : rest) -> do
    t <- emitEval 0 c d 0 v u
    materialize (VL (Lin (M.insert t 1 (M.fromList rest)) k))
materialize (VQ k a b (Lin m c0)) = do
  let ca = fr (M.findWithDefault 0 a m)
      cb = if a == b then 0 else fr (M.findWithDefault 0 b m)
      rest = M.filter ((/= 0) . fr) (M.delete a (M.delete b m))
  if M.null rest
    then emitEval k ca cb c0 a b
    else do
      t <- emitEval k ca cb c0 a b
      materialize (VL (Lin (M.insert t 1 rest) 0))

newVar :: Hint -> C Var
newVar h = do
  s <- get
  let v = sNext s
  put s {sNext = v + 1, sHints = (v, h) : sHints s}
  pure v

gate :: Gate -> C ()
gate g = modify $ \s -> s {sGates = g : sGates s}

-- | out = m a b + l a + r b + k, shared with an identical earlier one when optimizing
emitEval :: Integer -> Integer -> Integer -> Integer -> Var -> Var -> C Var
emitEval m l r k a b = do
  let (l', r', a', b') = if a <= b then (fr l, fr r, a, b) else (fr r, fr l, b, a)
      key = (fr m, l', r', fr k, a', b')
  st <- get
  case M.lookup key (sCse st) of
    Just v | sOpt st -> pure v
    _ -> do
      o <- newVar (HEval (fr m) l' r' (fr k) a' b')
      gate (Gate (fr m) l' r' (fr (-1)) (fr k) a' b' o)
      modify $ \s -> s {sCse = M.insert key o (sCse s)}
      pure o

assertZero :: Val -> C ()
assertZero (VC c) = when (fr c /= 0) (throwError "this assertion is always false")
assertZero (VArr _) = throwError "cannot assert on an array"
assertZero (VL (Lin m k)) = case M.toList (M.filter ((/= 0) . fr) m) of
  ((v, c) : (u, d) : rest@(_ : _ : _)) -> do
    t <- emitEval 0 c d 0 v u
    assertZero (VL (Lin (M.insert t 1 (M.fromList rest)) k))
  ts -> do
    let pad = take 3 (ts ++ repeat (0, 0))
        [(v1, c1), (v2, c2), (v3, c3)] = pad
    gate (Gate 0 (fr c1) (fr c2) (fr c3) (fr k) v1 v2 v3)
assertZero (VQ k a b (Lin m c0)) = do
  let ca = fr (M.findWithDefault 0 a m)
      cb = if a == b then 0 else fr (M.findWithDefault 0 b m)
      rest = M.filter ((/= 0) . fr) (M.delete a (M.delete b m))
  (vr, cr) <- case M.toList rest of
    [] -> pure (0, 0)
    [(v, c)] -> pure (v, c)
    _ -> do
      t <- materialize (VL (Lin rest 0))
      pure (t, 1)
  gate (Gate k ca cb (fr cr) (fr c0) a b vr)

assertNonZero :: Val -> C ()
assertNonZero (VC c) = when (fr c == 0) (throwError "this assertion is always false")
assertNonZero v = do
  x <- materialize v
  i <- newVar (HInv x)
  gate (Gate 1 0 0 0 (fr (-1)) x i 0)

-- ---- statements ----

lookupName :: String -> C Val
lookupName n = do
  sc <- gets sScopes
  case [v | Just v <- map (M.lookup n) sc] of
    (v : _) -> pure v
    [] -> throwError ("unknown name " ++ n)

withScope :: M.Map String Val -> C a -> C a
withScope m act = do
  modify $ \s -> s {sScopes = m : sScopes s}
  r <- act
  modify $ \s -> s {sScopes = drop 1 (sScopes s)}
  pure r

execBlock :: [Stmt] -> C (Maybe Val)
execBlock [] = pure Nothing
execBlock (s : ss) = do
  r <- exec s
  case r of
    Just v -> pure (Just v)
    Nothing -> execBlock ss

exec :: Stmt -> C (Maybe Val)
exec (Stmt pos k) = run k `catchError` \e ->
  throwError (if take 5 e == "line " then e else "line " ++ show (sourceLine pos) ++ ", column " ++ show (sourceColumn pos) ++ ": " ++ e)
  where
    run (Let n e) = do
      v <- eval e
      sc <- gets sScopes
      case sc of
        (top : rest) -> do
          when (M.member n top) (throwError (n ++ " is already defined here"))
          modify $ \s -> s {sScopes = M.insert n v top : rest}
        [] -> throwError "internal: no scope"
      pure Nothing
    run (Assign n e) = do
      v <- eval e
      sc <- gets sScopes
      let (before, after) = break (M.member n) sc
      case after of
        [] -> throwError ("assignment to an undefined name " ++ n)
        (m : rest) -> modify $ \s -> s {sScopes = before ++ M.insert n v m : rest}
      pure Nothing
    run (Assert a r b) = do
      d <- join2 (\x y -> scaleV (-1) y >>= addV x) a b
      case r of
        REq -> assertZero d
        RNeq -> assertNonZero d
      pure Nothing
    run (For i lo hi body) = do
      l <- constant "a loop bound" =<< eval lo
      h <- constant "a loop bound" =<< eval hi
      loop l h
      where
        loop c h
          | c >= h = pure Nothing
          | otherwise = do
              r <- withScope (M.singleton i (VC c)) (execBlock body)
              maybe (loop (c + 1) h) (pure . Just) r
    run (Do e) = Nothing <$ evalCall e
    run (Return e) = Just <$> eval e
    join2 f a b = do
      x <- eval a
      y <- eval b
      f x y
    evalCall (Call f args) = callFunc f args
    evalCall e = Just <$> eval e

constant :: String -> Val -> C Integer
constant _ (VC c) = pure c
constant what _ = throwError (what ++ " must be known at compile time")

eval :: Expr -> C Val
eval (Num n) = pure (VC n)
eval (Name n) = lookupName n
eval (Neg e) = eval e >>= scaleV (-1)
eval (Bin o a b) = do
  x <- eval a
  y <- eval b
  case o of
    Add -> addV x y
    Sub -> scaleV (-1) y >>= addV x
    Mul -> mulV x y
    Div -> intOp "/" div x y
    Mod -> intOp "%" mod x y
  where
    intOp name f x y = do
      p <- constant ("an operand of " ++ name) x
      q <- constant ("an operand of " ++ name) y
      when (q == 0) (throwError "division by zero")
      pure (VC (f p q))
eval (Index a i) = do
  arr <- eval a
  k <- constant "an index" =<< eval i
  case arr of
    VArr vs
      | k >= 0 && k < fromIntegral (length vs) -> pure (vs !! fromIntegral k)
      | otherwise -> throwError ("index " ++ show k ++ " out of range (length " ++ show (length vs) ++ ")")
    _ -> throwError "indexing something that is not an array"
eval (ArrayLit es) = VArr <$> mapM eval es
eval (Comprehension e i lo hi) = do
  l <- constant "a range bound" =<< eval lo
  h <- constant "a range bound" =<< eval hi
  VArr <$> forM [l .. h - 1] (\c -> withScope (M.singleton i (VC c)) (eval e))
eval (Call f args) = callFunc f args >>= maybe (throwError (f ++ " returns no value")) pure

callFunc :: String -> [Expr] -> C (Maybe Val)
callFunc "sum" [a] = do
  arr <- eval a
  case arr of
    VArr vs -> Just <$> foldM addV (VC 0) vs
    _ -> throwError "sum takes an array"
callFunc "len" [a] = do
  arr <- eval a
  case arr of
    VArr vs -> pure (Just (VC (fromIntegral (length vs))))
    _ -> throwError "len takes an array"
callFunc "interp" [x, lo, ys] = do
  xv <- eval x
  l <- constant "interp's start" =<< eval lo
  yv <- eval ys
  pts <- case yv of
    VArr vs -> mapM (constant "interp's values") vs
    _ -> throwError "interp takes an array of constants"
  -- Horner on the interpolating polynomial
  let cs = interpolate [(l + i, y) | (i, y) <- zip [0 ..] pts]
  Just <$> foldM (\acc c -> mulV acc xv >>= addV (VC c)) (VC 0) (reverse cs)
callFunc f args = do
  funcs <- gets sFuncs
  fn <- maybe (throwError ("unknown function " ++ f)) pure (M.lookup f funcs)
  unless (length args == length (fParams fn)) $
    throwError (f ++ " takes " ++ show (length (fParams fn)) ++ " arguments")
  vals <- mapM eval args
  d <- gets sDepth
  when (d > 200) (throwError "calls nested too deep (recursion?)")
  modify $ \s -> s {sDepth = d + 1}
  -- a function sees its parameters only, not the caller's locals
  saved <- gets sScopes
  modify $ \s -> s {sScopes = [M.fromList (zip (fParams fn) vals)]}
  r <- execBlock (fBody fn)
  modify $ \s -> s {sScopes = saved, sDepth = d}
  pure r

-- | the coefficients (lowest first) of the polynomial through the points, mod r
interpolate :: [(Integer, Integer)] -> [Integer]
interpolate pts = map fr (foldr1 (zipWith (+)) [map (* (fr y * inv (denom i))) (basis i) | (i, (_, y)) <- zip [0 ..] pts])
  where
    xs = map fst pts
    others i = [xj | (j, xj) <- zip [0 :: Int ..] xs, j /= i]
    basis i = foldl mulLinear [1] (others i)
    denom i = fr (product [xs !! i - xj | xj <- others i])
    -- p(X) (X - c)
    mulLinear p c = zipWith (+) (0 : p) (map (* (-c)) p ++ [0])
