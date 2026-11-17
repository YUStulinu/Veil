-- | The compiler's output: PLONK gates over variables, the witness program
-- (hints), and the text format the Rust prover reads. Also a witness
-- evaluator, so the compiler's own tests can check circuits.
module Circuit where

import Data.List (intercalate)
import qualified Data.Map.Strict as M

-- | the BLS12-381 scalar field
modulus :: Integer
modulus = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001

fr :: Integer -> Integer
fr x = x `mod` modulus

-- | modular inverse (0 for 0), by Fermat
inv :: Integer -> Integer
inv 0 = 0
inv x = powMod (fr x) (modulus - 2)

powMod :: Integer -> Integer -> Integer
powMod _ 0 = 1
powMod b e
  | even e = let h = powMod b (e `div` 2) in h * h `mod` modulus
  | otherwise = b * powMod b (e - 1) `mod` modulus

type Var = Int

data Hint
  = HInput Int
  | HEval Integer Integer Integer Integer Var Var -- out = qM a b + qL a + qR b + qC
  | HInv Var
  deriving (Show, Eq)

-- | qM a b + qL a + qR b + qO c + qC = 0
data Gate = Gate {qM, qL, qR, qO, qC :: Integer, wa, wb, wc :: Var}
  deriving (Show, Eq)

data InputDecl = InputDecl {inName :: String, inPublic :: Bool, inLen :: Int}
  deriving (Show)

data Circuit = Circuit
  { cName :: String,
    cVars :: Int,
    cInputs :: [InputDecl],
    cPublic :: [Var],
    cHints :: [(Var, Hint)],
    cGates :: [Gate]
  }
  deriving (Show)

-- | constants as decimal; large ones that are small negatives as "-k"
showFr :: Integer -> String
showFr x
  | y > modulus `div` 2 = show (y - modulus)
  | otherwise = show y
  where
    y = fr x

emit :: Circuit -> String
emit c =
  unlines $
    ["veil-circuit 1", "name " ++ cName c, "vars " ++ show (cVars c)]
      ++ [unwords ["input", inName d, if inPublic d then "public" else "private", show (inLen d)] | d <- cInputs c]
      ++ [unwords ("public" : show (length (cPublic c)) : map show (cPublic c))]
      ++ map hint (cHints c)
      ++ map gate (cGates c)
  where
    hint (o, HInput k) = unwords ["hint", show o, "input", show k]
    hint (o, HInv a) = unwords ["hint", show o, "inv", show a]
    hint (o, HEval m l r k a b) = unwords (["hint", show o, "eval"] ++ map showFr [m, l, r, k] ++ [show a, show b])
    gate g = unwords ("gate" : map showFr [qM g, qL g, qR g, qO g, qC g] ++ map show [wa g, wb g, wc g])

-- | all variable values from the flattened inputs
witness :: Circuit -> [Integer] -> M.Map Var Integer
witness c inputs = foldl step (M.singleton 0 0) (cHints c)
  where
    ins = M.fromList (zip [0 ..] (map fr inputs))
    step w (o, h) = M.insert o (value w h) w
    value _ (HInput k) = ins M.! k
    value w (HInv a) = inv (w M.! a)
    value w (HEval m l r k a b) = let (x, y) = (w M.! a, w M.! b) in fr (m * x * y + l * x + r * y + k)

-- | the indices of the gates that fail
failing :: Circuit -> M.Map Var Integer -> [Int]
failing c w = [i | (i, g) <- zip [0 ..] (cGates c), holds g /= 0]
  where
    holds g = let (a, b, x) = (w M.! wa g, w M.! wb g, w M.! wc g) in fr (qM g * a * b + qL g * a + qR g * b + qO g * x + qC g)

summary :: Circuit -> String
summary c =
  intercalate
    ", "
    [ show (length (cGates c)) ++ " gates",
      show (length (cPublic c)) ++ " public inputs",
      show (cVars c) ++ " variables",
      show (length [() | (_, HInv _) <- cHints c]) ++ " inverse hints"
    ]
