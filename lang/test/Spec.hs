-- | Tests of the compiler: gate counts of the lowering rules, the Sudoku
-- circuit on right and wrong solutions, error messages, and a randomized
-- differential test (random expressions, compiled with and without
-- optimization, against direct evaluation mod r).
module Main (main) where

import Circuit
import Compile
import Control.Monad (forM_, unless, when)
import Data.IORef
import Data.List (isInfixOf)
import qualified Data.Map.Strict as M
import Parser
import System.Exit (exitFailure)

compileSrc :: Bool -> String -> Either String Circuit
compileSrc opt src = parseProgram "test" src >>= compile (Options opt)

holds :: Circuit -> [Integer] -> Bool
holds c ins = null (failing c (witness c ins))

main :: IO ()
main = do
  failures <- newIORef (0 :: Int)
  let check name ok = unless ok $ do
        putStrLn ("FAIL: " ++ name)
        modifyIORef failures (+ 1)
      gates src = either (const (-1)) (length . cGates) (compileSrc True src)

  -- lowering rules: what one gate can hold
  check "a product is one gate" $ gates "circuit t(private x, private y, public z) { assert x * y == z; }" == 1
  check "(x-1)(x-2) = 0 is one gate" $ gates "circuit t(private x) { assert (x - 1) * (x - 2) == 0; }" == 1
  check "three variables fit one gate" $ gates "circuit t(private a, private b, private c) { assert a + 2 * b - c == 7; }" == 1
  check "five variables need three gates" $ gates "circuit t(private a[5]) { assert sum(a) == 0; }" == 3
  check "x != y is two gates" $ gates "circuit t(private x, private y) { assert x != y; }" == 2
  check "constants fold away" $ gates "circuit t(private x) { let k = 3 * 4 - 12; assert x * k + x == x; }" == 0

  -- Sudoku
  src <- readFile "examples/sudoku.veil"
  inputs <- readFile "examples/sudoku.inputs"
  let rows = [map read (words (drop 1 (dropWhile (/= ':') l))) | l <- lines inputs] :: [[Integer]]
      [puzzle, solution] = rows
  forM_ [True, False] $ \opt -> do
    let name = if opt then "optimized" else "-O0"
    case compileSrc opt src of
      Left e -> check ("sudoku compiles " ++ name ++ ": " ++ e) False
      Right c -> do
        putStrLn ("sudoku " ++ name ++ ": " ++ summary c)
        check ("sudoku solution " ++ name) $ holds c (puzzle ++ solution)
        let swap = take 2 solution
            given = reverse swap ++ drop 2 solution
            dup = take 2 solution ++ [5] ++ drop 3 solution
            outOfRange = map (\d -> if d == 1 then 10 else d) solution
        check ("a changed given " ++ name) $ not (holds c (puzzle ++ given))
        check ("a repeated digit " ++ name) $ not (holds c (puzzle ++ dup))
        check ("digits outside 1..9 " ++ name) $ not (holds c (puzzle ++ outOfRange))
        check ("an empty puzzle accepts the solution " ++ name) $ holds c (replicate 81 0 ++ solution)
  let o = either (const 0) (length . cGates) (compileSrc True src)
      u = either (const 0) (length . cGates) (compileSrc False src)
  check "optimization saves gates" (o < u)

  -- interpolation
  case compileSrc True "circuit t(private x, public y) { assert interp(x, 1, [1, 2, 4, 8]) == y; }" of
    Left e -> check ("interp: " ++ e) False
    Right c -> do
      check "interp at 3" $ holds c [3, 4]
      check "interp at 4" $ holds c [4, 8]
      check "interp wrong" $ not (holds c [3, 5])

  -- errors carry a position and a reason
  let err src' needle = case compileSrc True src' of
        Left e -> check ("error mentions " ++ show needle ++ ": " ++ e) (needle `isInfixOf` e)
        Right _ -> check ("expected an error: " ++ needle) False
  err "circuit t(private a[3]) {\n  assert a[3] == 0;\n}" "line 2"
  err "circuit t(private a[3]) { assert a[3] == 0; }" "out of range"
  err "circuit t(private x) { assert y == 0; }" "unknown name y"
  err "circuit t(private x) { assert 1 == 2; }" "always false"
  err "circuit t(private x) { for i in 0..x { assert x == 0; } }" "compile time"
  err "def f(a) { return f(a); } circuit t(private x) { assert f(x) == 0; }" "too deep"
  err "circuit t(private x) { let x = 1; }" "already defined"
  check "parse errors are reported" $ either (const True) (const False) (compileSrc True "circuit t( { }")

  -- randomized differential test
  seed <- newIORef (12345 :: Integer)
  let rand k = do
        s <- readIORef seed
        let s' = (s * 6364136223846793005 + 1442695040888963407) `mod` (2 ^ (64 :: Int))
        writeIORef seed s'
        pure (fromIntegral ((s' `div` 65536) `mod` fromIntegral k) :: Int)
      genExpr :: Int -> IO (String, M.Map String Integer -> Integer)
      genExpr 0 = do
        t <- rand 4
        case t of
          0 -> do
            k <- rand 20
            pure (show k, const (fromIntegral k))
          _ -> do
            v <- rand 3
            let n = ["a", "b", "c"] !! v
            pure (n, (M.! n))
      genExpr d = do
        t <- rand 5
        if t == 0
          then genExpr 0
          else do
            (s1, f1) <- genExpr (d - 1)
            (s2, f2) <- genExpr (d - 1)
            o <- rand 3
            pure $ case o of
              0 -> ("(" ++ s1 ++ " + " ++ s2 ++ ")", \e -> f1 e + f2 e)
              1 -> ("(" ++ s1 ++ " - " ++ s2 ++ ")", \e -> f1 e - f2 e)
              _ -> ("(" ++ s1 ++ " * " ++ s2 ++ ")", \e -> f1 e * f2 e)
  forM_ [1 .. 300 :: Int] $ \i -> do
    (s, f) <- genExpr 4
    vals <- mapM (\_ -> fromIntegral <$> rand 1000) [1 .. 3 :: Int]
    let env = M.fromList (zip ["a", "b", "c"] vals)
        out = fr (f env)
        prog = "circuit t(private a, private b, private c, public out) { assert " ++ s ++ " == out; }"
    forM_ [True, False] $ \opt -> case compileSrc opt prog of
      Left e
        | "always false" `isInfixOf` e -> pure ()
        | otherwise -> check ("random " ++ show i ++ ": " ++ e ++ " in " ++ s) False
      Right c -> do
        check ("random " ++ show i ++ " holds: " ++ s) (holds c (vals ++ [out]))
        -- a wrong output fails unless the expression does not involve out (then 0 gates)
        when (not (null (cGates c))) $
          check ("random " ++ show i ++ " rejects: " ++ s) (not (holds c (vals ++ [out + 1])))
  n <- readIORef failures
  if n == 0 then putStrLn "all compiler tests passed" else putStrLn (show n ++ " failures") >> exitFailure
