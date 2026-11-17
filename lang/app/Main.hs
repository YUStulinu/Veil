-- | veilc: compiles a Veil circuit to the gate list the Rust prover reads.
--
-- > veilc circuit.veil [-o out.circuit] [-O0] [--check inputs.txt]
--
-- --check runs the witness program on the inputs ("name: v1 v2 ..." lines)
-- and reports the first failing gate, without any cryptography.
module Main (main) where

import Circuit
import Compile
import qualified Data.Map.Strict as M
import Parser
import System.Environment (getArgs)
import System.Exit (exitFailure)
import System.IO (hPutStrLn, stderr)

data Args = Args {aIn :: Maybe FilePath, aOut :: Maybe FilePath, aOpt :: Bool, aCheck :: Maybe FilePath}

parseArgs :: [String] -> Args -> Either String Args
parseArgs [] a = Right a
parseArgs ("-o" : f : r) a = parseArgs r a {aOut = Just f}
parseArgs ("-O0" : r) a = parseArgs r a {aOpt = False}
parseArgs ("--check" : f : r) a = parseArgs r a {aCheck = Just f}
parseArgs (f : r) a | take 1 f /= "-" = parseArgs r a {aIn = Just f}
parseArgs (f : _) _ = Left ("unknown option " ++ f)

die' :: String -> IO a
die' m = hPutStrLn stderr m >> exitFailure

main :: IO ()
main = do
  args <- getArgs
  a <- either die' pure (parseArgs args (Args Nothing Nothing True Nothing))
  file <- maybe (die' "usage: veilc circuit.veil [-o out.circuit] [-O0] [--check inputs.txt]") pure (aIn a)
  src <- readFile file
  prog <- either die' pure (parseProgram file src)
  c <- either (die' . ((file ++ ": ") ++)) pure (compile (Options (aOpt a)) prog)
  hPutStrLn stderr (cName c ++ ": " ++ summary c)
  maybe (pure ()) (\o -> writeFile o (emit c)) (aOut a)
  case aCheck a of
    Nothing -> pure ()
    Just f -> do
      txt <- readFile f
      ins <- either die' pure (readInputs c txt)
      case failing c (witness c ins) of
        [] -> hPutStrLn stderr "all gates hold"
        (g : _) -> die' ("gate " ++ show g ++ " fails: " ++ show (cGates c !! g))

readInputs :: Circuit -> String -> Either String [Integer]
readInputs c txt = concat <$> mapM one (cInputs c)
  where
    table = M.fromList [(takeWhile (/= ':') l, map read (words (drop 1 (dropWhile (/= ':') l)))) | l <- lines txt, ':' `elem` l]
    one d = case M.lookup (inName d) table of
      Just vs | length vs == inLen d -> Right vs
      Just vs -> Left (inName d ++ ": " ++ show (length vs) ++ " values, expected " ++ show (inLen d))
      Nothing -> Left ("missing input " ++ inName d)
