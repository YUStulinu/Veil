-- | The parser (Parsec). Grammar, informally:
--
-- > program  := (def | circuit)*            exactly one circuit
-- > def      := "def" name "(" names ")" block
-- > circuit  := "circuit" name "(" params ")" block
-- > param    := ("public" | "private") name ("[" int "]")?
-- > block    := "{" stmt* "}"
-- > stmt     := "let" name "=" expr ";" | name "=" expr ";"
-- >           | "assert" expr ("==" | "!=") expr ";"
-- >           | "for" name "in" expr ".." expr block
-- >           | "return" expr ";" | expr ";"
-- > expr     := the usual + - * / % with unary minus, calls f(..), indexing a[i],
-- >             array literals [a, b] and comprehensions [e for i in lo..hi]
--
-- Comments run from // to the end of the line.
module Parser (parseProgram) where

import Control.Monad (void)
import Data.Functor.Identity (Identity)
import Syntax
import Text.Parsec
import Text.Parsec.Expr
import qualified Text.Parsec.Token as T

type P = Parsec String ()

lexer :: T.TokenParser ()
lexer =
  T.makeTokenParser
    T.LanguageDef
      { T.commentStart = "/*",
        T.commentEnd = "*/",
        T.commentLine = "//",
        T.nestedComments = False,
        T.identStart = letter <|> char '_',
        T.identLetter = alphaNum <|> char '_',
        T.opStart = oneOf "+-*/%=!.",
        T.opLetter = oneOf "=.",
        T.reservedNames = ["def", "circuit", "public", "private", "let", "assert", "for", "in", "return"],
        T.reservedOpNames = ["+", "-", "*", "/", "%", "=", "==", "!=", ".."],
        T.caseSensitive = True
      }

ident :: P String
ident = T.identifier lexer

reserved :: String -> P ()
reserved = T.reserved lexer

op :: String -> P ()
op = T.reservedOp lexer

sym :: String -> P ()
sym = void . T.symbol lexer

integer :: P Integer
integer = T.natural lexer

parseProgram :: FilePath -> String -> Either String Program
parseProgram file src = either (Left . show) Right (parse (T.whiteSpace lexer *> program <* eof) file src)

program :: P Program
program = do
  items <- many (Left <$> def <|> Right <$> circuit)
  case [c | Right c <- items] of
    [(name, ps, body)] -> pure (Program [f | Left f <- items] name ps body)
    [] -> fail "no circuit"
    _ -> fail "more than one circuit"

def :: P Func
def = do
  reserved "def"
  Func <$> ident <*> T.parens lexer (T.commaSep lexer ident) <*> block

circuit :: P (String, [Param], [Stmt])
circuit = do
  reserved "circuit"
  (,,) <$> ident <*> T.parens lexer (T.commaSep lexer param) <*> block

param :: P Param
param = do
  pub <- (True <$ reserved "public") <|> (False <$ reserved "private")
  name <- ident
  size <- optionMaybe (T.brackets lexer integer)
  pure (Param name pub size)

block :: P [Stmt]
block = T.braces lexer (many stmt)

stmt :: P Stmt
stmt = do
  pos <- getPosition
  Stmt pos <$> kind
  where
    semi = sym ";"
    kind =
      (reserved "let" *> (Let <$> ident <*> (op "=" *> expr)) <* semi)
        <|> (reserved "assert" *> (Assert <$> expr <*> rel <*> expr) <* semi)
        <|> ( reserved "for" *> do
                i <- ident
                reserved "in"
                lo <- expr
                op ".."
                hi <- expr
                For i lo hi <$> block
            )
        <|> (reserved "return" *> (Return <$> expr) <* semi)
        <|> try (Assign <$> ident <*> (op "=" *> expr) <* semi)
        <|> (Do <$> expr <* semi)
    rel = (REq <$ op "==") <|> (RNeq <$ op "!=")

expr :: P Expr
expr = buildExpressionParser table term
  where
    table :: OperatorTable String () Identity Expr
    table =
      [ [Prefix (Neg <$ op "-")],
        [Infix (Bin Mul <$ op "*") AssocLeft, Infix (Bin Div <$ op "/") AssocLeft, Infix (Bin Mod <$ op "%") AssocLeft],
        [Infix (Bin Add <$ op "+") AssocLeft, Infix (Bin Sub <$ op "-") AssocLeft]
      ]

term :: P Expr
term = do
  base <- atom
  idxs <- many (T.brackets lexer expr)
  pure (foldl Index base idxs)

atom :: P Expr
atom =
  T.parens lexer expr
    <|> (Num <$> integer)
    <|> array
    <|> do
      name <- ident
      maybe (Name name) (Call name) <$> optionMaybe (T.parens lexer (T.commaSep lexer expr))

array :: P Expr
array = T.brackets lexer $ do
  first <- optionMaybe expr
  case first of
    Nothing -> pure (ArrayLit [])
    Just e ->
      ( do
          reserved "for"
          i <- ident
          reserved "in"
          lo <- expr
          op ".."
          Comprehension e i lo <$> expr
      )
        <|> (ArrayLit . (e :) <$> many (sym "," *> expr))
