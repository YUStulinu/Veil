-- | The abstract syntax of the Veil circuit language.
module Syntax where

import Text.Parsec.Pos (SourcePos)

data Expr
  = Num Integer
  | Name String
  | Index Expr Expr
  | Bin Op Expr Expr
  | Neg Expr
  | Call String [Expr]
  | ArrayLit [Expr]
  | Comprehension Expr String Expr Expr -- [e for i in lo..hi]
  deriving (Show, Eq)

data Op = Add | Sub | Mul | Div | Mod
  deriving (Show, Eq)

data Rel = REq | RNeq
  deriving (Show, Eq)

data Stmt = Stmt SourcePos StmtKind
  deriving (Show)

data StmtKind
  = Let String Expr
  | Assign String Expr
  | Assert Expr Rel Expr
  | For String Expr Expr [Stmt] -- for i in lo..hi (hi excluded)
  | Do Expr -- a call for its constraints
  | Return Expr
  deriving (Show)

data Param = Param
  { pName :: String,
    pPublic :: Bool,
    pSize :: Maybe Integer -- Nothing: a scalar
  }
  deriving (Show)

data Func = Func
  { fName :: String,
    fParams :: [String],
    fBody :: [Stmt]
  }
  deriving (Show)

data Program = Program
  { progFuncs :: [Func],
    progName :: String,
    progParams :: [Param],
    progBody :: [Stmt]
  }
  deriving (Show)
