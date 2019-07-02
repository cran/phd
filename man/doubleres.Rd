\name{doubleres}
\alias{doubleres}
\title{Permutation test based on double residualization}
\description{
  Provides a class of tests for testing in high-dimensional linear models. The tests are robust against heteroscedasticity and non-normality. They often provide good type I error control even under anti-sparsity. 
}
\usage{
  doubleres(y,X,X1,nperm=2E4,lambda="lambda.min",flip="FALSE",nfolds=10)
}
\arguments{
  \item{y}{The values of the outcome.}
  \item{X}{The design matrix. 
   If the covariate of interest is included in \code{X},  it should be included in the first column. 
   If it is not included in \code{X}, then specify \code{X1}. The data do not need to be standardized, 
   since this is automatically done by this function. Do not include a columns of 1's.}
  \item{X1}{n-vector with the (1-dimensional) covariate of interest.
   \code{X1} should only be specified if the covariate of interest is not already included in \code{X}.}
  \item{nperm}{The number of random permutations (or sign-flipping maps) used by the test}
  \item{lambda}{The penalty used in the ridge regressions. Default is \code{"lambda.min"}, 
  which means that the penalty is obtained using cross-validation. One can also enter \code{"lambda.1se"},
  which is an upward-conservative estimate of the optimal lambda.}
  \item{flip}{Default is "FALSE", which means that permutation is used. 
  If "TRUE", then sign-flipping is used.}
  \item{nfolds}{The number of folds used in the cross-validation (in case lambda is determined using cross-validation).}
}

\value{A two-sided p-value.}

\examples{

set.seed(5193)
n=30

X <- matrix(nr=n,nc=60,rnorm(n*60)) 
y <- X[,1]+X[,2]+X[,3] + rnorm(n,mean=0)   #H0: first coefficient=0. So H0 is false

doubleres(y, X, nperm=2000, lambda=100,flip="FALSE")   
}

