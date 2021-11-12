\name{hsemfit}
\alias{hsemfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting Hierarchical Structural Equation Models using
h-likelihood Approach
}

\description{The hsemfit is used to fit a hierarchical structural equation models (HSEMs) allowing
different models for multivariate reponse variables. A variety of distributions and link
functions for both response and the random effects are allowed. To call the fitting function
hsemfit, models for the mean and dispersion must be specified by hsemmodleing object
preferably created by calling the hsemmodeling function.}

\usage{
hsemfit(RespDist = "gaussian", BinomialDen = NULL, 
DataMain, MeanModel,DispersionModel = NULL, 
PhiFix = NULL, LamFix = NULL, structure = "correlated", 
mord = 0, dord = 1, convergence = 1e-05, 
Init_Corr = NULL, EstimateCorrelations = TRUE)
}
\arguments{
  \item{RespDist}{The distribution of the response is set by the option \code{RespDist}. The user can set
it to: "gaussian" (default), "binomial", "poisson", or "gamma".}
  \item{BinomialDen}{When \code{RespDist}="binomial", one should use the option \code{BinomialDen} to specify
the denominator for the binomial distribution. This should be "NULL" (default)
or a numeric vector of length equal to the length of \code{DataMain}. When specified
as \code{BinomialDen}=NULL and \code{RespDist}="binomial", the denominator is 1.
}
  \item{DataMain}{The option \code{DataMain} determines the data frame to be used (non-optional).}
  \item{MeanModel}{For the mean model, this option requries \code{DGHLMMODELING} object which
should specified by the option Model="mean".}
  \item{DispersionModel}{For the overdispersion model, this option requries \code{DGHLMMODELING} object
which should be specified by the option Model="dispersion".}
  \item{PhiFix}{The option for overdispersion parameters (phi) to be estimated or maintaned constant.
Specifying defaults such as \code{PhiFix} =NULL implies that phi is to be estimated. If not, phi is fixed at a value
specified by \code{PhiFix}.}
  \item{LamFix}{The option for random-effect variance (lambda) to be estimated or maintaned constant.
Specifying defaults such as \code{LamFix} =NULL implies that lambda is to be estimated. If not, lambda is fixed at a value
specified by \code{LamFix}.
}
  \item{structure}{The option \code{structure} determines structure of random effects. When \code{structure}="correlated" (or "shared"),
  correlated (or shared) random-effects model is specified.}
  \item{mord}{The option \code{mord} specifies the order of Laplace approximation to the marginal
likelihood for fitting the mean parameters. The choice is either 0 or 1 (default).
}
  \item{dord}{The option \code{dord} specifies the order of adjusted
profile likelihood for fitting the dispersion parameters. The choice is either 1
(default) or 2.
}
  \item{convergence}{Setting this option determines the criterion for convergence, which is computed
as the absolute difference between the values of all the estimated parameters in
the previous and current iterations. The default criterion is 1e-06.
}
  \item{Init_Corr}{Setting initial values of correlation (or shared parameters) between random effects}
  \item{EstimateCorrelations}{Correlation are estimated or fixed when \code{EstimateCorrelations}=TRUE (default) or \code{EstimateCorrelations}=FALSE.}
}
\value{
\item{res}{the output class resulted from jointfit_correlated.}
}
\examples{
data(ml2)

MM1<-hsemmodeling(Model="mean",Link="identity",
   LinPred=urge~urge1+dep1+(1|id)+(urge1|id)+(dep1|id),
   RandDist=c("gaussian","gaussian","gaussian"))
DM1<-hsemmodeling(Model="dispersion",Link = "log", 
   LinPred=phi~(1|id),RandDist=c("gaussian"))
MM2<-hsemmodeling(Model="mean",Link="identity",
   LinPred=dep~urge1+dep1+(1|id)+(urge1|id)+(dep1|id),
   RandDist=c("gaussian","gaussian","gaussian"))
DM2<-hsemmodeling(Model="dispersion",Link = "log", 
   LinPred=phi~(1|id),RandDist=c("gaussian"))
res<-hsemfit(RespDist=c("gaussian","gaussian"),DataMain=list(ml2,ml2),
   structure="independent",MeanModel=list(MM1,MM2),
   DispersionModel=list(DM1,DM2))
}


