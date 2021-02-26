\name{gvarclust}
\alias{gvarclust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Utilize subgrouping features for chain graphical VAR models.
}
\description{
Estimates the subgrouped chain graphical VAR model through LASSO estimation coupled with extended Bayesian information criterion for choosing the optimal tuning parameters. The estimation procedure makes use of community detection procedures to identify homogeneous subgroups of individuals prior to fitting a chain graphical VAR model to the data.
}
\usage{
gvarclust(data, vars, idvar, beepvar, dayvar, search = TRUE,
                    scale = TRUE, maxit.in, centerWithin = TRUE, gamma =
                    0.5, verbose = TRUE, redo = FALSE, subjectNetworks =
                    TRUE, lambda_min_kappa_fixed = 0.001,
                    lambda_min_beta_fixed = 0.001, lambda_min_kappa =
                    0.05, lambda_min_beta = lambda_min_kappa,
                    lambda_min_glasso = 0.01, dir = getwd(), btwnNetworks
                    = FALSE, mlgvclust = TRUE, indplot = FALSE, elim.groups, g.ID = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A matrix or data frame containing repeated measures (rows) on a set of variables (columns). Must not contain missing data.
}
\item{vars}{
Vectors of variables to include in the analysis
}
  \item{idvar}{
ID variable.
}
  \item{beepvar}{
Idk, honestly.
}
  \item{dayvar}{
The day variable indicator  
}
  \item{search}{
Set to TRUE to search the adjacency matrix using a conductance optimization procedure. Set to FALSE if you would like to use the default procedure used by S-GIMME (Gates et al., 2017) which subtracts the minimum value. Conductance is more liberal with cluster identification.
}
\item{scale}{
Logical, should responses be standardized before estimation?
}
  \item{maxit.in}{
Maximum number of iterations in the inner loop (computing beta)
}
\item{centerWithin}{
Logical, should subject data be within-person centered before estimating fixed effects?
}
  \item{gamma}{
The EBIC hyper-parameter. Set to 0 to use regular BIC.
}
\item{verbose}{
Do you want all the info?
}
\item{redo}{
Test function currently under works.
}
\item{subjectNetworks}{
  Do you want N = 1 networks?
}
\item{lambda_min_kappa_fixed}{
Minimum lambda value for the precision matrix in the fixed effects network
}
\item{lambda_min_beta_fixed}{
Minimum lambda value for the regression matrix in the fixed effects network
}
  \item{lambda_min_kappa}{
      Multiplier of maximal tuning parameter for kappa
  }
    \item{lambda_min_beta}{
      Multiplier of maximal tuning parameter for beta
  }
\item{lambda_min_glasso}{
idk
}
\item{dir}{
The directory where you want your files saved. Should be a path and a folder will be created called graphical VAR
}
\item{btwnNetworks}{
Do you want the between subjects networks generated? Set to FALSE by default to save on estimation time.
}
\item{mlgvclust}{
  Test function.
}
\item{indplot}{
 Do you want individual plots? Set to FALSE due to issues plotting individuals with empty lagged and contemporaneous networks.
}
\item{elim.groups}{
 Do you want the similarity matrices to remove the group-level information?
}
\item{g.ID}{
 Is there a confirmatory group structure? Currently does not function.
}
  \item{\dots}{
Arguments sent to \code{\link{graphicalVAR}}
}
}
\details{
Let y_t denote the vector of centered responses of a subject on a set of items on time point t. The graphical VAR model, using only one lag, is defined as follows:

y_t = Beta y_{y-1} + epsilon_t

In which epsilon_t is a vector of error and is independent between time points but not within time points. Within time points, the error is normally distributed with mean vector 0 and precision matrix (inverse covariance matrix) Kappa. The Beta matrix encodes the between time point interactions and the Kappa matrix encodes the within time point interactions. We aim to find a sparse solution for both Beta and Kappa, and do so by applying the LASSO algorithm as detailed by Rothman, Levina and Zhu (2010). The LASSO algorithm uses two tuning parameters, lambda_beta controlling the sparsity in Beta and lambda_kappa controlling the sparsity in Kappa. We estimate the model under a (by default) 50 by 50 grid of tuning parameters and choose the tuning parameters that optimize the extended Bayesian Information Criterion (EBIC; Chen and Chen,2008).

After estimation, the Beta and Kappa matrices can be standardized as described by Wild et al. (2010). The Kappa matrix can be standardized to partial contemporaneous correlations (PCC) as follows:

PCC(y_{i,t}, y_{j,t}) = - kappa_{ij} / ( sqrt{kappa_{ii} kappa_{jj}})

Similarly, the beta matrix can be standardized to partial directed correlations (PDC):

PDC(y_{i,t-1}, y_{j,t}) = beta_{ji} / sqrt{sigma_{jj} kappa_{ii} + beta_{ji}^2}

In which sigma is the inverse of kappa. Note that this process transposes the beta matrix. This is done because in representing a directed network it is typical to let rows indicate the node of origin and columns the node of destination.

Set \code{lambda_beta = 0} argument and \code{lambda_kappa = 0} for unregularized estimation.
}
\value{
A \code{graphicalVAR} object, which is a list containing:

\item{PCC}{The partial contemporaneous correlation network}
\item{PDC}{The  partial directed correlation network}
\item{beta}{The estimated beta matrix}
\item{kappa}{The estimated kappa matrix}
\item{EBIC}{The optimal EBIC}
\item{path}{Results of all tested tuning parameters}
\item{labels}{A vector containing the node labels}
}
\references{
Chen, J., & Chen, Z. (2008). Extended Bayesian information criteria for model selection with large model spaces. Biometrika, 95(3), 759-771.

Fentaw Abegaz and Ernst Wit (2013). Sparse time series chain graphical
models for reconstructing genetic networks. Biostatistics. 14, 3:
586-599.

Fentaw Abegaz and Ernst Wit (2014). SparseTSCGM: Sparse time series chain graphical
models. R package version 2.1.1. http://CRAN.R-project.org/package=SparseTSCGM

Rothman, A.J., Levina, E., and Zhu, J. (2010). Sparse multivariate
regression with covariance estimation. Journal of Computational and
Graphical Statistics. 19: 947-962.

Wild, B., Eichler, M., Friederich, H. C., Hartmann, M., Zipfel, S., &
Herzog, W. (2010). A graphical vector autoregressive modelling
approach to the analysis of electronic diary data. BMC medical
research methodology, 10(1), 28.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

\examples{
# Simulate model:
Mod <- randomGVARmodel(4,probKappaEdge = 0.8,probBetaEdge = 0.8)

# Simulate data:
Data <- graphicalVARsim(100,Mod$beta,Mod$kappa)

# Estimate model:
Res <- graphicalVAR(Data, gamma = 0, nLambda = 10)

# Plot results:
layout(t(1:2))
plot(Mod, "PCC", layout = "circle")
plot(Res, "PCC", layout = "circle")

plot(Mod, "PDC", layout = "circle")
plot(Res, "PDC", layout = "circle")
}
