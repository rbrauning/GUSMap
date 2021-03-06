% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rfEst.R
\name{$computeMap}
\alias{$computeMap}
\alias{rf_est_FS}
\title{FS method: Compute ordered linkage maps}
\usage{
FSobj$computeMap(chrom=NULL, init_r=0.01, ep=0.001, method="optim", sexSpec=FALSE, 
                    err=TRUE, mapped=TRUE, nThreads=1)
}
\arguments{
\item{init_r}{A numeric value giving the initial values for the recombination fractions. Each 
recombination fraction parameter is set to the same initial value.}

\item{ep}{A numeric value giving the initial value for the sequencing error parameter.}

\item{sexSpec}{Logical value. If \code{TRUE}, sex-specific recombination fractions are
are estimated.}

\item{method}{A character string specifying whether optimization should be performed using
direct maximization (\code{optim}) or via the Expectation-Maximum (EM) algorithm (\code{EM}).}

\item{nThreads}{An integer value giving the number of threads to use in computing the likelihood in parallel.}

\item{chrom}{A integer vector giving the indices of the chromosomes (or linkage groups) to be computed.}

\item{err}{Locical value. If \code{TRUE}, the sequencing error parameter is estimated. Otherwise, the
sequenicng error parameter is fixed to the value of the \code{ep} argument.}

\item{mapped}{Locial value. If \code{TRUE}, the maps are computed using the marker order given 
in the combined linkage groups. Otherwise, the maps are computed using the original marker order  
given by the genomic assembly.}
}
\description{
Method for inferring parental phase (e.g., ordered parental genotype pair (OPGP)) and
estimating recombination fractions in full-sib families.
}
\details{
This function infers the parental phase (or ordered parental genotype pair (OPGP)) and
estimates adjacent recombination fractions using the hidden Markov model (HMM) approach as
described in \insertCite{bilton2018genetics1;textual}{GUSMap}.

The optimization of the likelihood for the HMM is performed using either the Expectation-Maximumization (EM) algorithm
(\code{method="EM"}) or using direct numeric optimization via the \code{\link{optim}} function (\code{method="optim"}).
The likelihood computations (and computation of derivatives if required) are scaled using 
forward and backward recursion to avoid overflow issues and are performed in C. These computations 
are also parallelization via the OpenMP package, where the argument \code{nThreads} specifies
how many threads to use. Be careful not to set \code{nThreads} to more than the number of threads available
on your computer (or bad things will happen). In addition, if the package is complied without OpenMP, then this 
parallelization has no effect and the likelihood is computed in serial.

If \code{mapped = TRUE}, then combined linkage groups must have been formed from the \code{\link{$addBIsnps}} function
first (and preferably ordered from the \code{\link{$orderLG}} function).
}
\examples{
#### Case 1: Compute linkage map from linkage groups
## Simulate some sequencing data
set.seed(6745)
config <- list(list(sample(c(1,2,4), size=30, replace=TRUE)))
F1data <- simFS(0.01, config=config, meanDepth=10, nInd=50)
## Compute 2-point recombination fractions
F1data$rf_2pt(nClust=1)
## create and order linkage groups
F1data$createLG()
F1data$addBIsnps()
F1data$orderLG(ndim=5)

## Compute the linkage map
F1data$computeMap()

#### Case 2: Compute map using original assembly order
F1data$computeMap(mapped = FALSE)
}
\references{
\insertRef{bilton2018genetics1}{GUSMap}
}
\seealso{
\code{\link{FS}}
}
\author{
Timothy P. Bilton and Chris Scott
}
