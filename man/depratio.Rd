\name{depratio}
\alias{depratio}
\alias{print.depratio}
\alias{plot.depratio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Empirical estimates of the dependence ratios}
\description{
  Calculates the observed values of the adjacent dependence ratios
  from the data.
}
\usage{
depratio(formula, data, subset, ord = 2, boot.ci = FALSE, n.boot = NULL,
         ci.width=0.95) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{the syntax is of form \code{y~cluster(id)+Time(time)},
    where \code{id} denotes the cluster indicator, and \code{Time}
    denotes the order along which the adjacent dependence ratios will be
  calculated.}
  \item{data}{optional data frame containing the variables in the formula}
  \item{subset}{an optional vector specifying a subset of observations
    from the data}
  \item{ord}{order of the dependence ratios to be calculated. The
    default is 2}
  \item{boot.ci}{logical argument specifying whether bootstrap
    confidence intervals will be calculated for the empirical dependence
    ratio estimates}
  \item{n.boot}{number of bootstrap replicates}
  \item{ci.width}{width of the confidence interval. Default is 0.95}
}
\value{An object of class \code{depratio}. Generic functions
  \code{print} and \code{plot} are also available.

  An object of class \code{depratio} is a list containing at least the
  following two components:
  \item{tau}{matrix of the observed dependence ratios}
  \item{freq}{matrix of the frequencies of events for the numerator of
  the observed dependence ratios}
}

\seealso{
  \code{\link{drm}}, \code{\link{cluster}}, \code{\link{Time}}
  }

\examples{
## calculate and plot the observed 2nd order dependence ratios
## for the marijuana data:
data(marijuana)
dr.male <- depratio(y~cluster(id)+Time(age), data=marijuana,
                    subset=sex=="male")
dr.male
plot(dr.male)

## confirm that the 1st order Markov assumption is adequate
## for the madras data:
data(madras)

dr2 <- depratio(symptom~cluster(id)+Time(month), data=madras)
dr3 <- depratio(symptom~cluster(id)+Time(month), ord=3, data=madras)
dr <- rbind(dr2$tau[-length(dr2$tau)]*dr2$tau[-1], dr3$tau)

matplot(1:ncol(dr), t(dr))

}
\keyword{arith}% at least one, from doc/KEYWORDS
