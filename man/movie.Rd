\name{movie}
\alias{movie}
\docType{data}
\title{Movie reviews from a Variety magazine}
\description{
  This dataset contains reviews of 93 movies by four critics
  in the Variety magazine from April 1995 through March 1997.}
\usage{data(movie)}
\format{
  This data frame contains the following columns:
  \describe{
    \item{movie}{
      Identification number of a movie.
    }
    \item{y}{
      review response, 1: Positive;
      2: Mixture of positive and negative; 3: Negative.
    }
    \item{critic}{
      an ordered factor: critics with levels
      \code{medved} < \code{siskel} < \code{ebert} < \code{lyons}
    }
  }
}
\source{
  Hartzel, J., Agresti, A., Caffo, B. (2001).
  Multinomial logit random effects models. \emph{Statistical Modelling},
  1: 81-102
}
}
\references{
Jokinen J, McDonald JW, Smith PWF. Meaningful regression and
association models for clustered ordinal data.
\emph{Sociological Methodology} 2006; 36:173-199.
}
\keyword{datasets}
