\name{obese}
\alias{obese}
\docType{data}
\title{Obesity among Muscatine children}
\description{
Data of the completers concerning the presence or abscence of obesity in five age
cohorts in Muscatine, Ohio, measured repeatedly in -77, -79 and -81. 
}
\usage{data(obese)}
\format{
  A data frame with 5310 observations on the following 5 variables.
  \describe{
    \item{id}{identification number of a child}
    \item{obese}{a binary indicator of obesity}
    \item{sex}{a factor with levels \code{male} \code{female}}
    \item{year}{year of measurement}
    \item{age}{a factor indicating age cohorts from \code{(5,7]} to \code{(17,19]}}
  }
}

\source{
  Woolson RF and Clarke WR. Analysis of categorical incomplete
  longitudinal data.
  \emph{Journal of the Royal Statistical Society Series A},
  1984; 147:87-99.

http://www.biostat.harvard.edu/~fitzmaur/ala/obesity.txt
}
\references{
Ekholm A, Skinner C. The Muscatine children's obesity data
reanalysed using pattern mixture models. \emph{Applied
  Statistics} 1998; 47:251-263.

Ekholm A, McDonald JW, Smith PWF.  Association  models
for  a  multivariate  binary  response. \emph{Biometrics} 2000;
56:712-718.}

\keyword{datasets}
