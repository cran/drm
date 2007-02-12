\name{marijuana}
\alias{marijuana}
\docType{data}
\title{Longitudinal data set on teenage marijuana use}
\description{
  This data set is from US National Youth Survey reporting
  a repeated questionnaire on marijuana use of 237 teenagers
  yearly from ages 13 to 17 on an ordinal scale.
}
\usage{data(marijuana)}
\format{
This data frame contains the following columns:
\describe{
\item{id}{
identification number of a teenager.
}
\item{y}{
  marijuana response, 1: never;
  2: less than once a month; 3: more than once a month.
}
\item{sex}{
a factor with levels
\code{female} 
\code{male} 
}
\item{age}{
age as a factor with levels
\code{13} 
\code{14} 
\code{15} 
\code{16} 
\code{17} 
}
}
}
\source{
  Lang, JB, McDonald, JW, Smith, PWF. (1999). Association
  marginal modelling of multivariate categorical responses: A maximum
  likelihood approach. \emph{JASA}, 94, 1161-1171.  
}
\references{
Jokinen J, McDonald JW, Smith PWF. Meaningful regression and
association models for clustered ordinal data.
\emph{Sociological Methodology} 2006; 36:173-199.
}
     
\keyword{datasets}
