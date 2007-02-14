\name{cluster}
\alias{cluster}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Identify clusters}
\description{
  This is a special function used in the \code{drm} or \code{depratio}
  syntax. 
     It identifies correlated groups of observations, and is used on
     the right hand side of a formula.
}
\usage{
cluster(x)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{A character, factor, or numeric variable.}
}
\details{
The function's only action is semantic, to mark a variable as the
cluster indicator.
}
\value{x}
}
\seealso{\code{\link{drm}}, \code{\link{depratio}}, \code{\link{Time}}}

\examples{
formula(wheeze~cluster(id))
}
\keyword{misc}
