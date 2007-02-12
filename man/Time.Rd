\name{Time}
\alias{Time}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Identify ordering of the responses}
\description{
  This is a special function used in the \code{drm} or \code{depratio}
  syntax. It identifies the ordering in correlated groups of observations,
and is used on the right hand side of a formula.
}
\usage{
Time(x)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{A factor, or numeric variable.}
}
\details{
The function's only action is semantic, to mark a variable as the
Time indicator.
}
\seealso{\code{\link{drm}}, \code{\link{depratio}}, \code{\link{cluster}}}

\examples{
data(wheeze)
drm(wheeze~cluster(id)+Time(age),data=wheeze,dep="M", print=0)
}

\keyword{misc}%-- one or more ...
