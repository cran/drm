\name{profiles.drm}
\alias{profiles.drm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creates a support matrix for `drm' for profile probability calculation}
\description{
  Creates a support matrix for the calculation of all possible
  profile probabilities using the ML estimates obtained from
  \code{\link{drm}}.
}
\usage{
profiles.drm(n.categories, n.repetitions, structure = "exchangeable")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n.categories}{Number of respose variable categories}
  \item{n.repetitions}{Number of response repetitions within a cluster}
  \item{structure}{Association structure indicator of the model to be
    fitted: "exchangeable" if the structure is exchangeable, otherwise
    either "M" or "M2"}
}
\details{Creation of this support object is useful if the cluster size
  is large and if several models are fitted to the same dataset. Avoids
  creating the same matrix over and over again for each model fit.
  See examples for \code{\link{drm}}.
}
\value{If \code{structure}="exchangeable", returns a
  matrix with row and column size equal to n.categories to the power of
  n.repetitions. For temporal structures "M" and "M2", the result is an
  array. 
}
\seealso{\code{\link{drm}}}
\examples{
## use marijuana-data
data(marijuana)

## create matrix for profiles:
W.mari <- profiles.drm(3, 5, "M")

## fit the model using W.mari:
fit <- drm(y~age+cluster(id)+Time(age), data=marijuana,
           subset=sex=="male",pmatrix="W.mari", save.profiles=FALSE,
           dep=list("M",tau12~function(a=1,b=0) a+b*c(0:3),
                    tau21~function(a=1,b=0) a+b*c(0:3)),print=0)

## a speedy way of saving the profile probabilities to the object:
fit <- drm(y~age+cluster(id)+Time(age), data=marijuana,
           subset=sex=="male",pmatrix="W.mari", save.profiles=TRUE,
           dep=list("M",tau12~function(a=1,b=0) a+b*c(0:3),
                    tau21~function(a=1,b=0) a+b*c(0:3)),print=0,
           start=coef(fit))

}
\keyword{array}% at least one, from doc/KEYWORDS
