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
## use madras-data
data(madras)

## create matrix for profiles (takes a couple of minutes):
W.madras <- profiles.drm(2, 12, "M")

## create four-level covariate, combining age and sex:
madras$age.sex <- factor(paste(madras$age,madras$sex,sep="."))

## fit the model using W.madras (takes less than a couple of minutes):
fit <- drm(symptom~age+sex+month+month:age+month:sex+cluster(id)+Time(month),
           data=madras, Ncond=FALSE, save.profiles=FALSE, pmatrix="W.madras",
           dep=list("NM",nu1~nu1:age.sex,
                    tau~function(a0=0,a1=0) 1+a0*exp(a1*c(0:10))), print=1)

## a speedy way of saving the profile probabilities to the object:
fit <- drm(symptom~age+sex+month+month:age+month:sex+cluster(id)+Time(month),
           data=madras, Ncond=FALSE, save.profiles=TRUE, pmatrix="W.madras",
           dep=list("NM",nu1~nu1:age.sex,
                    tau~function(a0=0,a1=0) 1+a0*exp(a1*c(0:10))),
           print=0, start=coef(fit))

}
\keyword{array}% at least one, from doc/KEYWORDS
