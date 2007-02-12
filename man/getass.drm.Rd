\name{getass.drm}
\alias{getass.drm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Support function for drm: specification of the association model}
\description{
A support function called by \code{\link{drm}} that parses from the argument
\code{dep} the covariates and functional forms for the association
model. This function should not be used directly by the user.
}
%- maybe also 'usage' for other objects documented here.

\details{
In the argument \code{dep}, the user needs to specify the covariates and
functional forms on the association parameters. The following list
describes the notation and concise interpretation of the association
parameters. For more details, see the reference below.

\item{\code{nu}}{For "N"-structure: the proportion of susceptibles in the
  population, i.e. proportion of units that can have a value greater
  than the smallest category. Covariates can be specified. See examples
  below.}

\item{\code{nu1, nu2, ..., nuk}}{For "L"-structure: the proportion in the population
  in the latent class 1, 2, ..., k. The proportion in the latent class 0
  can be obtained by calculating 1-\code{nu1}-\code{nu2}-... Note that
  for binary responses,  the number of latent classes can be specified
  with an argument \code{Lclass} (Default 2). Covariates can be
  specified. See examples below.}

\item{\code{kappa0, kappa1,...,kappak-1} (binary response)}{For "L"-structure:
  the ratio of success probabilities in the latent class 0, 1, ..., k-1
  compared to the success probability in the latent class k. Covariates
  can be specified. See examples below.}

\item{\code{kappa0, kappa1,...,kappak-1} (multicategorical response)}{For
  "L"-structure with two latent classes: the ratio of category
  probabilities in categories 0, 1,..., k-1 in the latent class 0
  compared to the corresponding probabilities in the latent class 1.}

\item{\code{xi1, xi0} (binary response)}{For "B"-structure: the shape
  parameters of the Beta-distribution (sometimes also noted as xi1=p and
  xi0=q). Covariates can be specified. See examples below.}

\item{\code{xi0, xi1, ..., xik} (multicategorical response)}{For "D"-structure:
  the shape parameters of the Dirichlet distribution}

\item{\code{tau} (binary responses)}{For "M"-structure: adjacent second order
  dependence ratio. If the number of repeated measurements is greater
  than two, the adjacent \code{tau}'s are assumed to be equal
  (i.e. stationarity of the dependence ratios). In order to specify
  equalities or functional forms (i.e. non-stationary overlapping
  dependence ratios), see examples below.}

\item{\code{tau12, tau13, tau123} (binary responses)}{For "M2"-structure:
  adjacent second order dependence ratio (\code{tau12}), adjacent third
  order dependence ratio (\code{tau123}) and the second order dependence
  ratio between first and third response (\code{tau13}). If the number
  of repeated measurements is greater than three, the \code{tau}'s are
  assumed to be equal (i.e. stationarity of the dependence
  ratios). Equalities and functional forms can be specified. See
  examples below.}

\item{\code{tau11, tau12, ..., tau21, tau22, ...,taukk} (multicategorical
  responses)}{For "M"-structure: adjacent second order
  dependence ratios for categories 1, 2, ..., k, where the smallest
  response value is regarded as the baseline, denoted by 0. Equalities
  and functional forms can be specified as for the binary "M"-structures
  above. See examples below.} 
}

\references{
  Jokinen J. Fast estimation algorithm for likelihood-based
  analysis of repeated categorical responses.
  \emph{Computational Statistics and Data Analysis} 2006; 51:1509-1522.}


\seealso{\code{\link{drm}}}
\examples{
### Example of functional forms:
## non-stationary second order Markov structure
## initial values of the dependence ratios are set to 1.
data(wheeze)
assoc <- list("M2",
              tau12 ~ function(a78=1, a89=1, a910=1)c(a78, a89, a910),
              tau123 ~ function(a789=1, a8910=1)c(a789, a8910),
              tau13 ~ function(a79=1, a810=1)c(a79, a810))

fit1 <- drm(wheeze~I(age>9)+smoking+cluster(id)+Time(age),
            data=wheeze, dep=assoc, print=0)

### Example of other parameter restrictions:
## fixing parameters to a known value: ~tau12==1, ~tau21==1
## setting parameters to equal: ~tau11==tau22
data(marijuana)
assoc <- list("M", ~tau12==1, ~tau21==1, ~tau11==tau22)

fit2 <- drm(y~age+cluster(id)+Time(age), data=marijuana,
            subset=sex=="female", dep=assoc, print=0)

## setting all parameters to equal:

assoc <- list("M", ~tau11==tau12, ~tau11==tau21, ~tau11==tau22)

fit3 <- drm(y~age+cluster(id)+Time(age), data=marijuana,
            subset=sex=="female", dep=assoc, print=0)

### Example of covariates for the association parameters:
## allow the latent class prevalence and the corresp. conditional
## probabilities vary by sex. Note: covariate needs to be a factor.

assoc <- list("L",
              nu1 ~ nu1:factor(sex),
              kappa0 ~ kappa0:factor(sex),
              kappa1 ~ kappa1:factor(sex))

fit4 <- drm(y~age+cluster(id)+Time(age), data=marijuana,
            dep=assoc, print=0)

### Example how to derive conditional probabilities from marginals
## Fit a model with three latent classes:

latent3 <- drm(wheeze~I(age>9)+smoking+cluster(id),data=wheeze,
               dep="L",Lclass=3, print=0)

## calculate conditional probabilities:
## pr(Y=1|L=2) = pr(Y=1)/(nu2+kappa1*nu1+kappa0*(1-nu1-nu2))
## pr(Y=1|L=1) = kappa1*pr(Y=1|L=2)
## pr(Y=1|L=0) = kappa0*pr(Y=1|L=2)

est <- coef(latent3)
psi2 <- latent3$fitted.marginals/
        (est["nu2"]+est["kappa1"]*est["nu1"]+
         est["kappa0"]*(1-est["nu1"]-est["nu2"]))
psi1 <- psi2*est["kappa1"]
psi0 <- psi2*est["kappa0"]

## check the model validity, i.e. require that 0 < psi_i <1:
range(cbind(psi0,psi1,psi2))

}
\keyword{utilities} % at least one, from doc/KEYWORDS
