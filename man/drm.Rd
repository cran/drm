\name{drm}
\alias{drm}
\alias{print.drm}
\alias{summary.drm}
\alias{print.summary.drm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Combined regression and association models for clustered categorical responses}
\description{
  \code{drm} fits a combined regression and association model for longitudinal
  or otherwise clustered categorical responses using dependence ratio as 
  a measure of the association.}
\usage{
drm(formula, family=binomial, data=sys.parent(), weights, offset,
subset=NULL, na.action, start=NULL, link="cum", dep="I", Ncond=TRUE,
Lclass=2, dropout=FALSE, drop.x=NULL, save.profiles=TRUE, pmatrix=NULL,
print.level=2, iterlim=200, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{formula}{a formula expression as for other regression models.
       In  addition  the  cluster term has to be specified in the
       expression by \code{cluster()} and if using temporal  association
       structure the temporal term has to be specified by \code{Time()}.
       See examples below and the documentation of \code{lm} and \code{formula}
       for further details.}
  \item{family}{a description of the link function to be  used  in  the
    model for a binary response. Default is logit link. 
    See \code{family} for details. For an ordinal response, link
    is defined for the cumulative probabilities when
    \code{link}-argument is set to "cum". See \code{link} below.}
  \item{data}{an optional data frame containing the variables in  the
       model.}
  \item{weights}{an optional vector of  weights  to  be  used  in  the
       fitting  process.  Only  equal  weights within cluster are
       allowed.}
  \item{offset}{this can be used to specify an a priori known component
    to be included in the linear predictor during fitting.}
  \item{subset}{an optional vector specifying a subset of
       observations to be used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when
       the  data  contain  NAs.  The  default is \code{na.include} after
       which the analysis assumes missing data mechanism at random (MAR) if
       \code{dropout=FALSE}, and not at random if
       \code{dropout = TRUE}. See \code{dropout} below.}
  \item{start}{an optional vector of starting values for the
       parameters. If not specified by the user, the starting values
       are estimated from \code{glm}-procedure assuming independence}
 \item{link}{this can be used to specify alternative link functions
       for nominal and ordinal responses. By default "cum", after which
       the link is specified through \code{family = binomial(link=?)} for the
       cumulative probabilities. Alternative links include adjacent
       category logit "acl" and baseline category logit "bcl" (baseline
       category being the last category). For "bcl", the regression
       parameters are estimated for each logit level. For a binary response,
       this argument is ignored.}
 \item{dep}{\code{dep} defines the association structure. The default is
       independence  "I".  Other  singular  options  are for  the
       exchangeable association: Necessary  factor  "N",  Latent
       categorical factor "L", Latent Beta-distributed propensity "B"
       (binary response), Latent Dirichlet-distributed propensities "D"
       (multicategorical response), and for the
       temporal association: first order  Markov  "M", and second order
       Markov "M2" (binary response). By default,
       Markov structure for the adjacent 2-way dependence ratios is
       assumed to be stationary. Superpositions  of  these  structures  
       can be imposed,  such  as  "NL", "NB","ND", "NM", "LM", "NLM","NM2".
       See [3-7] for further details.
       Parameter restrictions, covariates and  functional  forms  for  the
       association parameters can also be specified. In that case
       the \code{dep}-argument must be a list. See examples below. For
       the interpretation of the association parameters, see the
       documentation of the support function \code{\link{getass.drm}}.}
  \item{Ncond}{logical argument defining whether the regression model is marginal or
    conditional when the association is "N". The default is \code{TRUE}, 
    i.e. the regression estimation is conditional on {N=1}. If
    covariates are used for the "N"-association, it is advisable to set
    \code{Ncond=FALSE}, since otherwise the interpretation of the
    regression parameters is less clear.}
  \item{Lclass}{Number of latent classes in the population when the association is
    "L". Default is 2. Available only for binary response. Note that in the
  current implementation, the conditional probabilities are not
  calculated for \code{Lclass}>2. For checking the validity of the
  model, the user needs to check whether the estimated conditional
  probabilities fall within 0 and 1. See example in \code{\link{getass.drm}} for
  parameter interpretation and how to calculate the conditional probabilities}
  \item{dropout}{logical argument. For monotone missing patterns in longitudinal studies,
    this argument allows to impose a selection model (see [8] for
    details) on top of regression and association model to investigate the sensitivity of
    the results due to missingness. The model formula notation is:
    \code{logit(hz(drop.cur)) = (Intercept)d+response.cur+response.prev} ,
    where \code{response.cur} denotes the effect of current, possibly missing
    response value and \code{response.prev} denotes the effect of previous
    response value. MCAR, MAR and MNAR-models can be
    specified by imposing restrictions on selection model parameters
    in \code{dep}-argument as for the association parameters. See
    \code{dep} above and examples below. If the response is a factor, the effect of
    the factor levels are estimated contrasting to the lowest level.}    
  \item{drop.x}{an optional covariate vector for the selection model.
    The covariate's previous value (notation: \code{covariate.prev}) is used in
    the selection model.}
  \item{save.profiles}{logical argument defining whether the fitted
    values of all possible profiles is saved. If \code{FALSE}, only the
    indicator vector (-1 for a negative, 1 for a positive profile) over all
    units will be saved. If the cluster size is large, using
    \code{save.profiles=TRUE} may result to a very large object.}
  \item{pmatrix}{a character object specifying the name of the matrix
    for all possible profiles, created using
    \code{\link{profiles.drm}}. If the cluster size is large,
    this speeds up the estimation in case several models are fitted.
    See examples below.}
  \item{print.level}{level of printing during  numerical  optimisation.
        The default is 2. See \code{\link{nlm}} for further details.}
  \item{iterlim}{maximum iteration limit for the numerical maxisimisation.
    See \code{\link{nlm}} for further details.} 
  \item{\dots}{other arguments passed to \code{nlm}, e.g. controlling the
       convergence. See \code{\link{nlm}} for further details.}
}
\details{\code{drm} gives maximum likelihood estimates for the combined
       regression  and association model by decomposing a
       joint probability of responses in a cluster to univariate 
       marginal or cumulative probabilities and dependence ratios of all
       orders. See [1] and [5] for further details. The dimensionality
       of the association part is reduced by imposing a model for
       the  association  structure with \code{dep}-argument. See
       \code{\link{getass.drm}} and [3-7] for details. Furthermore, a
       selection model can be added on top of regression and association
       model. See examples below and [5] and [8] for details.}

    \value{\code{drm} returns an object of class \code{drm}.
      The function \code{summary} (i.e., \code{summary.drm}) can be
      used to obtain or print a summary of the results. The generic
      accessor function \code{coefficients} can be used to extract
      coefficients.

      An object of class \code{drm} is a list containing at least the
      following components:

      \item{coefficients}{a named vector of regression, and possibly
	association and selection model coefficients.}

      \item{cov.scaled}{a variance-covariance matrix of the parameter
	estimates.}

      \item{fitted.marginals}{the fitted values for the univariate
      means, obtained by transforming the linear predictors by the
      inverse of the link function.}

    \item{fitted.conditionals}{in case of "L"-structure, the fitted
      values for the conditional univariate means, otherwise NULL. Not
      yet implemented for \code{Lclass}>2; see also
      \code{\link{getass.drm}}} 

    \item{fitted.profiles}{the fitted response profile probabilities
      within each cluster, calculated by using the maximum likelihood
      estimates from the model. See also \code{save.profiles} above}

    \item{deviance}{minus twice the maximised log-likelihood.}

     \item{aic}{Akaike's Information Criterion, minus twice the maximised
       log-likelihood plus twice the number of coefficients. Not
       available if the likelihood is weighted with the dropout
       probabilities.} 

     \item{niter}{the number of iterations that \code{nlm} used.}

     \item{code}{convergence code from \code{nlm}. See \code{nlm} for details.}

     \item{call}{the matched call.}

     \item{terms}{the `terms' object used.}
}
\references{1. Ekholm A, Smith PWF, McDonald  JW.  Marginal  regression
       analysis  of  a  multivariate  binary response. \emph{Biometrika}
       1995; 82(4):847-854.

       2. Ekholm A, Skinner C. The Muscatine children's obesity data
       reanalysed using pattern mixture models. \emph{Applied
	 Statistics} 1998; 47:251-263.

       3. Ekholm A, McDonald JW, Smith PWF.  Association  models
       for  a  multivariate  binary  response. \emph{Biometrics} 2000;
       56:712-718.

       4. Ekholm A, Jokinen J, Kilpi T. Combining regression and
       association modelling on longitudinal data on bacterial carriage.
       \emph{Statistics in Medicine} 2002; 21:773-791.
       
       5. Ekholm A, Jokinen J, McDonald JW, Smith PWF. Joint regression
       and association modelling of longitudinal ordinal data.
       \emph{Biometrics} 2003; 59:795-803.

       6. Jokinen J, McDonald JW, Smith PWF. Meaningful regression and
       association models for clustered ordinal data.
       \emph{Sociological Methodology} 2006; 36:173-199.

       7. Jokinen J. Fast estimation algorithm for likelihood-based
       analysis of repeated categorical responses.
       \emph{Computational Statistics and Data Analysis} 2006; 51:1509-1522.

       8. Diggle PJ, Kenward MJ. Informative dropout in longitudinal
       data analysis. \emph{Applied Statistics} 1994;  43: 49-94.}

\author{Jukka Jokinen, \email{jukka.jokinen@helsinki.fi}}
\section{WARNING}{The maximum likelihood estimates may sometimes lead
  to negative fitted probabilities. In this case, both generic
  print-methods warn about this. In this case, the model is considered
  to be wrongly specified and model specification should be changed.}

\seealso{ \code{\link{getass.drm}}, \code{\link{nlm}},
\code{\link{cluster}}, \code{\link{Time}} 
  \code{\link{profiles.drm}}, \code{\link{depratio}}}

\examples{
######################################################
## Examples for binary responses
###########################################
## Wheeze among Steubenville (see [3]):
## Latent Beta-distributed propensity
data(wheeze)
fit1 <- drm(wheeze~I(age>9)+smoking+cluster(id),data=wheeze,dep="B", print=0)

## Obesity among Muscatine children (see [2]):
## Analysis for completers: M2 for girls
data(obese)
fit2 <- drm(obese~age+cluster(id)+Time(year), subset=sex=="female",
            dep="M2",data=obese)

\dontrun{
## Muscatine children continued (see [3]):
## LM for boys and girls separately
fit3 <- drm(obese~age+cluster(id)+Time(age), subset=sex=="male",
            dep="LM",data=obese)

fit4 <- drm(obese~age+cluster(id)+Time(age), subset=sex=="female",
            dep="LM",data=obese)
}
############################################
## Examples for ordinal responses
############################################
## Movie critic example (see [6]):
## Latent Dirichlet propensities with baseline category link.
data(movie)

options(contrasts=c("contr.treatment","contr.treatment"))
fit5 <- drm(y~critic+cluster(movie), data=movie, dep="D", link="bcl")

## Longitudinal dataset on teenage marijuana use (see [6]):
## Superposition of structures N, L and M for the girls.
data(marijuana)

fit6 <- drm(y~age+cluster(id)+Time(age), data=marijuana,
            subset=sex=="female", dep=list("NLM", ~kappa1==1,
            ~kappa2==0, ~tau12==1, ~tau21==1, ~tau11==tau22))

## Parameter restrictions with functions using M-structure for the boys.
## Plot the second order dependence ratios:
plot(depratio(y~cluster(id)+Time(age), data=marijuana,
     subset=sex=="male"))

## fit the model in [6]:
fit7 <- drm(y~age+cluster(id)+Time(age), data=marijuana,
            subset=sex=="male", dep=list("M", 
            tau12~function(a=1,b=0) a+b*c(0:3),
            tau21~function(a=1,b=0) a+b*c(0:3)))

\dontrun{
##############################################
## Covariates for the association (see [7]):
##############################################
data(madras)

## plot empirical 2nd order dependence ratios with bootstrap CI's
tau.madras <- depratio(symptom~cluster(id)+Time(month), data=madras,
                       boot.ci = TRUE, n.boot = 1000)
plot(tau.madras, log="y", ylim=c(1,40), plot.ci=TRUE)

## create matrix for profiles:
W.madras <- profiles.drm(n.categories=2, n.repetitions=12, "M")

## create four-level covariate, combining age and sex:
madras$age.sex <- factor(paste(madras$age,madras$sex,sep="."))

## fit the model in [7], Section 4:
fit8 <- drm(symptom~age+sex+month+month:age+month:sex+cluster(id)+Time(month),
            data=madras, Ncond=FALSE, save.profiles=FALSE, pmatrix="W.madras",
            dep=list("NM",nu~nu:age.sex,
                     tau~function(a0=0,a1=0) 1+a0*exp(a1*c(0:10))), print=2)

###################################################
## Dropout model on top of regression & association 
###################################################
## Continue with the madras data.
## fit a model without the dropout model:
fit9 <- drm(symptom~age+sex+month+month:age+month:sex+cluster(id)+Time(month),
            data=madras, save.profiles=FALSE, pmatrix="W.madras", print=0,
            dep=list("NM", tau~function(a0=0,a1=0) 1+a0*exp(a1*c(0:10))))

## A dropout model assuming MCAR for the thought disorders:

mcar <- drm(symptom~age+sex+month+month:age+month:sex+cluster(id)+Time(month),
            data=madras, save.profiles=FALSE, pmatrix="W.madras",
            dep=list("NM", tau~function(a0=0,a1=0) 1+a0*exp(a1*c(0:10)),
                     ~symptom.cur==0,~symptom.prev==0),
            dropout=TRUE, start=c(coef(fit9), -4))

## A dropout model assuming MAR; including sex as a covariate:

mar <- drm(symptom~age+sex+month+month:age+month:sex+cluster(id)+Time(month),
           data=madras, save.profiles=FALSE, pmatrix="W.madras",
           dep=list("NM", tau~function(a0=0,a1=0) 1+a0*exp(a1*c(0:10)),
                    ~symptom.cur==0), dropout=TRUE, drop.x=sex,
           start=c(coef(mcar),0,0))

## A dropout model assuming MNAR and sex as a covariate:

mnar <- drm(symptom~age+sex+month+month:age+month:sex+cluster(id)+Time(month),
            data=madras, save.profiles=FALSE, pmatrix="W.madras",
            dep=list("NM", tau~function(a0=0,a1=0) 1+a0*exp(a1*c(0:10))),
            dropout=TRUE, drop.x=sex, start=c(coef(mcar),0,0,0))

## print out coefficients and std.errors:
coef(summary(mnar))
}
## std.error of `symptom.cur' all over the place; too few dropouts
## for a comprehensive evaluation of the dropout mechanism
}
\keyword{models}







