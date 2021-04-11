\name{IPW.meandifference}
\alias{IPW.meandifference}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimating the treatment effect from a partially nested design using the inverse-propensity-weighted mean difference estimation approach
%%  ~~function to do ... ~~
}
\description{
IPW.meandifference is used to estimate the treatment effect from a partially nested design (where the treatment arm has multiple treatment clusters, and  the control arm has no clusters) with the inverse-propensity-weighted mean difference estimation approach adapted to the partially nested design feature.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
IPW.meandifference(
  Y,
  Trt,
  clus,
  Lyz,
  Ly
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Y}{An outcome variable}
\item{Trt}{Treatment assignment indicator (1 for treatment and 0 for control). The treatment arm has multiple treatment clusters, and  the control arm has no clusters.
}
\item{clus}{Observed treatment cluster assignment. clus = 0 for the control arm}
\item{Lyz}{A matrix containing pre-treatment covariates. Lyz can contain two types of covariates. The first type of covariates affect both the treatment assignment and outcome directly. The second type of covariates do not affect the outcome in a given treatment cluster directly, but affect both the treatment assignment and treatment cluster assignment.}
\item{Ly}{A matrix containing the pre-treatment covariates that affect both the treatment assignment and outcome directly. Ly is the first type of covariates contained in Lyz.}
}

\details{
The sandwisch type standard error estimation does not account for the clustering in the treatment arm.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
IPW.meandifference returns a list "ipw10.md" containing the following components:
%%  If it is a LIST, use list(meanDiff,se_sw, z.wald, mu1, mu0, ACov)
%%  \item{comp1 }{Description of 'comp1'}
\item{meanDiff }{the treatment effect estimate, i.e., the mean difference between the treatment and control arm. }
\item{se_sw }{the sandwich-type standard error estimate of the treatment effect estimate. }
\item{z.wald }{the Wald statistic (i.e., meanDiff/se_sw). }
}
%%  \references{
%% ~put references to the literature/web site here ~ }
%%  \author{
%%  ~~who you are~~
%%  }
%%  \note{
%%  ~~further notes~~ }

%% ~Make other sections like Warning with \section{Warning }{....} ~
%%  \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~}




\examples{
%\dontrun{
\donttest{
%##-- ==>  Define data, use random,
data(dat_obs)
head(dat_obs)
IPW.meandifference(Y=dat_obs$Y,
                   Trt=dat_obs$Trt,
                   clus=dat_obs$clus,
                   Lyz=dat_obs[,c('Ly','Lz')],
                   Ly=dat_obs[,c('Ly')] )

%}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~design }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~models }% __ONLY ONE__ keyword per line
