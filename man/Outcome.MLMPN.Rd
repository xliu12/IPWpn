\name{Outcome.MLMPN}
\alias{Outcome.MLMPN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating the treatment effect from a partially nested design with the multilevel modeling based outcome models
}
\description{
Outcome.MLMPN is used to estimate the treatment effect from a partially nested design (where the treatment arm has multiple treatment clusters, and  the control arm has no clusters) with the  multilevel modeling based outcome models (where a multilevel model is specified for the treatment outcome and a regression model is specified for the control outcome).
}
\usage{
Outcome.MLMPN(Y,
  Trt,
  clus,
  general_confounders,
  clustering_confounders)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{outcome}{a string indicating the name of the column containing the outcome variable.}
\item{Treatment_assignment}{a string indicating the name of the column containing the treatment assignments (1 for treatment and 0 for control). (In the partially nested design, the treatment group has multiple clusters, and  the control group has no cluster.)
%%     ~~Describe \code{x} here~~
}
\item{cluster_assignment}{a string indicating the name of the column containing the observed cluster assignments. For an individual in the treatment group, the cluster assignment is the cluster in which the individual is assigned to receive the treatment (e.g., an integer in 1, 2, ,..., J, where J is the total number of clusters in the treatment group). For an individual in the control group, set the cluster assignment as 0.}
%\item{Lyz}{A matrix containing pre-treatment covariates. Lyz can contain two types of covariates. The first type of covariates affect both the treatment assignment and outcome directly. The second type of covariates do not affect the outcome in a given treatment cluster directly, but affect both the treatment assignment and treatment cluster assignment.}
\item{clustering_confounders}{A character vector containing the names of the clustering confounders. The clustering confounders are those baseline confounders that directly affect both the treatment assignment  and cluster assignment in the treatment group, but do not directly influence the outcome  once the treatment assignment and cluster assignment are given.}
\item{general_confounders}{A character vector containing the names of the general confounders. The general confounders are those baseline confounders that directly affect both the treatment assignment and outcome within a cluster of the treatment group. (When one is uncertain about whether a covariate is a general confounder or a clustering confounder, the covariate should be classified as a general confounder to be conservative, because a general confounder can directly affect the outcome, cluster assignment, and treatment assignment.)
}
\item{data}{A data frame containing all variables.}
}
}
\details{
The multilevel model for the treatment outcome includes a random cluster-specific intercept. The covariates do not have random cluster-specific slopes. 
}
\value{
Outcome.MLMPN returns a list containing the following components:
%%  If it is a LIST, use list(gammDiff,se_sw, z.wald)
%%  \item{comp1 }{Description of 'comp1'}
\item{Estimate }{the average treatment effect estimate. }
\item{SE }{the standard error estimate of the average treatment effect estimate. }
\item{z.wald }{the Wald statistic (i.e., Estimate/SE). 
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
%\dontrun{
\donttest{
%##-- ==>  Define data, use random,
data(dat_obs)
head(dat_obs)

res= Outcome.MLMPN( outcome = 'Y', Treatment_assignment = 'Trt', cluster_assignment = 'K',
                general_confounders = c('C'), clustering_confounders = c('W') , data = dat_obs)
res

%}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ ~design }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~models }% __ONLY ONE__ keyword per line
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
