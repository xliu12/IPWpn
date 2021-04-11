\name{dat_obs}
\alias{dat_obs}
\docType{data}
\title{
Observed data in a partially nested design.
}
\description{
A dataset simulated from a partially nested design containing two types of pre-treatment confounders, treatment assignments,
treatment cluster assignments, and outcome scores.
}
\usage{data("dat_obs")}
\format{
  A data frame with 200 observations on the following 5 variables.
  \describe{
    \item{\code{Trt}}{a numeric vector}
    \item{\code{clus}}{a numeric vector}
    \item{\code{Y}}{a numeric vector}
    \item{\code{Ly}}{a numeric vector}
    \item{\code{Lz}}{a numeric vector}
  }
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
%%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%%}
%%\references{
%%  ~~ possibly secondary sources and usages ~~
%%}
\examples{
data(dat_obs)
}
\keyword{datasets}
