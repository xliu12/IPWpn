#' A simulated observed dataset from a partially nested design.
#'
#' A dataset simulated from a partially nested design
#' containing two types of pre-treatment confounders, treatment assignments,
#' treatment cluster assignments, and outcome scores.
#'
#'
#' @format A data frame with 200 rows and 5 variables:
#' \describe{
#'   \item{Trt}{Treatment assignment indicator (1 for treatment and 0 for control).}
#'   \item{K}{Observed cluster assignment. clus = 0 for the control arm.}
#'   \item{Y}{An outcome variable}
#'   \item{C}{The general confounder that affects both the treatment assignment and outcome directly.}
#'   \item{W}{The clustering confounder that does not affect the outcome in a given cluster directly, but affects both the treatment assignment and cluster assignment in the treatment group.}
#' }
#' @source \url{https://github.com/xliu12/IPWpn}
"dat_obs"
