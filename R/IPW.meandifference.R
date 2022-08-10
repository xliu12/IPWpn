
suppressMessages(library(tidyverse))
library(MplusAutomation)

options(dplyr.print_max = 1e9)
options(tibble.width = Inf)

IPW.meandifference <- function(
  Y,
  Trt,
  clus,
  Lyz,
  Ly
){ # the proposed

# PS estimation ----
# EstPS <- function(
#   # dat,
#   Trt, # 1=Treatment,0=Control
#   Lyz,
#   Ly
# ){
#
  Lyz = as.matrix(Lyz)
  Ly = as.matrix((Ly))
  glm1 = glm(Trt~Lyz-1, family = "binomial")
  glm0 = glm(Trt~Ly-1, family = "binomial")
  ps1=glm1$fitted.values
  ps0=glm0$fitted.values

#   dat_ps = cbind(Trt, ps1, ps0,Lyz,Ly)
#
#   return(dat_ps)
# }
#
#
# # IPW MeanDiff ----
#
# IPW.MeanDiff <- function(
#   Y,
#   Trt,
#   ps1,
#   ps0,
#   Lyz,
#   Ly
#
# ){

  sandwich.se=TRUE
  n=length(Y)

  Y = drop(Y)
  Trt = drop(Trt)
  ps1 = drop(ps1)
  ps0 = drop(ps0)

  # proposed: IPW10_MeanDiff
  mu1 = sum(Y*Trt/ps1)/sum(Trt/ps1)
  mu0 = sum(Y*(1-Trt)/(1-ps0))/sum((1-Trt)/(1-ps0))
  meanDiff = mu1-mu0

  if(sandwich.se){
    # sandwich-type standard error estimate accounting for uncertainty in ps
    S1i = (Trt-ps1)*Lyz
    S0i = (Trt-ps0)*Ly
    H1i = (Y-mu1)*Trt/ps1
    H0i = (Y-mu0)*c(1-Trt)/(1-ps0)
    Ui=cbind(S1i, S0i, H1i, H0i)
    B=t(Ui)%*%Ui

    A11 = t(Lyz)%*%diag(ps1*(1-ps1))%*%Lyz
    A22 = t(Ly)%*%diag(ps0*(1-ps0))%*%Ly
    A31 = t(H1i)%*%diag((1-ps1))%*%Lyz
    A33 = t(Trt)%*%(1/ps1)
    A42 = t(H0i)%*%diag(ps0)%*%Ly
    A44 = t(1-Trt)%*%(1/(1-ps0))

    A = rbind(
      cbind(A11, matrix(0,nrow(A11),ncol(A22)), matrix(0,nrow(A11),ncol(A33)), matrix(0,nrow(A11),ncol(A44)) )
      , cbind(matrix(0,nrow(A22),ncol(A11)), A22, matrix(0,nrow(A22),ncol(A33)), matrix(0,nrow(A22),ncol(A44)) )
      , cbind(A31, matrix(0,nrow(A31),ncol(A22)), A33, matrix(0,nrow(A31),ncol(A44)) )
      , cbind(matrix(0,nrow(A42),ncol(A11)), A42, matrix(0,nrow(A42),ncol(A33)), A44 )
    )

    tryAinv = try(solve(A))

    if( ( is.na(A[1,1]) | ( nrow(tryAinv)!=nrow(A) ) ) ){
      ACov=NA
      se_sw = NA
      z.wald = NA
    }
    if( ( (!is.na(A[1,1])) & ( nrow(tryAinv)==nrow(A) ) )  ){
      ACov = solve(A)%*%B%*%solve(A)
      mu.contrast = c(rep(0,ncol(as.matrix(Lyz))+ncol(as.matrix(Ly))), 1,-1)
      se_sw = as.numeric(sqrt( t(mu.contrast)%*%ACov%*%mu.contrast ))
      z.wald = meanDiff/se_sw
    }

  }

  if(!sandwich.se){
    se_sw = sqrt( sum( (Y-mu1)^2*Trt/ps1 )/(sum( Trt/ps1 ))^2 +
                    sum( (Y-mu0)^2*(1-Trt)/(1-ps0) )/(sum( (1-Trt)/(1-ps0) ))^2 )
    se_sw = as.numeric(se_sw)
    ACov=NA
    z.wald = meanDiff/se_sw
  }
  out=list(
    meanDiff,se_sw, z.wald
  )
  names(out)=as.character(expression(
    meanDiff,se_sw, z.wald
  ))
#
#   return(out)
# }

# IPW.meandifference <- function(
#   Y,
#   Trt,
#   clus,
#   Lyz,
#   Ly
#   ){
# the proposed ----
  sandwich.se=TRUE
  # dat_ps = EstPS(Trt = Trt,
  #                Lyz = Lyz,
  #                Ly = Ly)
  # sandwich.se = sandwich.se
  #
  # ipw10.md=IPW.MeanDiff(
  #   Y = Y,
  #   Trt = Trt,
  #   ps1 = dat_ps[,"ps1"],
  #   ps0 = dat_ps[,"ps0"],
  #   Lyz = Lyz,
  #   Ly = Ly,
  #   sandwich.se = sandwich.se
  # )
  ipw10.md=out
  return(ipw10.md)
}



