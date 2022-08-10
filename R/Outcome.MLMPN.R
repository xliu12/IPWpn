#’ @title
#’
#’ @description
#’
#’ @param  outcome
#’ @param  Treatment_assignment
#’ @param  cluster_assignment
#’ @param  general_confounders
#’ @param  clustering_confounders
#’ @param  data
#’
#’ @return A list.
#’ @import MplusAutomation, tidyverse
#’ @examples
#’ rnorm(10)
#’
#’ @export Outcome.MLMPN


options(dplyr.print_max = 1e9)
options(tibble.width = Inf)

########### adjust for the raw score of W and the raw score of C ################
Outcome.MLMPN <- function(
    outcome ='Y',
    Treatment_assignment='Trt',
    cluster_assignment='K',
    general_confounders,
    clustering_confounders,
    data
){ # the conventional

  Y=data[ , c(outcome)]; Trt=data[ , c(Treatment_assignment)];
  clus=data[ , c( cluster_assignment)];
  Lyz=data[ ,c( general_confounders, clustering_confounders)]
  Ly=data[ , c(general_confounders)]

  Y = drop(Y)
  clus = drop(clus)
  Trt = drop(Trt)

  Ly = as.matrix((Ly));if(ncol(Ly)>0){ colnames(Ly)=paste0('Ly', 1:ncol(Ly)) }
  Lyz = as.matrix(Lyz); if(ncol(Ly)>0){ colnames(Lyz)[1:ncol(Ly)]=paste0('Ly', 1:ncol(Ly)) }
  if( (ncol(Lyz)>ncol(Ly)) & (ncol(Ly)>0) ){ colnames(Lyz)=c(paste0('Ly', 1:ncol(Ly)), paste0('Lz', 1:(ncol(Lyz)-ncol(Ly)) )) }
  if( (ncol(Lyz)>ncol(Ly)) & (ncol(Ly)==0) ){ colnames(Lyz)=c(  paste0('Lz', 1:(ncol(Lyz)-ncol(Ly)) )) }

  dat=data.frame(Y=Y,Trt=Trt,clus=clus, Lyz)

  dat$clus = factor(dat$clus) # Declare ic as a factor

  dat1=dat[dat$Trt==1,]
  dat0=dat[dat$Trt==0,]

  if( (ncol(Lyz)>ncol(Ly)) &  (ncol(Ly)<=0) ){
    treat.outcome.mlm=paste0('Y ~ 1+', paste0('Lz', 1:(ncol(Lyz)-ncol(Ly)), collapse  = '+'), '+ ( 1 | clus)' )
  }
  if( (ncol(Lyz)<=ncol(Ly)) &  ( ncol(Ly)>0) ) {
    treat.outcome.mlm=paste0('Y ~ 1+', paste0('Ly', 1:ncol(Ly), collapse  = '+'), '+ ( 1 | clus)' )
  }
  if( (ncol(Lyz)>ncol(Ly)) &  (ncol(Ly)>0) ){
    treat.outcome.mlm=paste0('Y ~ 1+', paste0('Ly', 1:ncol(Ly), collapse  = '+'), '+',paste0('Lz', 1:(ncol(Lyz)-ncol(Ly)), collapse  = '+'), '+ ( 1 | clus)' )
  }
  treatfit = lmer(as.formula( treat.outcome.mlm ), data = dat1 )
  is_singular=isSingular(treatfit)
  # treatfit = summary(treatfit)
  gamm1 = fixef(treatfit)[1]
  # gamm1_var = diag(vcov(treatfit))[1] # or summary(treatfit)$coef[1,2]
  gamm1_var = summary(treatfit)$coef[1,2]^2
  gamm1_se = sqrt(gamm1_var)

  varcomp = as.data.frame(VarCorr(treatfit))
  varu1 = varcomp[1,c("vcov")]
  vare1 = varcomp[varcomp$grp=="Residual",c("vcov")]

  # estimates for the control arm T==0

  # controlfit = lm(Y ~ 1 +  Ly, data=dat0)
  if( ncol(Ly) ==0){
    control.outcome.reg= paste0('Y ~ 1')
  }
  if( ncol(Ly) >0){
    control.outcome.reg= paste0('Y ~ 1+', paste0('Ly', 1:ncol(Ly), collapse  = '+') )
  }
  controlfit = lm( as.formula(control.outcome.reg ), data=dat0)
  # controlfit = summary(controlfit)
  gamm0 = coefficients(controlfit)[1] #controlfit$coefficients[1,1]
  gamm0_var = diag(vcov(controlfit))[1]
  gamm0_se = sqrt(gamm0_var)
  vare0 = sigma(controlfit)^2

  # Combine the estimates from the treatment and control groups # and perform the test
  gammDiff = gamm1 - gamm0
  se_sw = sqrt( gamm1_var+gamm0_var )
  se_sw = as.numeric(se_sw)
  z.wald = gammDiff/se_sw

  # CIup=gammDiff + qnorm(0.975)*se_sw
  # CIlow=gammDiff - qnorm(0.975)*se_sw

  out=list(
    gammDiff, se_sw, z.wald
    # gamm1, varu1, vare1, gamm0, vare0, se_mplus, ACov
  )
  names(out)=as.character(expression(
    Estimate, SE, z.wald
    # gamm1, varu1, vare1, gamm0, vare0, se_mplus, ACov
  ))
  outcome.mlmpn=out
  return(out)
}
