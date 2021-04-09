#’ @title
#’
#’ @description
#’
#’ @param  Y
#’ @param  Trt
#’ @param  clus
#’ @param  Lyz
#’ @param  Ly
#’
#’ @return A list.
#’ @import MplusAutomation, tidyverse
#’ @examples
#’ rnorm(10)
#’
#’ @export IPW.MLMPN

suppressMessages(library(tidyverse))
library(MplusAutomation)

options(dplyr.print_max = 1e9)
options(tibble.width = Inf)


IPW.MLMPN <- function(
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

# IPW MLM-PN ----

# IPW.MLM_PN_sw <- function(
#   Y,
#   Trt,
#   clus,
#   ps1,
#   ps0,
#   Lyz,
#   Ly
#   ){

  sandwich.se=TRUE

  pid=Sys.getpid()

  Lyz = as.matrix(Lyz)
  Ly = as.matrix((Ly))

  Y = drop(Y)
  clus = drop(clus)
  Trt = drop(Trt)
  ps1 = drop(ps1)
  ps0 = drop(ps0)

  w = Trt/ps1 + (1-Trt)/(1-ps0)
  dat1=data.frame(Y,Trt,clus,Lyz,w)[Trt==1,]
  dat0=data.frame(Y,Trt,clus,Ly,w)[Trt==0,]

  # estimates for the treatment arm T==1
  # datfname = 'dat1.txt'
  # write.table(dat1, file = file.path(tempdir(), datfname)
  #             ,quote = F,row.names = F,col.names = F,na="999")
  datfname=paste(pid,"dat1.txt",sep = "")
  write.table(dat1, datfname,quote = F,row.names = F,col.names = F,na="999")

  # mplus using w as sampling weights
  cat(sep = " ",
      "
      Data:
      file =", datfname,";
      VARIANCES=NOCHECK;

      Variable:
      names = ", colnames(dat1), ";
      usevariables = Y;
      cluster = clus;
      weight = w;
      wtscale = unscaled;

      Analysis:
      type = twolevel;
      estimator = MLR;

      Model:
      %WITHIN%
      Y;
      %BETWEEN%
      Y;
      [Y] (gamma1);

      Output:
      NOCHISQUARE TECH1 TECH2 TECH3 TECH4;

      ",
      # file=file.path(tempdir(),"wmlm.inp")
      file = paste(sep="",pid,"wmlm.inp")
  )

  runModels(paste(sep="",pid,"wmlm.inp"), logFile = NULL)
  mplus_out=readModels(paste(sep="",pid,"wmlm.out"), what = c("all"))

  file.remove(paste(sep="",pid,"wmlm.inp"))
  file.remove(paste(sep="",pid,"wmlm.out"))
  file.remove(paste(pid,"dat1.txt",sep = ""))


  if(length(mplus_out$errors)==0){
    mplus_param=mplus_out$parameters$unstandardized[c(2,3,1),]
    theta=mplus_param$est
    if(length(theta)!=3){
      theta=rep(NA,3)
    }
    names(theta)=c("gamm1","varu1","vare1")
    se_mplus=mplus_param$se
    if( length(se_mplus) != 3 ){
      se_mplus=rep(NA,3)
    }
    names(se_mplus)=paste("se_",names(theta),sep = "")
    gamm1=theta[1]
    varu1=theta[2]
    vare1=theta[3]
  }
  if(length(mplus_out$errors)>0){
    theta=rep(NA,3)
    names(theta)=c("gamm1","varu1","vare1")
    se_mplus=rep(NA,3)
    names(se_mplus)=paste("se_",names(theta),sep = "")
    gamm1=theta[1]
    varu1=theta[2]
    vare1=theta[3]
  }

  # estimates for the control arm T==0
  gamm0 = sum(Y*(1-Trt)/(1-ps0)) / sum((1-Trt)/(1-ps0))
  vare0 = sum( (Y-gamm0)^2*(1-Trt)/(1-ps0) ) / sum( (1-Trt)/(1-ps0) )

  if(sandwich.se){
    # sandwich-type standard error estimate accounting for uncertainty in ps
    S1i = as.matrix((Trt-ps1)*Lyz)
    S0i = as.matrix((Trt-ps0)*Ly)

    # Q1
    wz1i = tibble(Y, Trt, clus, w) %>%
      filter(Trt==1) %>%
      group_by(clus) %>%
      mutate( wDtotal=w*(Y-gamm1) ) %>%
      mutate( wSwithin=w*( Y-sum(w*Y)/sum(w) )^2 )

    wz1.sum = wz1i %>%
      filter(Trt==1) %>%
      summarise_all( ~sum(.) ) %>%
      rename_at(vars(-clus), ~paste("z.sum_",., sep = "") ) %>%
      mutate(
        Q1gamm1.z1 = z.sum_wDtotal / (vare1+z.sum_w*varu1)
      ) %>%
      mutate(
        Q1varu1.z1 = -z.sum_w/(2*(vare1+z.sum_w*varu1)) + (Q1gamm1.z1)^2/2
      ) %>%
      mutate(
        Q1vare1.z1 = -(z.sum_w-1)/(2*vare1) + (z.sum_wSwithin)/(2*vare1^2) + Q1varu1.z1/z.sum_w
      )

    Q1z1 = select_at(wz1.sum, vars(contains("Q1")) )

    # A33

    a33 = wz1.sum %>%
      mutate(
        A33.11 = -z.sum_w/(vare1+z.sum_w*varu1)
      )  %>%
      mutate(
        A33.12 = A33.11*Q1gamm1.z1
      ) %>%
      mutate(
        A33.13 = A33.12/z.sum_w
      ) %>%
      mutate(
        A33.21 = A33.12
      ) %>%
      mutate(
        A33.22 = A33.11^2/2 + (Q1gamm1.z1)^2*A33.11
      ) %>%
      mutate(
        A33.23 = A33.22/z.sum_w
      ) %>%
      mutate(
        A33.31 = A33.13
      ) %>%
      mutate(
        A33.32 = A33.23
      ) %>%
      mutate(
        A33.33 = (z.sum_w-1)/(2*vare1^2) - z.sum_wSwithin/(vare1^3) + A33.22/(z.sum_w^2)
      )

    A33=matrix(as.numeric( summarise_at(a33, vars(contains("A33")), sum) ),3,3)

    # A31
    A31=matrix(0,3,ncol(as.matrix(Lyz)))

    for(alph in 1:ncol(Lyz)){
      oz1.sum = bind_cols(wz1i, tibble(L=Lyz[Trt==1,alph]) ) %>%
        filter(Trt==1) %>%
        group_by(clus) %>%
        mutate( dw.dalph=w*(1-1/w)*L ) %>%
        mutate( odds=dw.dalph) %>%
        mutate( oDtotal=odds*(Y-gamm1) ) %>%
        mutate( oSwithin=odds*(Y-Y-sum(w*Y)/sum(w))^2 ) %>%
        summarise_all( ~sum(.) ) %>%
        rename_at(vars(-clus), ~paste("z.sum_",., sep = "") )

      z1.sum =
        bind_cols(Q1z1, oz1.sum) %>%
        # bind_cols(wz1.sum, oz1.sum) %>%
        mutate(
          a31.1 = -z.sum_oDtotal/(vare1+z.sum_w*varu1) + z.sum_odds*z.sum_wDtotal*varu1/(vare1+z.sum_w*varu1)^2
        ) %>%
        mutate(
          a31.2 = -z.sum_odds*vare1/2/(vare1+z.sum_w*varu1)^2 + Q1gamm1.z1*a31.1
        ) %>%
        mutate(
          a31.3 = -z.sum_odds/(2*vare1) + z.sum_oSwithin/(2*vare1^2) - z.sum_odds*Q1varu1.z1/(z.sum_w^2) + a31.2/z.sum_w
        )

      A31[,alph] = as.numeric( summarise_at(z1.sum, vars(contains("a31")), sum) )
    }


    # control arm T==0
    # Q0
    Q0i = tibble(Y,Trt,w) %>% filter(Trt==0) %>%
      mutate(
        Q0gamm0i=(Y-gamm0)*w/vare0
      ) %>%
      mutate(
        Q0vare0i=(Y-gamm0)^2*w/(2*vare0^2)-w/(2*vare0)
      ) %>%
      select_at( vars(contains("Q0")) )

    # A44
    A44i = tibble(Y,Trt,w) %>% filter(Trt==0) %>%
      mutate(
        A44.11 = -w/vare0
      ) %>%
      mutate(
        A44.12 = -(Y-gamm0)*w/vare0^2
      ) %>%
      mutate(
        A44.21 = A44.12
      ) %>%
      mutate(
        A44.22 = w/(2*vare0^2) -(Y-gamm0)^2*w/vare0^3
      )

    A44 = matrix(as.numeric(summarise_at(A44i, vars(contains("A44")),sum)),2,2)


    # A42
    A42 = matrix(0, 2, ncol(as.matrix(Ly)))
    for( alph in 1:ncol(as.matrix(Ly)) ){
      A42i = tibble(Y,Trt,w,L=as.matrix(Ly)[,alph]) %>%
        filter(Trt==0) %>%
        mutate(
          odds = w*(1-1/w)*L
        ) %>%
        mutate(
          A42.1= -odds*(Y-gamm0)/vare0
        ) %>%
        mutate(
          A42.2= -odds/(2*vare0) + odds*(Y-gamm0)^2/(2*vare0^2)
        )

      A42[,alph] = as.numeric(summarise_at(A42i, vars(contains("A42")), sum) )
    }


    # independent units: treatment clusters and control individuals (each as a singleton cluster)
    S = cbind(S1i, S0i)
    Sz = aggregate(S, by=list(clus=clus), sum)
    Sz1 = Sz[Sz$clus!=0,-1]
    Sz0 = S[Trt==0,]
    colnames(Sz0) = colnames(Sz1)
    Sj = rbind(Sz1,Sz0)

    Q1j= rbind( as.matrix(Q1z1), matrix(0, nrow(Q0i), ncol(Q1z1)) )
    Q0j= rbind( matrix(0, nrow(Q1z1),ncol(Q0i)), as.matrix(Q0i) )

    Uj=as.matrix(cbind(Sj, Q1j, Q0j))
    B=t(Uj)%*%Uj

    A11 = t(Lyz)%*%diag(ps1*(1-ps1))%*%Lyz
    A22 = t(Ly)%*%diag(ps0*(1-ps0))%*%Ly

    A = rbind(
      cbind(A11, matrix(0,nrow(A11),ncol(A22)), matrix(0,nrow(A11),ncol(A33)), matrix(0,nrow(A11),ncol(A44)) )
      , cbind(matrix(0,nrow(A22),ncol(A11)), A22, matrix(0,nrow(A22),ncol(A33)), matrix(0,nrow(A22),ncol(A44)) )
      , cbind(A31, matrix(0,nrow(A31),ncol(A22)), A33, matrix(0,nrow(A31),ncol(A44)) )
      , cbind(matrix(0,nrow(A42),ncol(A11)), A42, matrix(0,nrow(A42),ncol(A33)), A44 )
    )

    gammDiff = gamm1-gamm0
    tryAinv = try(solve(A))

    if( ( is.na(A[1,1]) | ( nrow(tryAinv)!=nrow(A) ) ) ){
      ACov=NA
      se_sw = NA
      z.wald = NA
    }
    if( ( (!is.na(A[1,1])) & ( nrow(tryAinv)==nrow(A) ) ) ){
      ACov = solve(A)%*%B%*%solve(A)
      gamm.contrast = c(rep(0,ncol(as.matrix(Lyz))+ncol(as.matrix(Ly))), 1,0,0,-1,0)
      se_sw = sqrt( t(gamm.contrast)%*%ACov%*%gamm.contrast )
      se_sw = as.numeric(se_sw)
      z.wald = gammDiff/se_sw
    }
  }

  if(!sandwich.se){
    gammDiff = gamm1-gamm0
    se_sw = sqrt( (se_mplus[1])^2+vare0/sum( (1-Trt)/(1-ps0) ) )
    se_sw = as.numeric(se_sw)
    z.wald = gammDiff/se_sw
    ACov = NA
  }

  out=list(
    gammDiff, se_sw, z.wald
    #, gamm1, varu1, vare1, gamm0, vare0, se_mplus, ACov
  )
  names(out)=as.character(expression(
    gammDiff, se_sw, z.wald
    #, gamm1, varu1, vare1, gamm0, vare0, se_mplus, ACov
  ))


#   return(out)
#
# }



# IPW.MLMPN <- function(
#   Y,
#   Trt,
#   clus,
#   Lyz,
#   Ly
# ){
# the proposed ----
  # sandwich.se=TRUE
  # dat_ps = EstPS(Trt = Trt,
  #                Lyz = Lyz,
  #                Ly = Ly)
  # sandwich.se = sandwich.se
  #
  # out10=IPW.MLM_PN(
  #   Y = Y,
  #   Trt = Trt,
  #   clus = clus,
  #   ps1 = dat_ps[,"ps1"],
  #   ps0 = dat_ps[,"ps0"],
  #   Lyz = Lyz,
  #   Ly = Ly,
  #   sandwich.se = sandwich.se
  # )
  ipw10.mlmpn=out
  return(ipw10.mlmpn)
}

