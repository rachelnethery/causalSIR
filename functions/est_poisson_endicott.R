est_poisson<-function(X,Y_all,offs,C,nsamp=ncol(as.matrix(Y_all)),pburn=.25,tau=.001){
  
  ##########################
  #### DEFINE CONSTANTS ####
  ##########################
  
  ## values used in proposal dist'n ##
  Cinv<-solve(C)
  
  Y_all<-as.matrix(Y_all)
  X<-cbind(1,X)
  nburn<-ceiling(nsamp*pburn)
  p<-ncol(X)
  Nm<-nrow(Y_all)
  
  #######################################
  #### SET STARTING AND PRIOR VALUES ####
  #######################################
  
  ## priors ##
  beta0<-rep(0,p)
  sigma0<-rep(.5,p)
  
  ## starting values ##
  beta<-matrix(0,nrow=p,ncol=1)
  accept<-0
  
  ## proposal variance ##
  propvar<-(tau*diag(p))%*%C%*%(tau*diag(p))
  
  #################
  #### STORAGE ####
  #################
  beta_save<-matrix(NA,nrow=nsamp-nburn,ncol=p)
  
  ##################################
  #### START METROPOLIS SAMPLER ####
  ##################################
  
  for (g in 1:nsamp){
    
    ## define Y ##
    if (ncol(as.matrix(Y_all))>1){
      Y<-Y_all[,g,drop=F]
    } else{
      Y<-Y_all
    }
    
    ## propose a new beta vector from MVN, call it beta_star ##
    beta_star<-matrix(rmvnorm(n=1,mean=beta,sigma=propvar),ncol=1)
    
    ## compute the acceptance ratio ##
    log_ar<-(sum(dnorm(beta_star,mean=beta0,sd=sqrt(sigma0),log=T))+sum(dpois(Y,exp(X%*%beta_star+as.matrix(log(offs))),log=T)))-
      (sum(dnorm(beta,mean=beta0,sd=sqrt(sigma0),log=T))+sum(dpois(Y,exp(X%*%beta+as.matrix(log(offs))),log=T)))
      #t(beta_star)%*%solve(sigma0)%*%beta_star)-(2*t(beta_star)%*%solve(sigma0)%*%beta0)+sum((X%*%beta_star)*Y-exp(X%*%beta_star))-
      #((t(beta)%*%solve(sigma0)%*%beta)-(2*t(beta)%*%solve(sigma0)%*%beta0)+sum((X%*%beta)*Y-exp(X%*%beta)))
    ar<-exp(log_ar)
    
    ## decide whether to accept beta_star ##
    u<-runif(n=1,min=0,max=1)
    if (u<ar){
      beta<-beta_star
      if (g>nburn) accept<-accept+1
    }
    
    ## if past burn-in, save the new beta value ##
    if (g>nburn){
      beta_save[g-nburn,]<-c(beta)
    }
  }
  
  return(list('beta_postsamp'=beta_save,
              'acceptance'=accept/(nsamp-nburn),
              'ebeta_postsum'=cbind(apply(exp(beta_save),2,mean),
              apply(exp(beta_save),2,quantile,.025),
              apply(exp(beta_save),2,quantile,.975)
              )))
  
}
