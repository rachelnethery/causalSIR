pred_multinom<-function(X,Y,cID,offs,nsamp=10000,nburn=5000,tau=.75){
  
  ##############################################
  #### DATA MANAGEMENT AND DEFINE CONSTANTS ####
  ##############################################
  
  ## values used in proposal dist'n ##
  C<-vcov(glm(Y~X,family='poisson'))[-1,-1]
    
  Y<-matrix(Y[order(cID)],ncol=1)
  X<-X[order(cID),]
  
  p<-ncol(X)
  N<-length(unique(cID))
  n_i<-c(table(cID))
  Ntot<-nrow(Y)
  
  ###########################################
  ###### SET STARTING AND PRIOR VALUES ######
  ###########################################
  
  ## priors ##
  mu0<-matrix(0,nrow=p,ncol=1)
  sigma0<-diag(p)
  
  ## starting values ##
  beta<-matrix(0,nrow=p,ncol=1)
  accept<-0
  
  ## proposal variance ##
  propvar<-(tau*diag(p))%*%C%*%(tau*diag(p))
    
  #################
  #### STORAGE ####
  #################
  beta_save<-matrix(NA,nrow=nsamp,ncol=p)
  
  ##################################
  #### START METROPOLIS SAMPLER ####
  ##################################
  
  for (g in 1:(nburn+nsamp)){
    
    ## propose a new beta vector from MVN, call it beta_star ##
    beta_star<-matrix(rmvnorm(n=1,mean=beta,sigma=propvar),ncol=1)
    
    ## set up a matrix of the normalization factor ##
    nf_beta<-matrix(rep(tapply(c(exp(X%*%beta)*offs),cID,sum),times=n_i),ncol=1)
    nf_beta_star<-matrix(rep(tapply(c(exp(X%*%beta_star)*offs),cID,sum),times=n_i),ncol=1)
    
    ## compute the acceptance ratio ##
    log_ar<-(sum(Y*(X%*%beta_star+as.matrix(log(offs))-log(nf_beta_star)))+t(beta_star)%*%solve(sigma0)%*%beta_star-2*t(mu0)%*%solve(sigma0)%*%beta_star)-
      (sum(Y*(X%*%beta+as.matrix(log(offs))-log(nf_beta)))+t(beta)%*%solve(sigma0)%*%beta-2*t(mu0)%*%solve(sigma0)%*%beta)
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
  
  return(list(beta_save,accept/nsamp))
  
}