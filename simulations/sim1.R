#######################################################################################################
## CSIR ESTIMATION SIMULATIONS WITH NO X->T RELATIONSHIP, NO X->Y RELATIONSHIP, NO T->Y RELATIONSHIP ##
## created by Rachel Nethery                                                                         ##
#######################################################################################################

nsims<-5000

library(MatchIt)
library(mvtnorm)
source('../functions/cdc_sir.R')
source('../functions/est_poisson.R')

## read in the confounder dataset used to generate the data ##
Xan<-read.csv('sim_confounders.csv',header=T,stringsAsFactors=F)

names(Xan)<-c('MoneyFood','MoneySmoke','P65Plus','PMale','PWhite','Unemploy','Commute','Income','Pop')

## remove the population variable from the dataset ##
popan<-Xan$Pop

Xan<-Xan[,-ncol(Xan)]

## define constants ##
N<-nrow(Xan)
p<-ncol(Xan)

## set seed for reproducibility ##
set.seed(4)

## generate the exposure ##
E<-rbinom(n=N,size=1,prob=.24)

## sample 10 'exposed' CBGs to represent the community of interest ##
community<-sample(which(E==1),size=10)

## do matching ##
## subset to only our exposed community of interest and all unexposed communities ##
inds_matchin<-c(community,which(E==0))
Ematchin<-E[inds_matchin]
Xmatchin<-Xan[inds_matchin,]
matchin<-data.frame(Ematchin,Xmatchin)
row.names(matchin)<-1:nrow(matchin)

## apply 20:1 mahalanobis distance nearest neighbor matching procedure ##
matchout<-matchit(Ematchin~MoneyFood+
                    MoneySmoke+
                    P65Plus+PMale+PWhite+Unemploy+
                    Commute+Income,
                  data=matchin,method='nearest',distance='mahalanobis',ratio=20)
## see if the matching has achieved balance in the exposure groups ##
summary(matchout,standardize=T)

## output the matched dataset ##
all_matched<-match.data(matchout)
Xmatchout<-all_matched[,-c(ncol(all_matched)-1,ncol(all_matched))]
names(Xmatchout)[1]<-'E'

popanmatchin<-popan[inds_matchin]
popanmatchout<-popanmatchin[as.numeric(row.names(all_matched))]

## store results ##
cdc_results<-matrix(NA,nrow=nsims,ncol=3)
pr_results<-matrix(NA,nrow=nsims,ncol=3)
csir_results<-matrix(NA,nrow=nsims,ncol=3)

## run sims ##
for (i in 1:nsims){
  ## generate the outcomes ##
  Y<-rpois(n=N,lambda=0.0025*popan)
  Ymatchin<-Y[inds_matchin]
  Ymatchout<-Ymatchin[as.numeric(row.names(Xmatchout))]
  
  ## a. compute SIR and 95% CI using the CDC's method ##
  Yseer<-Y[-community]
  POPseer<-popan[-community]
  Ycom<-Y[community]
  POPcom<-popan[community]
  cdc_results[i,]<-cdc_sir(Ycom=Ycom,POPcom=POPcom,Yseer=Yseer,POPseer=POPseer)
  
  ## b. compute SIR and 95% CI using a poisson regression without matching ##
  Epr<-rep(0,N)
  Epr[community]<-1
  pr<-glm(Y~as.matrix(cbind(Epr,Xan)),offset=log(popan))
  pr_results[i,]<-exp(c(coef(pr)[2],confint(pr)[2,]))
  
  ## c. compute cSIR and 95% CI ##
  csir_temp<-est_poisson(X=as.matrix(Xmatchout),Y_all=Ymatchout,offs=popanmatchout,nsamp=30000,tau=.5)
  csir_results[i,]<-csir_temp[['ebeta_postsum']][2,]
}

save(csir_results,pr_results,cdc_results,file='sim1.RData')
