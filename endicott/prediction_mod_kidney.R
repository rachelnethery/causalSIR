library(mvtnorm)
library(coda)

source('../functions/pred_multinom_offset.R')

## read in NY prediction model data ##
preddat<-read.csv('PredModelData.csv',header=T,stringsAsFactors=F)

## construct set of predictors to be used in analysis ##
Xpred<-as.matrix(preddat[,c(28,31:35,38:42)])

## construct county ID variable ##
cID<-floor(preddat$geoid10/10000000)

fit1<-pred_multinom(X=Xpred,Y=preddat$observed_Kidney,cID=cID,offs=preddat$Pop,nsamp=200000,nburn=10000)

save(fit1,file='kidney_output.RData')

pdf('kidney_traceplots.pdf',width=12)
par(mfrow=c(3,4))
for (i in 1:11){
  traceplot(mcmc(fit1[[1]][,i])) 
}
dev.off()

print(paste0('Acceptance Rate=',fit1[[2]]))

print(data.frame('var'=colnames(Xpred),'est'=apply(fit1[[1]],2,mean,na.rm=T),
           'll'=apply(fit1[[1]],2,quantile,probs=.025),
           'ul'=apply(fit1[[1]],2,quantile,probs=.975)))

