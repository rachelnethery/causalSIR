###############################################################################################
## CODE TO PERFORM ENDICOTT SENSITIVITY ANALYSIS OMITTING MONEYFOOD, MONEYSMOKE, EXERCISE    ##
## CREATED BY RACHEL NETHERY                                                                 ##
###############################################################################################

#####################
## data management ##
#####################

library(MatchIt)
library(mvtnorm)
library(coda)
source('../functions/est_poisson_endicott.R')

## read in the new york cancer and confounder data ##
ny2010<-read.csv('NYdata.csv',stringsAsFactors=F)
ny2010$county<-floor(ny2010$geoid10/10000000)

## endicott cbgs ##
eecbg<-c(360070131003,360070137001,360070135001,360070137002,360070136003,'36007DOH0011','36007DOH0010')

## endicott dataset ##
eedat<-ny2010[which((ny2010$dohregion %in% eecbg)==T),]

## read in the CBG confounder and exposure data (for seer+) ##
XTseer<-read.csv('XTseer.csv',header=T,stringsAsFactors=F)

##############
## matching ##
##############

## find unexposed "controls" ##
controldat<-XTseer[which(XTseer$TRIexpose==0 & XTseer$SPFexpose==0),]

## remove rural areas from potential controls, as every CBG in endicott area is non-rural ##
controldat<-controldat[which(controldat$Rural==0),]

## remove the county that endicott is in ##
controldat<-controldat[-which(floor(controldat$ID/10000000)==36007),]

## make dataset for input into the matching function ##
matchin<-data.frame(rbind(eedat[,c(33:35,38:41)],
                          controldat[,c(7:9,12:15)]))

matchin$E<-c(rep(1,nrow(eedat)),rep(0,nrow(controldat)))
rownames(matchin)<-1:nrow(matchin)

for (M in c(3,5)){
  for (metho in c("mahalanobis","logit","GAMlogit")){
    
    ## run matching procedure ##
    matchout<-matchit(E~P65Plus+PMale+PWhite+Unemploy+
                        Commute+Income+Industry,
                      data=matchin,method='nearest',ratio=M,distance=metho,replace=F)
    
    ## output balance info ##
    temptab<-cbind(summary(matchout,standardize=T)[[3]][,4],summary(matchout,standardize=T)[[4]][,4])
    rownames(temptab)<-rownames(summary(matchout,standardize=T)[[3]])
    colnames(temptab)<-c("Before","After")
    write.csv(temptab,file=paste0('balance_',M,'_',metho,'_sa2.csv'))
    
    ## construct matched dataset ##
    matchdat<-match.data(matchout)
    matchdat<-matchdat[,-c(ncol(matchdat)-1,ncol(matchdat))]
    controlind<-as.numeric(row.names(matchdat))[(nrow(eedat)+1):nrow(matchdat)]-nrow(eedat)
    
    ## ids of control units ##
    controlid_all<-controldat$ID[controlind]
    
    ## separate control units within NY (where we havee CBG data) and outside NY (where we only have county data) ##
    controlid_ny<-controlid_all[which(floor(controlid_all/10000000000)==36)]
    controlid<-controlid_all[-which(floor(controlid_all/10000000000)==36)]
    
    ## add proper identifiers, population size, and outcome data (where available) into the matched dataset ##
    matchdat$ID<-c(eedat$dohregion,ny2010$dohregion[match(controlid_ny,ny2010$geoid)],controlid)
    matchdat$Pop<-c(eedat$Pop,ny2010$Pop[match(controlid_ny,ny2010$geoid)],controldat$Pop[match(controlid,controldat$ID)])
    matchdat$Kidney<-c(eedat$observed_Kidney,ny2010$observed_Kidney[match(controlid_ny,ny2010$geoid)],rep(NA,length(controlid)))
    matchdat$Bladder<-c(eedat$observed_Bladder,ny2010$observed_Bladder[match(controlid_ny,ny2010$geoid)],rep(NA,length(controlid)))
    
    ##########################################
    ## cancer predictions for matched areas ##
    ##########################################
    for (ctype in c("Kidney","Bladder")){
      
      ## load the posterior samples from the prediction model ##
      load(paste0(tolower(ctype),'_output.RData'))
      ps_beta<-fit1[[1]]
      
      ## load the SEER county cancer totals for 2005-2009 ##
      ctype_ct<-read.csv(paste0(ctype,'_County_Ct.csv'),stringsAsFactors = F)
      
      ## add the NY county counts ##
      nyctype1<-tapply(ny2010[[paste0('observed_',ctype)]],ny2010$county,sum)
      nyctype2<-cbind(as.numeric(names(nyctype1)),nyctype1)
      colnames(nyctype2)<-names(ctype_ct)
      rownames(nyctype2)<-NULL
      ctype_ct<-rbind(ctype_ct,nyctype2)
      
      ## extract confounder and pop data from all counties used in the matched data ##
      counties<-floor(as.numeric(controlid)/10000000)
      XTseer$county<-floor(XTseer$ID/10000000)
      all_countydat<-XTseer[which(XTseer$county %in% counties),]
      X_countydat<-all_countydat[,c(2,5:9,12:16)]
      
      ## do prediction of CBG cancer counts for each county in this dataset ##
      Ypred<-NULL
      for (i in unique(all_countydat$county)){
        
        ## county i predictors for the prediction model ##
        temp<-as.matrix(X_countydat[which(all_countydat$county==i),])
        ## ids of the matched units in county i ##
        temp2<-controlid[which(counties==i)]
        ## ids of all units in county i ##
        temp3<-all_countydat$ID[which(all_countydat$county==i)]
        ## locations of matched units in the county i data ##
        temp4<-which(temp3 %in% temp2)
        
        ## matrix to save ppd samples for the matched units in county i ##
        cbg_pred<-matrix(NA,nrow=length(temp2),ncol=nrow(ps_beta))
        rownames(cbg_pred)<-temp3[temp4]
        
        ## begin ppd sampling ##
        for (j in 1:nrow(ps_beta)){
          ## get estimated probabilities for each CBG ##
          pi<-exp(temp%*%matrix(ps_beta[j,]))/sum(exp(temp%*%matrix(ps_beta[j,])))
          ## draw posterior predictive sample from multinomial distribution ##
          foo<-rmultinom(n=1,size=ctype_ct[[paste0(ctype,'_Ct')]][which(ctype_ct$County==i)],prob=pi)
          ## save the predicted values for the matched units ##
          cbg_pred[,j]<-foo[temp4]
        }
        
        ## concatenate ppd samples from all matched units ##
        Ypred<-rbind(Ypred,cbg_pred)
      }
      
      ## add the NY outcome data to the predictions ##
      temp<-matrix(rep(matchdat[[ctype]][1:(nrow(eedat)+length(controlid_ny))],nrow(ps_beta)),
                   nrow=nrow(eedat)+length(controlid_ny),ncol=ncol(Ypred),byrow=F)
      rownames(temp)<-matchdat$ID[1:(nrow(eedat)+length(controlid_ny))]
      
      ## make sure ordering is the same in outcome and predictor data ##
      Ypred<-Ypred[order(as.numeric(rownames(Ypred))),]
      Ypred<-rbind(temp,Ypred)
      ## write the prediction info to an external dataset ##
      predsum<-data.frame(apply(Ypred,1,mean),apply(Ypred,1,quantile,.025),apply(Ypred,1,quantile,.975))
      write.csv(predsum, file=paste0(ctype,'_predsum_',M,'_',metho,'_sa2.csv'),row.names=F)
      temp<-matchdat[(nrow(eedat)+length(controlid_ny)+1):nrow(matchdat),]
      temp<-temp[order(as.numeric(temp$ID)),]
      matchdat[(nrow(eedat)+length(controlid_ny)+1):nrow(matchdat),]<-temp
      
      #####################################################
      ## fit the poisson regression to estimate the cSIR ##
      #####################################################
      
      X_est<-as.matrix(matchdat[,1:ncol(matchin)])
      
      C<-75*vcov(glm(round(apply(Ypred,1,mean))~X_est,family='poisson'))
      
      accept<-1
      tuneparam<-.02
      
      while(accept<.2 | accept>.4){
        if (accept<.2) tuneparam<-tuneparam-.001
        if (accept>.4) tuneparam<-tuneparam+.001
        fit_ctype<-est_poisson(X=X_est,Y_all=Ypred,offs=matchdat$Pop,tau=tuneparam,pburn=.75,C=C)
        accept<-fit_ctype[['acceptance']]
      }
      
      ## output traceplots ##
      pdf(paste0(ctype,'_traceplots_',M,'_',metho,'_sa2.pdf'),width=12)
      par(mfrow=c(2,6))
      for (i in 1:ncol(fit_ctype[['beta_postsamp']])){
        traceplot(mcmc(fit_ctype[['beta_postsamp']][,i]))
      }
      dev.off()
      
      ## output results ##
      if (M==3 & metho=="mahalanobis"){
        sink(paste0(ctype,'_results_sa2.txt'))
        cat('\n')
        cat("=============================\n")
        cat('M=',M,', Method=',metho,'\n')
        cat("=============================\n")
        
        cat('Acceptance=',fit_ctype[['acceptance']],'\n')
        cat('Model output:\n')
        print(data.frame(fit_ctype[['ebeta_postsum']],c(1,colnames(X_est))))
        cat('\n')
        sink()
      } else{
        sink(paste0(ctype,'_results_sa2.txt'),append=T)
        cat('\n')
        cat("=============================\n")
        cat('M=',M,', Method=',metho,'\n')
        cat("=============================\n")
        
        cat('Acceptance=',fit_ctype[['acceptance']],'\n')
        cat('Model output:\n')
        print(data.frame(fit_ctype[['ebeta_postsum']],c(1,colnames(X_est))))
        cat('\n')
        sink()
      }
    }
  }
}
