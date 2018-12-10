cdc_sir<-function(Ycom,POPcom,Yseer,POPseer){
  ## observed number of cases in the community of interest ##
  Ycom_tot<-sum(Ycom)
  ## rate of cancer in seer ##
  lambda<-sum(Yseer)/sum(POPseer)
  ## multiply lambda by the population in the community of interest to get expected value ##
  E<-lambda*sum(POPcom)
  SIR<-Ycom_tot/E
  SIRll<-qchisq(.025,df=2*Ycom_tot)/(2*E)
  SIRul<-qchisq(.975,df=2*(Ycom_tot+1))/(2*E)
  return(c(SIR,SIRll,SIRul))
}