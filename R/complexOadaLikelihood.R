
#Only works with static networks at the moment

complexOadaLikelihood <- function(parVect, complexNBDAdata,transmissionFunction=function(par,connectionToInformed) sum(par[1]*connectionToInformed),sumRateAcrossNetworks=T){

  complexNBDAdata->nbdadata

  #If there are multiple diffusions "borrow" the first diffusion to extract necessary parameters
  if(is.list(nbdadata)){
    totalLikelihood<-0
    for(i in 1:length(nbdadata)){
      subdata<-nbdadata[[i]];
      totalLikelihood<-totalLikelihood+complexOadaLikelihood(parVect,subdata,transmissionFunction=transmissionFunction,sumRateAcrossNetworks=sumRateAcrossNetworks)
    }
    return(totalLikelihood)
  }else{

           noILVasoc<- dim(nbdadata@asocILVdata)[2] #ILV effects on asocial learning
           noILVint<- dim(nbdadata@intILVdata)[2] #ILV effects on interaction (social learning)
           noILVmulti<- dim(nbdadata@multiILVdata)[2] #ILV multiplicative model effects

           if(nbdadata@asoc_ilv[1]=="ILVabsent"){noILVasoc<-0} #Ignore dummy ILV
           if(nbdadata@int_ilv[1]=="ILVabsent"){noILVint<-0} #Ignore dummy ILV
           if(nbdadata@multi_ilv[1]=="ILVabsent"){noILVmulti<-0} #Ignore dummy ILV

           #Number of parameters provided for transmission function
           noTransFuncPar<-length(parVect) -noILVasoc-noILVint-noILVmulti
           #Split the parameters between the transmission function and normal NBDA
           transFuncPar<-parVect[1:noTransFuncPar]
           otherPar<-parVect[-(1:noTransFuncPar)]

           #If the transmission function has an argument named connectionToInformed
           #Then this vector is provided to generate the rates
           #Otherwise, it is assumed that the arguments connectionToAll and statusOfOthers are present
           #And these are provided to generate the rate
           #The latter is more flexible but the former is quicker where it can be used

           if(sum(formalArgs(transmissionFunction)=="connectionToInformed")>0){

              #Calculate network connections x status- i.e. take out connections to uninformed
              connectionToInformed<-nbdadata@stMetric
              for(i in 1:dim(nbdadata@stMetric)[3]){
                connectionToInformed[,,i]<-nbdadata@stMetric[,,i]*nbdadata@statusOfOthers
              }

              if(sumRateAcrossNetworks){
                #The default is to assume that the transmission function is applied to each network and then added
                unscaled.st<-apply(connectionToInformed,c(1,3),transmissionFunction,par=transFuncPar)
                if(dim(nbdadata@assMatrix)[3]>1) unscaled.st<-apply(unscaled.st,1,sum)
              }else{

                #But if the user specifies sumRateAcrossNetworks=F, they can provide a function that
                #a N x n matrix as its network argument (n individuals, N networks) and returns a rate
                unscaled.st<-apply(connectionToInformed,1,transmissionFunction,par=transFuncPar)

              }
           }else{

            #Extract network connections and status
            connectionToAll<-nbdadata@stMetric
            statusOfOthers<-nbdadata@statusOfOthers
            connectionsAndStatus<-array(NA,dim=c(dim(connectionToAll),2))
            connectionsAndStatus[,,,1]<-connectionToAll
            for (j in 1:dim(connectionToAll)[3]){
                 connectionsAndStatus[,,j,2]<-statusOfOthers
            }

            if(sumRateAcrossNetworks){
                 #The default is to assume that the transmission function is applied to each network and then added
                 unscaled.st<-apply(connectionsAndStatus,c(1,3),transmissionFunction,par=transFuncPar)
                 if(dim(connectionToAll)[3]>1) unscaled.st<-apply(unscaled.st,1,sum)
            }else{
                 #But if the user specifies sumRateAcrossNetworks=F, they can provide a function that
                 #a N x n matrix as its network argument (n individuals, N networks) and returns a rate
                 #Cycle through individuals x event combinations and calculate the rate each time
                 unscaled.st<-apply(connectionsAndStatus,1,transmissionFunction,par=transFuncPar)
               }
             }


           #Transform data into an nbdaData object with stMetric as unscaled.st above (this will be constrained to have an s parameter of 1)
           transformedData<-nbdadata
           transformedData@stMetric=matrix(unscaled.st,ncol=1)
           #Strictly speaking this is the wrong class of object for this slot in a complexNBDAdata object,
           #but we should get away with it since we are going to use it as an nbdaData object, and it will only exist inside the function

           return(oadaLikelihood(par=c(1,otherPar),transformedData))
           }
  }

complexOadaGradient<- function(parVect, complexNBDAdata,transmissionFunction=function(par,connectionToInformed) sum(par[1]*connectionToInformed),sumRateAcrossNetworks=T){

  grad(complexOadaLikelihood,parVect,complexNBDAdata=complexNBDAdata,transmissionFunction=transmissionFunction,sumRateAcrossNetworks=sumRateAcrossNetworks)

}
