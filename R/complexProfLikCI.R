#Used by complexOadaProfLik
complexSplitParamsOadaLikelihood<-function(parVect, which, value,complexNBDAdata,
                                           transmissionFunction=function(par,connectionToInformed) sum(par[1]*connectionToInformed),sumRateAcrossNetworks=T){
  newParVect<-rep(NA,length(parVect)+1)
  newParVect[-which]<-parVect
  newParVect[which]<-value
  complexOadaLikelihood(parVect=newParVect, complexNBDAdata=complexNBDAdata,
                        transmissionFunction=transmissionFunction,sumRateAcrossNetworks=sumRateAcrossNetworks)
}

#Used by complexOadaProfLik
complexSplitParamsOadaGradient<-function(parVect, which, value,complexNBDAdata,
                                           transmissionFunction=function(par,connectionToInformed) sum(par[1]*connectionToInformed),sumRateAcrossNetworks=T){
  grad(complexSplitParamsOadaLikelihood,parVect,which=which,value=value,complexNBDAdata=complexNBDAdata,transmissionFunction=transmissionFunction,sumRateAcrossNetworks=sumRateAcrossNetworks)
}

#Finds the profile log likelihood for a specified parameter which, set to value. i.e. it optimises the other parameters to give the maximum likelihood
complexOadaProfLik<-function(value, which,complexNBDAdata,
                             transmissionFunction=function(par,connectionToInformed) sum(par[1]*connectionToInformed),sumRateAcrossNetworks=T,startValue=NULL, lower=NULL,noParTransFunct,
                             iterations=150){

    nbdadata<- complexNBDAdata

    #If a multiple diffusion is specified as a character vector, convert to a list (for backwards compatibility)
     if(is.character(nbdadata)){
       newNbdaData<-list()
       for(i in 1:length(nbdadata)){
         newNbdaData<-c(newNbdaData,list(eval(as.name(nbdadata[i]))))
       }
       nbdadata<-newNbdaData
     }
     #If there are multiple diffusions "borrow" the first diffusion to extract necessary parameters
     if(is.list(nbdadata)){
       nbdadataTemp<-nbdadata[[1]]
     }else{nbdadataTemp<-nbdadata}


     #calculate the number of each type of parameter
     noSParam <- noParTransFunct #s parameters
     noILVasoc<- dim(nbdadataTemp@asocILVdata)[2] #ILV effects on asocial learning
     noILVint<- dim(nbdadataTemp@intILVdata)[2] #ILV effects on interaction (social learning)
     noILVmulti<- dim(nbdadataTemp@multiILVdata)[2] #ILV multiplicative model effects

     if(nbdadataTemp@asoc_ilv[1]=="ILVabsent"){noILVasoc<-0} #Ignore dummy ILV
     if(nbdadataTemp@int_ilv[1]=="ILVabsent"){noILVint<-0} #Ignore dummy ILV
     if(nbdadataTemp@multi_ilv[1]=="ILVabsent"){noILVmulti<-0} #Ignore dummy ILV

     #Record asocialVar names
     if(nbdadataTemp@asoc_ilv[1]=="ILVabsent"){asocialVarNames<-NULL}else{asocialVarNames<-nbdadataTemp@asoc_ilv};
     if(nbdadataTemp@int_ilv[1]=="ILVabsent"){intVarNames<-NULL}else{intVarNames<-nbdadataTemp@int_ilv};
     if(nbdadataTemp@multi_ilv[1]=="ILVabsent"){multiVarNames<-NULL}else{multiVarNames<-nbdadataTemp@multi_ilv};

     #Set staring values if not specified by the user
     if(is.null(startValue)) {
       NewStartValue<-rep(0,noSParam+noILVasoc+noILVint+noILVmulti);
       if(is.null(lower)){
         newLower<-rep(-Inf,length(NewStartValue));
         newLower[1:noSParam]<-0;
       }
       upper<-rep(Inf,length(NewStartValue));
     }else{
       if(is.null(lower)){
         newLower<-rep(-Inf,length(startValue));
         newLower[1:noSParam]<-0;
      }
      upper<-rep(Inf,length(startValue));
     }

     #Remove "which" parameter from the startValue, upper and lower vectors
     if(is.null(startValue)) startValue<-NewStartValue[-which]
     if(is.null(lower)) lower<-newLower[-which]
     upper<-upper[-which]

     if(length(startValue)==0) {
       return(complexOadaLikelihood(parVect=value, complexNBDAdata=complexNBDAdata,
                                                             transmissionFunction=transmissionFunction,sumRateAcrossNetworks=sumRateAcrossNetworks,
                                                             noParTransFunct=1))
     }

     #Optimise other parameters while setting which to value
     try(fit1<-nlminb(start=startValue, objective= complexSplitParamsOadaLikelihood, gradient=complexSplitParamsOadaGradient,
                  transmissionFunction = transmissionFunction,sumRateAcrossNetworks=sumRateAcrossNetworks, which= which, value=value,
                  lower=lower, upper=upper,complexNBDAdata=nbdadata,control=list(iter.max=iterations)))
     #Return the prof log lik and the indicator as to whether the there was convergence
     return(c(fit1$objective,fit1$convergence))
}

#Plots the profile log likelihood for a specified parameter, which, over a specified range, from a specified complexOadaFit model
#The dashed line gives the cutoff point for the specified confidence interval (conf)
complexPlotProfLik<-function(which,model,range,resolution=20,inflation=1,conf=0.95,startValue=NULL,lower=NULL,iterations=150){

  if(model@type=="asocial"){
    print("Function does not yet work for asocial models")
    return(NULL)
  }

  complexNBDAdata<-model@nbdadata
  noParTransFunct<-model@noParTransFunct

  xVals<-seq(range[1],range[2],length=resolution)
  profLik<-matrix(NA,nrow=length(xVals),ncol=2)

  cutoff<-model@loglik+inflation*qchisq(conf,1)/2

  if(is.null(startValue)) startValue<-model@outputPar[-which]

  for(i in 1:length(xVals)){
    profLik[i,]<-complexOadaProfLik(which=which, value=xVals[i],complexNBDAdata=complexNBDAdata,startValue=startValue,
                                    lower=lower,noParTransFunct=noParTransFunct,transmissionFunction=model@transmissionFunction,iterations=iterations)
    plot(xVals,profLik[,1],type="l",xlim=range,ylim=c(model@loglik-(max(na.omit(profLik))-model@loglik)*0.03,max(na.omit(profLik))), col=profLik[,2]+1,
         xlab=model@varNames[which],ylab="Profile log-likelihood")
    abline(h= cutoff, lty=2)
  }

  return(data.frame(xVals,profLik=profLik[,1],convergence=profLik[,2]))

}

complexDistanceFromCutoff<-function(value,which,model,inflation=1,conf=0.95,startValue=NULL,lowerIn=NULL,iterations=150){

  if(model@type=="asocial"){
    print("Function does not yet work for asocial models")
    return(NULL)
  }

  complexNBDAdata<-model@nbdadata
  noParTransFunct<-model@noParTransFunct

  cutoff<-model@loglik+inflation*qchisq(conf,1)/2

  if(is.null(startValue)) startValue<-model@outputPar[-which]

  profLik<-complexOadaProfLik(which=which, value=value,complexNBDAdata=complexNBDAdata,startValue=startValue,
                                    lower=lowerIn,noParTransFunct=noParTransFunct,transmissionFunction=model@transmissionFunction,iterations=iterations)[1]
  return(abs(profLik-cutoff))

}

complexProfLikCI<-function(which,model,interval,inflation=1,conf=0.95,startValue=NULL,lowerIn=NULL,iterations=150){

  if(model@type=="asocial"){
    print("Function does not yet work for asocial models")
    return(NULL)
  }

  complexNBDAdata<-model@nbdadata
  noParTransFunct<-model@noParTransFunct

  cutoff<-model@loglik+inflation*qchisq(conf,1)/2

  if(is.null(startValue)) startValue<-model@outputPar[-which]

  model<-optimise(f=complexDistanceFromCutoff,interval=interval,which=which,model=model,startValue=startValue,
           lowerIn=lowerIn,iterations=iterations,inflation=inflation)

  c(model$minimum, model$objective)

}
