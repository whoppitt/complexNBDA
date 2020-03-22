#Define class of object for the fitted additive model
setClass("complexOadaFit",representation(nbdaMultiDiff="character",nbdadata="list",transmissionFunction="function",noParTransFunct="numeric",sumRateAcrossNetworks="logical",
                                         optimisation="list",optim="list",loglik="numeric",aic="numeric",aicc="numeric",varNames="character",
                                         hessian="matrix",outputPar="numeric",se="numeric",type="character",SLdom="logical"));


#Method for initializing addFit object- including model fitting
setMethod("initialize",
          signature(.Object = "complexOadaFit"),
          function (.Object, nbdadata,type,startValue,lower,method,interval,gradient,iterations,standardErrors,transmissionFunction,noParTransFunct,sumRateAcrossNetworks,SLdom,...)
          {

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
              if(is.null(startValue)) startValue<-rep(0,noSParam+noILVasoc+noILVint+noILVmulti);

              #Set vector of upper and lower values for each parameter. Unbounded for asocial learning variables
              if(is.null(lower)){
              lower<-rep(-Inf,length(startValue));
              lower[1:noSParam]<-0;
              }
              upper<-rep(Inf,length(startValue));
              if(is.null(interval)) interval<-c(0,999);

              #Optimise
              #All being done with nlminb at the moment
              fit1<-NULL
              if(gradient){
                try(fit1<-nlminb(start=startValue, objective= complexOadaLikelihood,gradient= complexOadaGradient, transmissionFunction = transmissionFunction,sumRateAcrossNetworks=sumRateAcrossNetworks,lower=lower, upper=upper,complexNBDAdata=nbdadata,control=list(iter.max=iterations)));
              }else{
                try(fit1<-nlminb(start=startValue, objective= complexOadaLikelihood,transmissionFunction = transmissionFunction,sumRateAcrossNetworks=sumRateAcrossNetworks,lower=lower, upper=upper,complexNBDAdata=nbdadata,control=list(iter.max=iterations)));
              }
              if(is.null(fit1)){
                print("Error in likeihood optimization");
                return(NULL)
              }

              #SEs set to numeric since a custom Hessian would be required
              standardErrors<-"Numeric"

              if(method=="both"){
                if (is.null(fit1)){
                  try(fit2<-optim(par=fit1$par,fn=complexOadaLikelihood,transmissionFunction = transmissionFunction,sumRateAcrossNetworks=sumRateAcrossNetworks,method="L-BFGS-B",gr=complexOadaGradient,hessian=T,lower=lower, upper=upper,complexNBDAdata=nbdadata,control=list(maxit=iterations)))
                }else{
                  try(fit2<-optim(par=startValue,fn=complexOadaLikelihood,transmissionFunction = transmissionFunction,sumRateAcrossNetworks=sumRateAcrossNetworks,method="L-BFGS-B",gr=complexOadaGradient,hessian=T,lower=lower, upper=upper,complexNBDAdata=nbdadata,control=list(maxit=iterations)))
                }
              }else{fit2<-as.list(NA)}


              #Record MLEs
              outputPar<-fit1$par;

              #Perform LRT for social transmission
              loglik<-fit1$objective;

              #Get sample size across all diffusions
              sampleSize<-sampSizeExtract_complex(list(nbdadata));


              #Calculate aic and aicc
              aic<-2*length(fit1$par)+2*loglik;
              aicc<-2*(length(fit1$par))*(sampleSize/(sampleSize-(length(fit1$par))-1))+2*loglik;

              #To prevent a low AICc when there are more parameters than data!
              if(is.nan(aic)|is.nan(aicc)){}else{
                if(aicc<aic) aicc<-Inf;
              }


              #		}

              #Extract names of variables
              parCounter<-0
              varNames<-paste(parCounter+(1:noSParam),"Transmission parameter",1:noSParam)
              parCounter<-parCounter+noSParam
              if(!is.null(asocialVarNames)){
                varNames<-c(varNames,paste(parCounter+(1:length(asocialVarNames)),"Asocial:",asocialVarNames))
                parCounter<-parCounter+length(asocialVarNames)
              }
              if(!is.null(intVarNames)){
                varNames<-c(varNames,paste(parCounter+(1:length(intVarNames)),"Social:",intVarNames))
                parCounter<-parCounter+length(intVarNames)
              }
              if(!is.null(multiVarNames)){
                varNames<-c(varNames,paste(parCounter+(1:length(multiVarNames)),"Social= asocial:",multiVarNames))
              }


              #Get hessian matrix and use it to get standard errors
              hessianMat<- hessian(func=complexOadaLikelihood,x=fit1$par,complexNBDAdata=nbdadata,transmissionFunction=transmissionFunction,sumRateAcrossNetworks=sumRateAcrossNetworks)

              if(is.null(hessianMat)){
                se<-rep(NaN,length(outputPar))
                hessianMat<-matrix(NA)
              }else{
                seTemp<-NULL
                if(det(hessianMat)==0|is.na(det(hessianMat))){
                  #If the hessian matrix is not positive definite, SEs cannot be calculated
                  se<-rep(NaN,length(outputPar));
                }else{
                  #Initialise varTemp to NaN so it is found if the hessianMat cannot be inverted
                  varTemp<-diag(hessianMat)
                  varTemp[]<-NaN
                  try(varTemp<-diag(solve(hessianMat)),silent=T);
                  varTemp[varTemp<0]<-NaN
                  seTemp<-sqrt(varTemp)
                  #Record SEs giving a zero for constrained values
                  if(is.null(seTemp)){
                    se<-rep(NaN,length(outputPar));
                  }else{
                    se<-seTemp;
                  }
                }
                if(!is.matrix(hessianMat))hessianMat<-matrix(hessianMat)
              }

              if(is.list(nbdadata)){
                callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = nbdadata, transmissionFunction=transmissionFunction,noParTransFunct=noParTransFunct,sumRateAcrossNetworks=sumRateAcrossNetworks,
                               optimisation=fit1,optim=fit2,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,SLdom=SLdom,...)
              }else{
                callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = list(nbdadata), transmissionFunction=transmissionFunction,noParTransFunct=noParTransFunct,sumRateAcrossNetworks=sumRateAcrossNetworks,
                               optimisation=fit1,optim=fit2,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,SLdom=SLdom,...)

              }

          }
            )

#Function for implementing the initialization
complexOadaFit<-function(complexNbdaData,type="social",startValue=NULL, lower=NULL,interval=c(0,999), method="nlminb", gradient=T,iterations=150,
                  transmissionFunction=function(par,connectionToInformed) sum(par[1]*connectionToInformed), noParTransFunct=1,sumRateAcrossNetworks=T){

  nbdadata<-complexNbdaData
    #If a multiple diffusion is specified as a character vector, convert to a list (for backwards compatibility)
    if(is.character(nbdadata)){
      newNbdaData<-list()
      for(i in 1:length(nbdadata)){
        newNbdaData<-c(newNbdaData,list(eval(as.name(nbdadata[i]))))
      }
      nbdadata<-newNbdaData
    }

      return(new("complexOadaFit",nbdadata= complexNbdaData,transmissionFunction=transmissionFunction,noParTransFunct=noParTransFunct,type= "social", startValue= startValue,
                 lower=lower,interval= interval,method= method,gradient=gradient,iterations=iterations, standardErrors="Numeric",SLdom=FALSE,sumRateAcrossNetworks=sumRateAcrossNetworks))
}


sampSizeExtract_complex<-function (complexnbdadata)
{
  complexnbdadata->nbdadata
  if (class(nbdadata) == "complexNBDAdata") {
    return(sum(nbdadata@status))
  }
  if (is.list(nbdadata)) {
    totalSS <- 0
    for (i in 1:length(nbdadata)) {
      totalSS <- totalSS + sampSizeExtract_complex(nbdadata[[i]])
    }
    return(totalSS)
  }
}
