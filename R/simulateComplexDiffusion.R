

simulateComplexDiffusion_OADA<-function(par,network,transmissionFunction=function(par,connectionToInformed) sum(par[1]*connectionToInformed),sumRateAcrossNetworks=T,
                                        demons=rep(0,dim(network)[1]), ilv="ILVabsent",asoc_coef=NULL, int_coef=NULL,
                                        label="simulatedDiffusion"){

  #Convert network to 3d array if needed
  if(length(dim(network))==2){
    network<-array(network,dim=c(dim(network),1))
  }

  #If provided, calculate the LP for the effect of ILVs
  if(ilv[1]!="ILVabsent"){
  # If there is just one asocial variable matrix for all events and times, then you will have a column matrix for each ILV, the length of the number of individuals
  ilv.dim <- dim(eval(as.name(ilv[1])))[1] # specify the dimensions of assoc.array
  # create asoc.array to hold the asocial variables, needs to be treatment  "constant", create a one-matrix array
  ilv.array <- array(dim=c(ilv.dim, 1, length(ilv)))
  dimnames(ilv.array) <- list(NULL, NULL, ilv)

  for(a in 1:length(ilv)){ # Loop through asocial variables - a loop
    ilv.array[,,a] <- eval(as.name(ilv[a])) # evaluate each one in turn
  }
  }
  asocialLP<-socialLP<-rep(0, dim(network)[1])
  if(ilv[1]!="ILVabsent"){
    asocialLP<-as.matrix(ilv.array[,1,])%*%asoc_coef
    socialLP<-as.matrix(ilv.array[,1,])%*%int_coef
  }

  #Get starting status
  status<-demons

  orderAcq<-totalRatesRecord<-rep(NA,dim(network)[1]-sum(demons))

  i<-0
  while(sum(status==0)>0){
    i<-i+1


    #If the transmission function has an argument named connectionToInformed
    #Then this vector is provided to generate the rates
    #Otherwise, it is assumed that the arguments connectionToAll and statusOfOthers are present
    #And these are provided to generate the rate
    #The latter is more flexible but the former is quicker where it can be used
    if(sum(formalArgs(transmissionFunction)=="connectionToInformed")>0){

      #Calculate network connections x status- i.e. take out connections to uninformed
      connectionToInformed<-network
      for(j in 1:dim(network)[3]){
        for(k in 1:dim(network)[1]){
          connectionToInformed[k,,j]<-network[k,,j]*status
        }
      }

      if(sumRateAcrossNetworks){
        #The default is to assume that the transmission function is applied to each network and then added
        unscaled.st<-apply(connectionToInformed,c(1,3),transmissionFunction,par=par)
        if(dim(network)[3]>1) unscaled.st<-apply(unscaled.st,1,sum)
      }else{

        #But if the user specifies sumRateAcrossNetworks=F, they can provide a function that
        #a N x n matrix as its network argument (n individuals, N networks) and returns a rate
        unscaled.st<-apply(connectionToInformed,1,transmissionFunction,par=par)

      }
    }else{

      #Extract network connections and status
      connectionToAll<-network
      statusOfOthers<-status

      connectionsAndStatus<-array(NA,dim=c(dim(connectionToAll),2))
      connectionsAndStatus[,,,1]<-connectionToAll
      for (j in 1:dim(connectionToAll)[3]){
        connectionsAndStatus[,,j,2]<-matrix(rep(statusOfOthers,dim(connectionToAll)[1]),nrow=dim(connectionToAll)[1],byrow = T)
      }

      if(sumRateAcrossNetworks){
        #The default is to assume that the transmission function is applied to each network and then added
        unscaled.st<-apply(connectionsAndStatus,c(1,3),transmissionFunction,par=par)
        if(dim(network)[3]>1) unscaled.st<-apply(unscaled.st,1,sum)
      }else{
        #But if the user specifies sumRateAcrossNetworks=F, they can provide a function that
        #a N x n matrix as its network argument (n individuals, N networks) and returns a rate
        #Cycle through individuals x event combinations and calculate the rate each time
        unscaled.st<-apply(connectionsAndStatus,1,transmissionFunction,par=par)
      }
    }

    #Calculate total (relative) rate
    totalRate <- (exp(asocialLP) + exp(socialLP) * unscaled.st)*(1-status)
    #Calculate the probability each will be the next to learn
    probs<- totalRate/sum(totalRate)

    #Simulate the next individual to learn
    nextToLearn<- which(rmultinom(n=1,size=1,prob=probs)==1)
    #Record in orderAcq vector
    orderAcq[i]<-nextToLearn
    #Update status
    status[nextToLearn]<-1
    #Record the total relative rate, since this can used later to generate times of acquisition
    totalRatesRecord[i]<-sum(totalRate)
  }

  #Convert network to 4d array if needed
    network<-array(network,dim=c(dim(network),1))

    simData<-complexNBDAdata(label=label, assMatrix=network, asoc_ilv = ilv,
                    int_ilv = ilv, orderAcq=orderAcq,
                    demons = demons,asocialTreatment = "constant")

  return(simData)

}

