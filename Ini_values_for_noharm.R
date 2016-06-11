################################################################
#
#	INIVALUES FOR NOHARM BASED ON UNIDIMENSIONAL
#				MODELS
#	
################################################################
#	N<-dim(LSAT)[1]

#	dat	=	Data
#	nc	=	Number of Clusters
#	d	= 	Dimension
#
# 	For now we suppose that d = nc



calculate_init_vals_unidim = function(data, size.cluster){
      
    
    ################################################################
    #			 Separate Clusters
    ################################################################
    
    
    #################################3
    # SEPARATE CLUSTER
    
    CLUSTER<- list()
    
    start = 1
    for(i in 1:length(size.cluster))
    {
      CLUSTER[[i]]<- data[,start:(start + (size.cluster[i] - 1 ))]
      start = start + size.cluster[i] 
    }
    
    
    #################################
    # REORGANIZE CLUSTERS
    
    for(i in 1:length(size.cluster))
    {
      ACP<-PCA(CLUSTER[[i]], graph=F)
      axe.1<- which(ACP$var$coord[,1]==max(ACP$var$coord[,1]))
      CLUSTER[[i]]<- cbind(CLUSTER[[i]][,axe.1], CLUSTER[[i]][,-axe.1])
    }
    
    
    
    ################################################################
    #		Recover Unidimensional IRT models to obtain
    #			Initial values for NOHARM
    ################################################################
    
    
    fit_ltm<- list()
    traits<- list()
    for(i in 1:length(size.cluster)) # == for(i in 1:nc)
    {
      fit_ltm[[i]]<- tpm(CLUSTER[[i]])
      traits[[i]]<- factor.scores(fit_ltm[[i]], resp.patterns = CLUSTER[[i]])
    }
    
    #################################################################
    #	Transform "b" (difficultie parameter) into "d"
    #################################################################
    AUX1<-lapply(fit_ltm,function(X) -coef(X)[,2]*coef(X)[,3])
    COEF<- lapply(fit_ltm, coef)
    
    for(i in 1:length(size.cluster))
    {
      COEF[[i]][,2] <- AUX1[[i]]
    }
    
    
    #################################################################
    #	Found Correlation Matrix and Transform A into 
    # 			A^*
    #################################################################
    
    
    THETA<- matrix(NA, ncol = length(size.cluster), nrow = nrow(data))
    A<- list()
    
    for(i in 1:length(size.cluster))
    {
      THETA[,i]<- traits[[i]]$score.dat$z1
      A[[i]]<- COEF[[i]][,3]
    }
    AUX2<- matrix(NA, ncol=length(size.cluster), nrow= sum(size.cluster))
    for(i in 1:length(size.cluster))
    {
      col=i
      if(i==1)
      {
        AUX2[,col] <- c(as.numeric(c(A[[i]])), rep(NA, length(AUX2[,1]) -  length(c(A[[i]]))))
      }
      if(i==length(size.cluster))
      {
        AUX2[,col] <- c(rep(NA, length(AUX2[,1]) -  length(c(A[[i]]))), as.numeric(c(A[[i]])))
      }
      if(i != 1 &&  i!= length(size.cluster))
      {
        AUX2[,col]<- c(rep(NA, sum(size.cluster[1:(col-1)])),
                       as.numeric(c(A[[i]])),
                       rep(NA, sum(size.cluster[(col+1):length(size.cluster)])))
      }
    }
    A_matrix<- ifelse(is.na(AUX2), 0,AUX2)
    SIGMA<- cov(THETA)
    CORR<- cor(THETA)
    # CORR<- matrix(c(1,.5,.5,1), ncol=2)
    # SIGMA<- matrix(c(.9,.8,.8,.9), ncol=2)
    
    CholCorr<- chol(CORR) 
    CholSigm<- chol(SIGMA)
    
    #t(chol(CORR))%*%chol(CORR)
    #t(chol(SIGMA))%*%chol(SIGMA)
    
    A_asterisco<- A_matrix%*%CholSigm%*%solve(CholCorr)
    A_asterisco[is.na(AUX2)] = 0
    
    #	A%*%t(THETA)
    #	A_asterisco%*%CholCorr%*%solve(CholSigm)%*%t(THETA)
    list(coefs = cbind(A_asterisco, rep(0, nrow(A_asterisco)),rep(0, nrow(A_asterisco))), corr = CORR)
}

