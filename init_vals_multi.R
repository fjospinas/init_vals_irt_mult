#Calculate init unidimensional vals
calculate_init_vals_unidim = function(data, size.cluster, verbose = FALSE){
  
  if(verbose) print("Start process to calculate multidimensional init vals.")
  
  #separate cluster
  cluster<- list()
  start = 1
  for(i in 1:length(size.cluster)){
    cluster[[i]]<- data[,start:(start + (size.cluster[i] - 1 ))]
    start = start + size.cluster[i] 
  }

  #REORGANIZE clusterS
  for(i in 1:length(size.cluster)){
    acp<-PCA(cluster[[i]], graph=F)
    axe.1<- which(acp$var$coord[,1]==max(acp$var$coord[,1]))
    cluster[[i]]<- cbind(cluster[[i]][,axe.1], cluster[[i]][,-axe.1])
  }
  
  #Recover Unidimensional IRT models to obtain
  #Initial values for NOHARM
  fit_ltm<- list()
  traits<- list()
  for(i in 1:length(size.cluster)){ # == for(i in 1:nc){
    fit_ltm[[i]]<- tpm(cluster[[i]])
    traits[[i]]<- factor.scores(fit_ltm[[i]], resp.patterns = cluster[[i]])
  }
  
  #	Transform "b" (difficultie parameter) into "d"
  aux1<-lapply(fit_ltm,function(x) -coef(x)[,2]*coef(x)[,3])
  coef<- lapply(fit_ltm, coef)
  for(i in 1:length(size.cluster)){
    coef[[i]][,2] <- aux1[[i]]
  }
  
  #	Found Correlation Matrix and Transform A into A^*
  #and build d and c
  theta<- matrix(NA, ncol = length(size.cluster), nrow = nrow(data))
  coef_d = coef_c = c()
  A<- list()
  for(i in 1:length(size.cluster)){
    theta[,i]<- traits[[i]]$score.dat$z1
    A[[i]]<- coef[[i]][,3]
    coef_d = c(coef_d, coef[[i]][,2])
    coef_c = c(coef_c, coef[[i]][,1])
  }
  
  aux2<- matrix(NA, ncol=length(size.cluster), nrow= sum(size.cluster))
  for(i in 1:length(size.cluster)){
    col=i
    if(i==1){
      aux2[,col] <- c(as.numeric(c(A[[i]])), rep(NA, length(aux2[,1]) -  length(c(A[[i]]))))
    }else if(i == length(size.cluster)){
      aux2[,col] <- c(rep(NA, length(aux2[,1]) -  length(c(A[[i]]))), as.numeric(c(A[[i]])))
    }else if(i != 1 &&  i!= length(size.cluster)){
      aux2[,col]<- c(rep(NA, sum(size.cluster[1:(col-1)])),
                     as.numeric(c(A[[i]])),
                     rep(NA, sum(size.cluster[(col+1):length(size.cluster)])))
      }
    }
  
  A_matrix<- ifelse(is.na(aux2), 0,aux2)
  sigma<- cov(theta)
  corr<- cor(theta)
  # corr<- matrix(c(1,.5,.5,1), ncol=2)
  # sigma<- matrix(c(.9,.8,.8,.9), ncol=2)
  
  Cholcorr<- chol(corr) 
  CholSigm<- chol(sigma)
  
  #t(chol(corr))%*%chol(corr)
  #t(chol(sigma))%*%chol(sigma)
  
  A_asterisco<- A_matrix%*%CholSigm%*%solve(Cholcorr)
  A_asterisco[is.na(aux2)] = 0
  
  #	A%*%t(theta)
  #	A_asterisco%*%Cholcorr%*%solve(CholSigm)%*%t(theta)
  
  if(verbose) print("Done")
  list(coefs = cbind(A_asterisco, coef_d, coef_c), corr = corr)
}

#Calculate init vals for multidimensional models in IRTPP
calculate_init_vals = function(dat, init_uni, dim_clust, corr, verbose = FALSE, probit = FALSE){
  #Input:
  # dat: (matrix) it's dichotomic dataset
  # init_uni: (matrix) it's initial values for coefs. 
  # dim_clust: array with clusters dimensions
  # corr: (matrix) it's intial value of correlation matrixmatrix
  # Verbose: (boolean) print message in start and finish process
  # probit: (boolean) if probit = False then the model is logit else the model is probit
  
  #Output:
  # coefs: (matrix) initial values to multidimensional IRT model
  # corr: (matrix) correlations between constructs 
  
  
  if(verbose) print("Start process to calculate multidimensional init vals.")
  
  #Control section
  if(length(unique(dim(corr))) > 1 || !isSymmetric(corr) || sum(corr > 1) > 0) stop("The matrix corr is not correlation matrix")
  if(!is.logical(probit)) stop("The value of probit is not logical")
  if(sum(dat %in% c(0,1)) + sum(is.na(dat)) != dim(dat)[1] * dim(dat)[2]) stop("The matrix dat is not binary")
  if(is.null(dim_clust) | length(dim_clust) != dim(corr)[1]) stop("The value of dim_clust is invalid")
  if(dim(init_uni)[2] - 2 != dim(corr)[2]) stop("Dimensions of init_uni matrix and corr matrix do not correspond")
  
  #Convert to probit
  if(!probit){
    sub_init = init_uni[,c(1:(ncol(init_uni) - 1))]
    sub_init = sub_init / 1.702
    init_uni[,c(1:(ncol(init_uni) - 1))] = sub_init
  }
  
  #Create Fpatt and Fval
  dim = dim(corr)[2]
  Fpatt = matrix(1, nrow = ncol(dat), ncol = dim)
  start_act = 1
  for(i in 1:length(dim_clust)){
    Fpatt[start_act,] = rep(0, dim)
    start_act = start_act + dim_clust[i]
  }
  
  Fval = init_uni[,1:dim]
  
  
  
  #Created Ppatt
  Ppatt = matrix(1, ncol = dim, nrow = dim)
  diag(Ppatt) = 0

  #Get noharm fit
  fit = sirt::noharm.sirt( dat = dat, Ppatt = Ppatt, Fpatt = Fpatt, Fval = Fval, Pval = corr)
  
  #build list for return
  if(!probit){
    coef_a = fit$loadings * 1.702
    coef_b = fit$final.constants * 1.702
  }else{
    coef_a = fit$loadings
    coef_b = fit$final.constants
  }
  
  coefs = cbind(coef_a, coef_b, init_uni[,ncol(init_uni)])
  colnames(coefs) = c(paste("a_", 1:dim, sep = ""), "d", "c")
  corr = fit$factor.cor
  colnames(corr) = c(paste("a_", 1:dim, sep = ""))
  rownames(corr) = c(paste("a_", 1:dim, sep = ""))
  ret = list(coefs = coefs, corr = corr)
  
  if(verbose) print("Done")
  
  ret
}

#Calculate init vals in two steps
init_vals_mult = function(dat, size.cluster, verbose, probit){
  suppressWarnings({
    fit_uni = calculate_init_vals_unidim(data = dat, size.cluster = size.cluster, verbose = verbose)
  })
  fit = calculate_init_vals(dat = dat, init_uni = fit_uni$coefs, dim_clust = size.cluster, corr = fit_uni$corr, verbose = verbose, probit = verbose)
  fit
}