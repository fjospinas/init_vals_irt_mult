#Calculate init vals for multidimensional models in IRTPP

#Input:
# dat: (matrix) it's dichotomic dataset
# init_uni: (matrix) it's initial values for coefs. 
# Pval: (matrix) it's intial value of correlation matrixmatrix
# dim: (integer) number of dimensions
# eps_cut_fpat: (numeric) If some value of a for any item in any dimensiónis < eps_cut_fpat, re calculate Fpatt and remove a of that item in that dimension 
# Verbose: (boolean) print message in start and finish process
# cut_fpat: (boolean) reload Fpatt

#Output:
# coefs: (matrix) initial values to multidimensional IRT model
# corr: (matrix) correlations between constructs 

calculate_init_vals = function(dat, init_uni, dim_clust, corr, verbose = FALSE, probit = FALSE){
  if(verbose) print("Start process to calculate init vals.")
  
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
 
  #Set lower bounds
  lower = init_uni[,ncol(init_uni)]

  if(TRUE %in% (lower < 0)){
    stop("The value of c should be >= 0")
  }
  
  #Get noharm fit
  fit = sirt::noharm.sirt( dat = dat, Ppatt = Ppatt, Fpatt = Fpatt, Fval = Fval, Pval = Pval, lower = lower)

  if(verbose) print("Done")
  
  #build list for return
  if(!probit){
    coef_a = fit$loadings * 1.702
    coef_b = fit$final.constants * 1.702
  }
  
  coefs = cbind(coef_a, coef_b, init_uni[,ncol(init_uni)])
  colnames(coefs) = c(paste("a_", 1:dim, sep = ""), "d", "c")
  corr = fit$factor.cor
  colnames(corr) = c(paste("a_", 1:dim, sep = ""))
  rownames(corr) = c(paste("a_", 1:dim, sep = ""))
  ret = list(coefs = coefs, corr = corr)
  ret
}
