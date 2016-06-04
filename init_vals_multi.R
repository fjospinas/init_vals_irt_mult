#Calculate init vals for multidimensional models in IRTPP

#Input:
# dat: (matrix) it's dichotomic dataset
# init_uni: (matrix) it's initial values for coefs. 
# Pval: (matrix) it's intial value of correlation matrixmatrix
# dim: (integer) number of dimensions
# eps_c: (numeric) epsilon to move c
# eps_cut_fpat: (numeric) If some value of a for any item in any dimensiónis < eps_cut_fpat, re calculate Fpatt and remove a of that item in that dimension 
# Verbose: (boolean) print message in start and finish process
# cut_fpat: (boolean) reload Fpatt

#Output:
# coefs: (matrix) initial values to multidimensional IRT model
# corr: (matrix) correlations between constructs 

calculate_init_vals = function(dat, init_uni, Pval, dim, eps_c = 0, eps_cut_fpat = NULL, verbose = FALSE, cut_fpat = FALSE){
  if(verbose) print("Start process to calculate init vals.")
  #Create Fpatt and Fval
  Fpatt = matrix(1, nrow = ncol(data), ncol = dim)
  Fval = matrix(0, nrow = ncol(data), ncol = dim)
  for(i in 1:nrow(init_uni)){
    Fval[i,init_uni[i,1]] = init_uni[i,2] 
  }
  
  #Created Ppatt
  Ppatt = matrix(1, ncol = dim, nrow = dim)
  diag(Ppatt) = 0
  
  #Set lower and upper bounds
  lower = init_uni[,ncol(init_uni)] - eps_c
  upper = rep(1,nrow(init_uni))
  
  #Get noharm fit
  fit = noharm.sirt( dat = dat, Ppatt = Ppatt, Fpatt = Fpatt, Fval = Fval, Pval = Pval, lower = lower, upper = upper)
  
  #if cut_pat == TRUE AND values < eps_cut_pat reload Fpatt and calcule fit 
  if(cut_fpat){
    if(is.null(eps_cut_fpat) || eps_cut_fpat <= 0 || eps_cut_fpat > 0.5){
      stop("Invalid value for: eps_cut_fpat")
    }
    Fpatt = ifelse(abs(fit$loadings.theta) < eps_cut_fpat,0,1)
    equal_zero = rowSums(Fpatt) == 0
    Fpatt[equal_zero,] = rep(1,dim)
    
    #If fit fail for "reorder" return previous fit
    fit = tryCatch({
      fit = noharm.sirt( dat = dat, Ppatt = Ppatt, Fpatt = Fpatt, Fval = Fval, Pval = Pval, lower = lower, upper = upper)
    }, error = function(e) {
      fit
    })
    
  }

  if(verbose) print("Done")
  
  #build list for return
  coefs = cbind(fit$loadings, fit$final.constants, init_uni[,ncol(init_uni)])
  colnames(coefs) = c(paste("a_", 1:dim), "d", "c")
  corr = fit$factor.cor
  colnames(corr) = c(paste("a_", 1:dim))
  rownames(corr) = c(paste("a_", 1:dim))
  ret = list(coefs = coefs, corr = corr)
  ret
}



