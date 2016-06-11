#Example
#rm(list = ls())

#load libraries
library(ltm)
library(FactoMineR)

#Generate dataset
size.cluster = c(25,25)
sim = simulate.dic.multi(size.cluster=size.cluster,dim.data =  length(size.cluster), model = "3PL", seed_data = 20, sample.size = 1000) 
dat = sim$data


#Calculate init vals
fit = init_vals_mult(dat = dat, size.cluster = size.cluster, verbose = TRUE, probit = FALSE)


#Compare results
coef_pob = sim$param.items
coef_pob = coef_pob[,c(2:3,1,4)]
colMeans(abs(fit$coefs - coef_pob)[,1:2])

fit$corr
fit
