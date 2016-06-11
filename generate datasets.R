#rm(list=ls())
########ARGUMENTOS:

#size.cluster: tamaño de cada cluster de items (dim. teorica)
#dim.data: dimension de los datos, (dim. practica)
#model: tipo de moodelo dicotomico: 2pl o 3pl
#sample.size: num de individuos

########DEVUELVE:

#1.  test 
#2.  parámetros con los que se generó, bajo la parámetrizacion
#en la que se distribuye la discriminacion, en el predictor lineal
#3. indices de cada cluster de items

#rm(list=ls())

#normal truncada y funcion que normaliza vectores:
source("C:/Users/Milder/Desktop/simdicomulti/rtnorm.r")
source("/home/pacho/Proyectos Personales/rtnorm.r")


simulate.dic.multi=function(size.cluster=c(25,25,25,25),
                            dim.data=3,sample.size=1000,model="2PL", seed_data = 1289){
  
  #tamaño del test(num. de items)
  K=sum(size.cluster)
  
#dimension teorica, numero de clusters
num.cluster=length(size.cluster)

#indices para especificar items por cada cluster
ind2=cumsum(size.cluster)
ind1=c(1,ind2[-length(ind2)]+1)

#vectores directores principales
set.seed(2)
dir.beta.1=matrix(runif(dim.data^2,min=0.08,max=0.28),
nrow = dim.data,ncol = dim.data)

if(num.cluster!=dim.data){
dir.beta.2=rbind(dir.beta.1,
matrix(rep(1,dim.data*(num.cluster-dim.data)),ncol=dim.data))
diag(dir.beta.2)=1
dir.beta=normalize(dir.beta.2)
}else{
diag(dir.beta.1)=1
dir.beta=normalize(dir.beta.1)
}


##ruido, variacion respecto a los vectores directores principales
noise <- seq(0.08, 0.25 ,length.out=num.cluster)

##betas, vectores directores de cada item
betas.1 = matrix(NA,K,dim.data)

# semilla para generar los vectores directores de cada ítem.
seed_beta <- 100L
set.seed(seed_beta)

#generando los vectores directores de cada item:
for(i in 1:length(ind1)){
betas.1[ind1[i]:ind2[i],]=
matrix(dir.beta[i,], size.cluster[i],dim.data,byrow=TRUE)+matrix(runif(size.cluster[i]*dim.data,-noise[i],noise[i]), size.cluster[i],dim.data)
}

# replace  negative componentes
betas.2=ifelse(betas.1<0,0,betas.1)
betas=normalize(betas.2)


##se fijan d vectores directores al origen de la dimension de los datos
#para resolver la falta de identificabilidad del modelo
for(i in 1:dim.data){
betas[ind1[i],i]=1 ; betas[ind1[i],-i]=0
}


##se calculan los vectores directores principales "reales"
dir_betas=matrix(NA,nrow=num.cluster,ncol=dim.data)

for(i in 1:nrow(dir_betas)){
  dir_betas[i,]=abs(eigen(t(betas[ind1[i]:ind2[i],])%*%betas[ind1[i]:ind2[i],])$vectors)[,1]
}


#######################################################################
# 2.item parameters
#######################################################################

# alpha params
# seed for reproductibility
seed_alpha <- 400L
set.seed(seed_alpha)
l_a=0.25
u_a = Inf


#se simulan las normas de los alphas
norm.alphas<-rtnorm(K, mean = 0, sd = 1.0, lower = l_a,  upper = u_a)

#se generan los alpha,(alpha=beta|alpha|)
alphas = betas * matrix(norm.alphas, K, dim.data, byrow=FALSE)


#se fijan d alphas para resolver la no identificabilidad.
for(i in 1:dim.data){
  alphas[ind1[i],i]=1 ; alphas[ind1[i],-i]=0
}


##slope parameter
seed_gamma <- 250L
set.seed(seed_gamma)
sd.gamma <-1
Gamma <- rnorm(K,0,sd.gamma)

##no identificcabilidad
for(i in 1:dim.data){
  Gamma[ind1[i]]=0 
}




# guessing parameter
seed_gamma <- 350L
set.seed(seed_gamma)
if(model=="3PL"){
c=runif(K,0,0.25)
}else{c=rep(0,K)}

for(i in 1:dim.data){
  c[ind1[i]]=0 
}




#######################################################################
# 3. latent traits
#######################################################################
# sample size

seed.5 <-500L
set.seed(seed.5)
theta <- matrix(rnorm(sample.size*dim.data,0,1),sample.size,dim.data)
# this line is to garantee the the covariance matrix is the identity
theta <- theta %*% solve((chol(cov(theta))))
theta.true <- theta


#######################################################################
# 3. test data
#######################################################################
# 3.1  linear predictor
######################################################################

eta  <- theta%*% t(alphas) -  matrix(Gamma,sample.size,K,byrow=TRUE)

######################################
#  3.2. probability matrix(logit model)
#####################################
P = c + (1-c)/(1+exp(-eta))

######################################
#  3.3 Test data
#####################################
# Uniform
set.seed(seed_data)
U <- runif(sample.size*K)
dim(P) <- sample.size*K
Y <- array(ifelse(P>U,1,0))
dim(Y) <- c(sample.size,K)
dim(P) <- c(sample.size,K)


retorno=list(data=Y,param.items=cbind(Gamma,alphas,c),indclust1=ind1,indclust2=ind2)
return(retorno)
}

