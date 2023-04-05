library(pracma)
library(matrixcalc)
library(fields)
library(wordspace)
library(Matrix)
#load source file
src.path <- "../src"

source(paste(src.path, "DataGenerationProcess", "synth_basis.R", sep="/"))
source(paste(src.path, "DataGenerationProcess", "synth_data.R", sep="/"))
source(paste(src.path, "DataGenerationProcess", "synth_graph.R", sep="/"))
source(paste(src.path, "DataGenerationProcess", "synth_linearop.R", sep="/"))

source(paste(src.path, "Estimation", "basis_estimation.R", sep="/"))
source(paste(src.path, "Utility", "utility.R", sep="/"))
source(paste(src.path, "Estimation", "cca_estimation.R", sep="/"))
source(paste(src.path, "Utility", "R2python.R", sep="/"))

###specify save path and filename
#path <- "../data21"
#cov_name <- "power"

### Reads in arguments passed in from command line
### args is a vector of strings 


args <- (commandArgs(TRUE))

for(i in 1:length(args)){
  # run each element of the vector as if passed directly to the console
  # have run.ind

  if(grepl("cov_name", args[[i]])) cov_name <- strsplit(args[[i]], split="=")[[1]][2]
  else if(grepl("path", args[[i]])) path <- strsplit(args[[i]], split="=")[[1]][2]
  else eval(parse(text = args[[i]]))
}

if(!exists("cov_name")){
    cov_name <- "tridiag1"
}
if(!exists("path")){
    path <- "../test_ss/ss_2"
}

###specify parameters
n <- 200
p <- 150

M <- 2
obs.time <- seq(0,1,1/50)
k.gen <- 9


for (km in c(5,7,9,11)){
for (n in c(100, 500)){
    for (p in c(50,100,150)){
        ## be careful for the choice of the number of basis function 
## fourier basis: km must be odd
## bspline basis km>4
print(paste("N=",n,"p=",p))
km.gen <- c(k.gen,k.gen)

Apinv_list <- list()
A_list <- list()
N_list <- list()
basis.m_list <- list()
true.basis_list <- list()
true.values_list <- list()

for(k in seq(1,km,by=2)){
graph.filename <- paste("graph_",cov_name, "_p", p, "_N", n, "estkm_", km,'est_k', k, sep="")
data.filename <-  paste("data_", cov_name, "_p", p, "_N", n, "estkm_", km,'est_k', k, sep="")

#generate data from the graph


if (cov_name == "power"){
    omega <- synth.omega.power(p, k.gen)
    omega <- solve(omega)
}
if (cov_name == "tridiag1"){
    omega <- synth.omega.tridiag1(p, k.gen)
}
if (cov_name == "tridiag2"){
    omega <- synth.omega.tridiag2(p, k.gen)
}
if (cov_name == "tridiag3"){
    omega <- synth.omega.tridiag3(p, k.gen)
}
if (cov_name == "tridiag4"){
    omega <- synth.omega.tridiag4(p, k.gen)

}

#ensure that the diagonal values are all 1

G.true <- matrix(0, p, p) # p by p adjacency matrix
for(i in 1:p){
  for(j in 1:p){
    if(sum(abs(omega[((i-1)*k.gen+1):(i*k.gen), ((j-1)*k.gen+1):(j*k.gen)])) > 0)
      G.true[i,j] <- 1
  }
}

cov <- solve(omega)


for(m in 1:M){
    Am <- synth.linear_op.sparse_orthogonal(k.gen, km.gen[m], k.gen, scale=.2)
    Am <- t(t(Am) %*% diag(.2*((1:k.gen))+1)) #this is to make the singular values distinct

    A_list[[m]] <- Am
    Apinv_list[[m]] <- solve(Am)
}

for(m in 1:M){
    
    #noise covariance
    N_list[[m]] <- diag(p*km.gen[m])*.05
    basis.m_list[[m]] <- synth.fourier.bases.m(obs.time, km.gen[m])
}

#convert to regression B
B_list <- utility.graph2B(omega,p)
# save true graphs
utlity.save_graphs(A_list, B_list, G.true, thre=1e-3, path, graph.filename)

#generate data 
data <- synth.data_from_graph(n, p, cov, basis.m_list, Apinv_list,N_list, dependent=TRUE, addnoise=FALSE)

##estimate A by CCA

#data modality 1

concate_data <- list()
for(m in 1:M){
        dm <- data[[m]][,1,]
    for(i in 2:p){
        dm <- rbind(dm, data[[m]][,p,])
    }
    concate_data[[m]] <- dm
}

#implement basis number selection
#sum_list <- list()
#for(m in 1:M){
#    for( b in c(3, 5, 7, 9, 11, 13, 15, 17)){
#        basis <- create.fourier.basis(rangeval=c(0, 1), nbasis=b)
#        sum.a <- 0
#        g.fdobj<-Data2fd(argvals=obs.time, t(concate_data[[m]]), basis)

#        for(i in 1:(dim(concate_data[[m]])[1])){
#            a <- sum(sqrt((eval.fd(obs.time, g.fdobj[i])- concate_data[[m]][i,])**2))/sum(concate_data[[m]][i,]**2)
#            sum.a <- sum.a+a
#        }
#        sum_list <- append(sum_list,sum.a)
#    }
#    print(sum_list)
#}





fourier.basis1 <- create.fourier.basis(rangeval=c(0,1), nbasis=km)
d1 <-Data2fd(argvals=obs.time, y=t(data[[1]][,1,]), basisobj=fourier.basis1)
#data modality 2
fourier.basis2 <- create.fourier.basis(rangeval=c(0,1), nbasis=km)
d2 <-Data2fd(argvals=obs.time, y=t(data[[2]][,1,]), basisobj=fourier.basis2)
cca.r.est <- estimate.cca.basis_expansion(d1, d2, k)

##save data 
est_A <- list()
est_A[[1]] <- cca.r.est$A1
est_A[[2]] <- cca.r.est$A2 

print(paste('difference A1 (before callibration)', norm(est_A[[1]] - A_list[[1]][1:k,1:km])))
print(paste('difference A2 (before callibration)', norm(est_A[[2]] - A_list[[2]][1:k,1:km])))

for(i in 1:k){
    if(norm(est_A[[1]][i,]+A_list[[1]][i,]) < norm(est_A[[1]][i,]-A_list[[1]][i,1:km])){
        est_A[[1]][i,] <- - est_A[[1]][i,] 
    }

}

for(i in 1:k){
    
    if(norm(est_A[[2]][i,]+A_list[[2]][i,]) < norm(est_A[[2]][i,]-A_list[[2]][i,1:km])){
        est_A[[2]][i,] <- - est_A[[2]][i,] 
    }
}

print(paste('difference A1', norm(est_A[[1]] - A_list[[1]][1:k,1:km])))

print(paste('difference A2', norm(est_A[[2]] - A_list[[2]][1:k,1:km])))
#this is for model 2, no longer needed anymore
est_D <- list()
est_D[[1]] <- cca.r.est$D1
est_D[[2]] <- cca.r.est$D2

X_score <- list()
basis_list <- list()
basis_list[[1]] <- fourier.basis1
basis_list[[2]] <- fourier.basis2
for (m in 1:M){
    data.arr <- array(data = NA, dim = c(n,p,km), dimnames = NULL)
    for(i in 1:p){
        obj <- Data2fd(argvals=obs.time, y=t(data[[m]][,i,]), basisobj=basis_list[[m]] )
        data.arr[,i,] <- t(obj$coefs)
    }
    X_score[[m]] <- data.arr
}

utility.save_data(est_A, est_D, X_score, data, path, data.filename)
    }
}

}

}