#Descritption: data generation process with noise model 1

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
#path <- "../test_ss/ss_10"
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
k.gen <- 9
M <- 2
obs.time <- seq(0,1,1/50)


N50 =  c(275, 305, 344, 393, 458, 550, 688, 917, 1376, 2752, 13761)
N100 = c(324, 360, 405, 462, 540, 648, 810, 1080, 1620, 3240,16200)
N150 = c(352, 391, 440, 503, 587, 705, 881, 1175, 1762, 3525, 17626)

N_list = list()
N_list$N50  = list(N50, 50)
N_list$N100 = list(N100, 100)
N_list$N150 = list(N150, 150)

for( li in N_list){
    for (n in li[[1]]){
        for (p in li[[2]]){
            ## be careful for the choice of the number of basis function 
            ## fourier basis: km must be odd
            ## bspline basis km>4
            print(paste("N=",n,"p=",p))
            km.gen <- c(9,9)

            Apinv_list <- list()
            A_list <- list()
            N_list <- list()
            basis.m_list <- list()
            true.basis_list <- list()
            true.values_list <- list()


            graph.filename <- paste("graph_",cov_name,"_p", p, "_N", n,sep="")
            data.filename <- paste("data_",cov_name,"_p", p, "_N", n,sep="")

            #generate data from the graph

            if (cov_name == "power"){
                omega <- synth.omega.power(p, k.gen)
                omega <- solve(omega)
            }
            if (cov_name == "tridiag1"){
                omega <- synth.omega.tridiag1_v2(p, k.gen)
            }
            if (cov_name == "tridiag2"){
                omega <- synth.omega.tridiag2_v2(p, k.gen)
            }
            if (cov_name == "tridiag3"){
                omega <- synth.omega.tridiag3_v2(p, k.gen)
            }

            #ensure that the diagonal values are all 1

            G.true <- matrix(0, p, p) # p by p adjacency matrix
            for(i in 1:p){
            for(j in 1:p){
                if(sum(abs(omega[((i-1)*k.gen+1):(i*k.gen), ((j-1)*k.gen+1):(j*k.gen)])) > 0)
                G.true[i,j] <- 1
            }
            }

            #compute the covariance matrix
            cov <- solve(omega)



            for(m in 1:M){
                Am <- synth.linear_op.sparse_orthogonal(k.gen, km.gen[m], k.gen, scale=.2)
                Am <- t(t(Am) %*% diag(.2*((1:k.gen))+1)) #this is to make the singular values distinct
                #Am <- diag(.2*((1:k.gen)))
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
            utility.save_graphs(A_list, B_list, G.true, thre=1e-3, path, graph.filename)

            #generate data 
            data <- synth.data_from_graph(n, p, cov, basis.m_list, Apinv_list,N_list, dependent=TRUE, addnoise=FALSE)

            ########################
            #  Estimate A by CCA   #
            ########################


            #data modality 1
            fourier.basis1 <- create.fourier.basis(rangeval=c(0,1), nbasis=km.gen[1])

            d1 <-Data2fd(argvals=obs.time, y=t(data[[1]][,1,]), basisobj=fourier.basis1)
            #data modality 2
            fourier.basis2 <- create.fourier.basis(rangeval=c(0,1), nbasis=km.gen[2])

            d2 <-Data2fd(argvals=obs.time, y=t(data[[2]][,1,]), basisobj=fourier.basis2)
            cca.r.est <- estimate.cca.basis_expansion(d1, d2, km.gen[1])

            ##save data 
            est_A <- list()
            est_A[[1]] <- cca.r.est$A1
            est_A[[2]] <- cca.r.est$A2 

            print(paste('difference A1 (before callibration)', norm(est_A[[1]] - A_list[[1]])))
            print(paste('difference A2 (before callibration)', norm(est_A[[2]] - A_list[[2]])))

            #calibrate A
            for(i in 1:k.gen){
                if(norm(est_A[[1]][i,]+A_list[[1]][i,]) < norm(est_A[[1]][i,]-A_list[[1]][i,])){
                    est_A[[1]][i,] <- - est_A[[1]][i,] 
                }

            }

            for(i in 1:k.gen){
                
                if(norm(est_A[[2]][i,]+A_list[[2]][i,]) < norm(est_A[[2]][i,]-A_list[[2]][i,])){
                    est_A[[2]][i,] <- - est_A[[2]][i,] 
                }
            }

            print(paste('difference A1 (after callibration)', norm(est_A[[1]] - A_list[[1]])))
            print(paste('difference A2 (after callibration)', norm(est_A[[2]] - A_list[[2]])))
            #this is for model 2, no longer needed anymore
            est_D <- list()
            est_D[[1]] <- cca.r.est$D1
            est_D[[2]] <- cca.r.est$D2
            
            #generate functional data
            X_score <- list()
            basis_list <- list()
            basis_list[[1]] <- fourier.basis1
            basis_list[[2]] <- fourier.basis2
            for (m in 1:M){
                data.arr <- array(data = NA, dim = c(n,p,km.gen[[1]]), dimnames = NULL)
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

