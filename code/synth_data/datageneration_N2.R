#Descritption: data generation process with noise model 2

library(pracma)
library(matrixcalc)
library(fields)
library(wordspace)

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
path <- "../data20_batch/ss_1"
#cov_name <- "power"

### Reads in arguments passed in from command line
### args is a vector of strings 


args <- (commandArgs(TRUE))

for(i in 1:length(args)){
  # run each element of the vector as if passed directly to the console
  # have run.ind

  if (grepl("cov_name", args[[i]])) cov_name <- strsplit(args[[i]], split="=")[[1]][2]
 
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
eig_ratio <- 0.001
obs.time <- seq(0,1,1/50)


for (n in c(50,100, 200, 500)){
    for (p in c(50,100,150)){
        ## be careful for the choice of the number of basis function 
        ## fourier basis: km must be odd
        ## bspline basis km>4
        km.gen <- c(9,9)

        Apinv_list <- list()
        A_list <- list()
        N_list <- list()
        basis.m_list <- list()
        true.basis_list <- list()
        true.values_list <- list()


        graph.filename <- paste("graph_",cov_name, "_p", p, "_N", n, "_noise3", sep="")
        data.filename <-  paste("data_", cov_name, "_p", p, "_N", n, "_noise3", sep="")

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
            Am <- synth.linear_op.sparse_orthogonal(k.gen, km.gen[m], 4, scale=.2)
            Am <- t(t(Am) %*% diag(.2*((1:k.gen))+1)) #this is to make the singular values distinct

            A_list[[m]] <- Am
            Apinv_list[[m]] <- solve(Am)

            #############################
            # Generate structured noise #
            #############################

            #calculate the eigenvalue of the covariance
            Acov <- kronecker(diag(p), solve(Am)) %*% cov%*% t(kronecker(diag(p), solve(Am)))
            eig_A <- eigen(Acov, symmetric = TRUE)
            min_eigA <- min(eig_A$values)
            max_eigA <- max(eig_A$values)
            

            
            x<- synth.linear_op.sparse_orthogonal(p*km.gen[m], p*km.gen[m], num_block=max(p/10, 3), scale=1)
            if(m == 1){
                rev.n <- p*km.gen[m]
                x <- x[rev.n:1,]

            }   

            noise_cov <- (x+t(x))/2


            eig_y <- eigen(noise_cov, symmetric = TRUE) 
            id <- which(eig_y$values > 0)
            #print(paste("rank:",length(id)))    
            y4 <- eig_y$vectors[,id] %*% diag(eig_y$values[id]) %*% t(eig_y$vectors[,id])


            eig_y4 <- eigen(y4, symmetric=TRUE)
            eig_y4$values <- eig_y4$values / max(eig_y4$values) * (eig_ratio*(max_eigA))

            y5 = eig_y4$vectors %*% diag(eig_y4$values) %*% t(eig_y4$vectors)
            N_list[[m]]  <- y5


            basis.m_list[[m]] <- synth.fourier.bases.m(obs.time, km.gen[m])
            #compute the inverse covariance
            omega_m <- solve(Acov + N_list[[m]])

            ###print the information of the graph of modality m
            Bm_list <- utility.graph2B(omega_m,p)
            print(paste("Display graph information of modality ",m))
            utility.print.graph.info(Bm_list,p, thre=1e-2)
        }

        #convert to regression B
        B_list <- utility.graph2B(omega,p)
        # save true graphs
        utility.save_graphs(A_list, B_list, G.true, thre=1e-2, path, graph.filename)

        #generate data 
        data <- synth.data_from_graph(n, p, cov, basis.m_list, Apinv_list,N_list, dependent=TRUE, addnoise=FALSE)

        ##estimate A by CCA

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

        print(paste('difference A1 (after calibration)', norm(est_A[[1]] - A_list[[1]])))

        print(paste('difference A2 (after calibration)', norm(est_A[[2]] - A_list[[2]])))
        #this is for model 2, no longer needed anymore
        est_D <- list()
        est_D[[1]] <- cca.r.est$D1
        est_D[[2]] <- cca.r.est$D2

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



