#Description:
#Generate observed data from latent graph, function basis, and operator A


library(mvtnorm)

### dependent data generation process
### data of different modality are dependent, where the underlying latent variable z is shared

synth.dependent_data <- function(n, p, cov, Apinv_list, N_list=NULL){ 
    # Input:
    #   n, number of samples, scaler
    #   p, number of nodes, scaler
    #   cov, covariance of the latent variables, (kp)x(kp) matrix
    #   Apinv_list, list of the inverse of Am, list of length M
    #   N_list, list of the covariance of the noise, list of length M  
    # Output:
    #   data, list of observed data
    kp <- dim(cov)[1]
    stopifnot(kp%%p == 0)

    k <- floor(kp/p)
    
    #generate latent samples
    sample_z <- rmvnorm(n, rep(0, kp), cov) # n by kp
    
    #generate observed data
    data_list <- list()
    M <- length(Apinv_list)
    for(m in 1:M){
        data_m <- matrix(, nrow = n, ncol = 0)
        # project data form latent space to observed space
        for(i in 1:p){
            idx <- ((i-1)*k + 1): (i*k)
            sub_m <- sample_z[, idx] %*% t(Apinv_list[[m]])
            data_m <- cbind(data_m, sub_m)
        }
        
        if(!is.null(N_list[[m]]) & (norm(N_list[[m]])>1e-5)){
            kmp <- dim(N_list[[m]])[1]
            #generate independent noise
            add_noise <- rmvnorm(n, rep(0, kmp), N_list[[m]])
            data_list[[m]] <- data_m + add_noise
        }
        else{
            data_list[[m]] <- data_m
        }

    }
    return(data_list)
}

###independent data generation process
###data of different modality are generated independently

synth.independent_data <- function(n, p, cov, Apinv_list, N_list=NULL){ 
    # Input:
    #   n, number of samples, scaler
    #   p, number of nodes, scaler
    #   cov, covariance of the latent variables, (kp)x(kp) matrix
    #   Apinv_list, list of the inverse of Am, list of length M
    #   N_list, list of the covariance of the noise, list of length M  
    # Output:
    #   data, list of observed data
    kp <- dim(cov)[1]
    stopifnot(kp%%p == 0)

    k <- floor(kp/p)
     # n by kp
    data_list <- list()
    M <- length(Apinv_list)
    for(m in 1:M){
        #generate latent samples
        sample_z <- rmvnorm(n, rep(0, kp), cov)
        data_m <- matrix(, nrow = n, ncol = 0)
        # project data from latent space to the observed space
        for(i in 1:p){
            idx <- ((i-1)*k + 1): (i*k)
            sub_m <- sample_z[, idx] %*% t(Apinv_list[[m]])
            data_m <- cbind(data_m, sub_m)
        }

        if(!is.null(N_list[[m]]) & norm(N_list[[m]])>1e-5){
            kmp <- dim(N_list[[m]])[1]
            #add independent noise
            add_noise <- rmvnorm(n, rep(0, kmp), N_list[[m]])
            data_list[[m]] <- data_m + add_noise
        }
        else{
            data_list[[m]] <- data_m
        }
    }
    return(data_list)
}

### generate functional data from the functional score and function basis
synth.functional_data <- function(p, score, basis.m, addnoise=TRUE){
    # Input:
    #   p, number of nodes, scaler
    #   score, functional score, n x (kmp) matrix
    #   basis.m, list of functional basis, tau x (km) matrix
    #   addnoise, add noise or not, Boolean
    # Output:
    #   data, matrix of functional data, n x p x tau matrix
    n <- dim(score)[1]
    km <- floor( dim(score)[2] / p)

    tau <- dim(basis.m)[1]
    h <- array(0, c(n, p, tau))
    for(i in 1:n){
        for(j in 1:p){
            h[i,j,] <- basis.m %*% matrix(score[i, ((j-1)*km+1) : (j*km)], ncol=1) 
            if(addnoise) h[i,j,] <- h[i,j,] + rnorm(tau, 0, 0.5)
        }
    }
    return(h)

}

### the complete process
synth.data_from_graph <- function(n, p, cov, basis.m_list, Apinv_list, N_list, dependent=FALSE, addnoise=FALSE){
    # Input:
    #   n, number of samples, scaler
    #   p, number of nodes, scaler
    #   cov, covariance of the latent variables, (kp)x(kp) matrix
    #   Apinv_list, list of the inverse of Am, list of length M
    #   N_list, list of the covariance of the noise, list of length M  
    #   basis.m_list, list of functional basis, tau x (km) matrix
    #   addnoise, add noise or not, Boolean
    # Output:
    #   data, list of observed functional data    
    if(dependent){
        data_list <- synth.dependent_data(n, p, cov, Apinv_list, N_list)
    } else {
        data_list <- synth.independent_data(n, p, cov, Apinv_list, N_list)
    }
    M <- length(data_list)
    obs <- list()
    for(m in 1:M){
        print(paste("Generating data from modality",m))
        obs[[m]] <- synth.functional_data(p, data_list[[m]], basis.m_list[[m]], addnoise=addnoise)
    }
    return(obs)
}

#synthesize true basis via functional pca
synth.true.basis <- function(obs.time, cov, pinvAm, basis.m){
    # Input:
    #   obs.time, observed time sequence
    #   cov, covariance of the latent variables, (kp)x(kp) matrix
    #   pinvAm, inverse of Am, kmxk matrix
    #   basis.m, functional basis, tau x (km) matrix
    # Output:
    #   r.list, dictionary of pca basis and associated eigenvalue
    k <- dim(pinvAm)[2]
    km <- dim(pinvAm)[1]
    p <- floor(dim(cov)[1] / k)
    sum.cov <- matrix(0,km, km)
    
    for(i in 1:p){
        for(j in 1:p){
            xidx <- (k*(i-1) + 1):(k*i)
            yidx <- (k*(j-1) + 1):(k*j) 
            sum.cov <- sum.cov + (pinvAm %*% cov[xidx, yidx] %*% t(pinvAm))
        }
    }
    sum.cov <- sum.cov / (p^2)
    eigen.result <- eigen(sum.cov, symmetric=TRUE)
    eigenvec <- eigen.result$vectors
    eigenval <- eigen.result$values
    true_basis <- t(eigenvec)  %*% t(basis.m)
    
    
    bspline.basis <- create.bspline.basis(rangeval=c(0,1), nbasis=km*3)
    fd.object <- Data2fd(argvals=obs.time, y=t(true_basis), basisobj=bspline.basis)
    
    r.list <- list("vectors" = fd.object, "values" = eigenval)
    return(r.list)
}

#synthesize true basis via functional pca with node aggregation
synth.true.sum.basis <- function(obs.time, cov, pinvAm, basis.m){

    # Input:
    #   obs.time, observed time sequence
    #   cov, covariance of the latent variables, (kp)x(kp) matrix
    #   pinvAm, inverse of Am, kmxk matrix
    #   basis.m, functional basis, tau x (km) matrix
    # Output:
    #   r.list, dictionary of pca basis and associated eigenvalue
    k <- dim(pinvAm)[2]
    km <- dim(pinvAm)[1]
    p <- floor(dim(cov)[1] / k)
    sum.cov <- matrix(0,km, km)
    if(p == 1){
        sum.cov <- pinvAm %*% cov %*% t(pinvAm)
    } else {
        for(i in 1:p){
            xidx <- (k*(i-1) + 1):(k*i)
            yidx <- (k*(i-1) + 1):(k*i) 
            sum.cov <- sum.cov + (pinvAm %*% cov[xidx, yidx] %*% t(pinvAm))
        }
        sum.cov <- sum.cov / p
    }
    
    eigen.result <- eigen(sum.cov, symmetric=TRUE)
    eigenvec <- eigen.result$vectors
    eigenval <- eigen.result$values
    true_basis <- t(eigenvec)  %*% t(basis.m)
    
    
    bspline.basis <- create.bspline.basis(rangeval=c(0,1), nbasis=km*3)
    fd.object <- Data2fd(argvals=obs.time, y=t(true_basis), basisobj=bspline.basis)
    
    r.list <- list("vectors" = fd.object, "values" = eigenval)
    return(r.list)
}

main <- function(){
    n <- 100
    p <- 50
    k <- 5
    cov <- diag(p*k)
    M <- 4
    pinvA_list <- list()
    N_list <- list()
    km_list <- c(1,3,5,7)
    for(m in 1:M){
        pinvA <- matrix(rnorm(k*km_list[m]), ncol=k)
        N <- runif(1) * diag(p*km_list[m])
        pinvA_list[[m]] <- pinvA
        N_list[[m]] <- N
    }
    dep_z <- synth.dependent_data(n, p, cov, pinvA_list, N_list)
    
    check1 <- dim(dep_z[[1]])[2] == p*km_list[1]
    check2 <- dim(dep_z[[2]])[2] == p*km_list[2]
    check3 <- dim(dep_z[[3]])[2] == p*km_list[3]
    check4 <- dim(dep_z[[4]])[2] == p*km_list[4]
    check <- check1 & check2 & check3 & check4
    print(paste("check dimension: ", check))

    indep_z <- synth.independent_data(n, p, cov, pinvA_list, N_list)
    
    check1 <- dim(indep_z[[1]])[2] == p*km_list[1]
    check2 <- dim(indep_z[[2]])[2] == p*km_list[2]
    check3 <- dim(indep_z[[3]])[2] == p*km_list[3]
    check4 <- dim(indep_z[[4]])[2] == p*km_list[4]
    check <- check1 & check2 & check3 & check4
    print(paste("check dimension: ", check))

    #check synthesis data
}



#test functions
#if (!interactive()) {
#  main()
#}