
#functional basis
library(fda)


#estimate the basis via functional pca
estimate.fpca.basis <- function(rangeval, obs.time, obs.data, nbasis, bs.nbasis=10){
    #obs.data: observation data, dimension  time*n
    #requirement: nbasis <= bs.nbasis
    bspline.basis <- create.bspline.basis(rangeval=rangeval, nbasis=bs.nbasis)
    fd.object <- Data2fd(argvals=obs.time, y=obs.data, basisobj=bspline.basis)
    pca.result <- pca.fd(fd.object, nharm=nbasis)
    return(list("vectors"=pca.result$harmonics, "values"=pca.result$values))
}

#compute the functional score
estimate.fscore <- function(rangeval, obs.time, obs.data, basis.m, bs.nbasis=10){
    #return matrix n*km
    #convert to functional object
    bspline.basis <- create.bspline.basis(rangeval=rangeval, nbasis=bs.nbasis)
    fd.object <- Data2fd(argvals=obs.time, y=obs.data, basisobj=bspline.basis)
    #compute functional score
    f.score <- inprod(fd.object, basis.m)
    return(f.score)
}

#estimate the basis via functional pca with data aggragation accross nodes
estimate.fpca.sum.basis <- function(rangeval, obs.time, data, num.basis, bs.nbasis=10){
    #sum of covariance operator 
    # data n*p*t
    if(length(dim(data)) > 2){
        p <- dim(data)[2]
    }
    else{
        p <- 1
    }
    
    cov.sum <- matrix(0, bs.nbasis, bs.nbasis)
    fourier.basis <- create.fourier.basis(rangeval=rangeval, nbasis = bs.nbasis)
    if(p == 1){
        fd.object <- Data2fd(argvals=obs.time, y=t(data), basisobj = fourier.basis)
        cov.sum <-  var.fd(fd.object)$coefs
    } else {
        for(i in 1:p){
            fd.object <- Data2fd(argvals=obs.time, y=t(data[,i,]), basisobj = fourier.basis)
            cov.sum <-  cov.sum +  var.fd(fd.object)$coefs
        }
        cov.sum <- cov.sum / p
    }
    
    cov.eigvector <- eigen(cov.sum, symmetric=TRUE)$vectors
    cov.eigenval  <- eigen(cov.sum, symmetric=TRUE)$values
    
    true_basis <- t(cov.eigvector)  %*% t(eval.basis(obs.time, fourier.basis))
    true_basis <- true_basis[1:num.basis,]
    fd.object <- Data2fd(argvals=obs.time, y=t(true_basis), basisobj=fourier.basis)
    
    r.list <- list("vectors" = fd.object, "values" = cov.eigenval, "cov" = cov.sum[1:num.basis,1:num.basis])
    return(r.list)


}