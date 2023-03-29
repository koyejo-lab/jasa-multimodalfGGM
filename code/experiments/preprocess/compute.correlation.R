library(pracma)
library(matrixcalc)
library(fields)

#library(wordspace)

#load source file
src.path <- "../../src"
source(paste(src.path, "Utility", "utility.R", sep="/"))
source(paste(src.path, "Estimation", "cca_estimation.R", sep="/"))
source(paste(src.path, "Utility", "R2python.R", sep="/"))

load('../data/resting_data/fmri_eeg_score_movie_run2.Rdata')
fmri <- data$fmri
print(paste("dimension of fmri data", dim(fmri)[1], dim(fmri)[2], dim(fmri)[3]))
eeg <- data$eeg
print(paste("dimension of eeg data", dim(eeg)[1], dim(eeg)[2], dim(eeg)[3]))


estimate.A <- function(obj1, obj2, ncan){
    #obj, obj2: km by N
    stopifnot(dim(obj1)[2] == dim(obj2)[2])
    n <- dim(obj1)[2]
    d1 <- obj1
    d2 <- obj2

    cov1  <- (d1 %*% t(d1)) / n
    cov2  <- (d2 %*% t(d2)) / n
   
    cov12 <- (d1 %*% t(d2)) / n

    
    inv_cov1.s <- solve( cov1 + 1e-8*eye(dim(cov1)[1]))

    #inv_cov1.s <- solve( cov1 )
    est.r <- eigen(inv_cov1.s , symmetric=TRUE)
    
    inv_cov1.s <- est.r$vectors %*% diag(sqrt(est.r$values)) %*% t(est.r$vectors)



    inv_cov2.s <- solve( cov2 + 1e-8*eye(dim(cov2)[1]) )

    est.r <- eigen(inv_cov2.s,  symmetric=TRUE)
    inv_cov2.s <- est.r$vectors %*% diag(sqrt(est.r$values)) %*% t(est.r$vectors)
    
    
    m <- (inv_cov1.s %*% cov12) %*% inv_cov2.s
    m <- m %*% t(m)
    
    eigen.m <- eigen(m, symmetric=TRUE)
    print("canonical correlation values:")
    print(eigen.m$values)
    inv_corr <- diag(1. / sqrt(eigen.m$values[1:ncan]))
 
    A1 <- inv_cov1.s  %*% eigen.m$vectors[,1:ncan] %*% inv_corr
    
    A2 <- inv_cov2.s %*% inv_cov2.s %*% t(cov12) %*% A1 

    return(list("A1"=t(A1), "A2"=t(A2)))

}

k <- min(dim(fmri)[3], dim(eeg)[3])


#aggregate sample
A1 <- list()
A2 <- list()
for(i in 1:dim(eeg)[2]){
    x <- estimate.A(t(fmri[,i,]), t(eeg[,i,]),k)
#    par(mfrow=c(1,2),bg="white")
#    print(dim(x$A1))
#    print(dim(x$A2))
    A1[[i]] <- x$A1
    A2[[i]] <- x$A2
#    image.plot(x$A1, col=viridis(128))
#    image.plot(x$A2, col=viridis(128))
}


fmri_arr <- matrix(0, dim(fmri)[1]*dim(fmri)[2], dim(fmri)[3])
eeg_arr <- matrix(0, dim(eeg)[1]*dim(eeg)[2], dim(eeg)[3])
for(i in 1:dim(fmri)[1]){
    for(j in 1:dim(fmri)[2]){
        fmri_arr[(i-1)*dim(fmri)[2]+j,] <- fmri[i,j,]
        eeg_arr[(i-1)*dim(fmri)[2]+j,] <- eeg[i,j,]
    }
}

x <- estimate.A(t(fmri_arr), t(eeg_arr),k)

A1[[dim(eeg)[2] + 1]] <- x$A1
A2[[dim(eeg)[2] + 1]] <- x$A2

source <- list()
source$A1 <- A1
source$A2 <- A2

save(source, file="../data/resting_data/Amatrix_movie_run2.Rdata")