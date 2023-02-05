

estimate.cca <- function(obj1,obj2, basis1, basis2, lambda1, lambda2, ncan=NULL){
    ###################
    # obj1: functional object1
    # obj2: functional object1
    # basis1: funcational basis object 1
    # basis2: functional basis object 2
    ###################
    
    D.matrix <- estimate.cca.matrix(obj1,obj2, basis1, basis2, lambda1, lambda2, ncan)

    svd.D <- svd(D.matrix)

    if(is.null(ncan)){
        ncan <- length(length(which(svd.D$d > 1e-8)))
    }
    #print(ncan)
    inv_score1 <-  diag(1 / sqrt(lambda1))
    inv_score2 <- diag(1 / sqrt(lambda2))
    corr <- svd.D$d[1:ncan]
    A1.m <- t(svd.D$u[,1:ncan] %*% inv_score1[1:ncan, 1:ncan])
    A2.m <- t(svd.D$v[,1:ncan] %*% inv_score2[1:ncan, 1:ncan])
    #A1.m <- t(svd.D$u[,1:ncan] )
    #A2.m <- t(svd.D$v[,1:ncan] )
    r.list <- list("corr"=corr, "A1"=A1.m, "A2"=A2.m)
    return(r.list)
}


estimate.cca.basis_expansion <- function(obj1, obj2, ncan){
    #obj, obj2:fd object
    stopifnot(dim(obj1$coefs)[2] == dim(obj2$coefs)[2])
    n <- dim(obj1$coefs)[2]
    d1 <- obj1$coefs#t(inprod(obj1, obj1$basis))
    d2 <- obj2$coefs#t(inprod(obj2, obj2$basis))

    cov1  <- (d1 %*% t(d1)) / n
    cov2  <- (d2 %*% t(d2)) / n
   
    cov12 <- (d1 %*% t(d2)) / n

    
    inv_cov1.s <- solve( cov1 + 1e-8*diag(dim(cov1)[1]) )

    #inv_cov1.s <- solve( cov1 )
    est.r <- eigen(inv_cov1.s , symmetric=TRUE)
    
    inv_cov1.s <- est.r$vectors %*% diag(sqrt(est.r$values)) %*% t(est.r$vectors)



    inv_cov2.s <- solve( cov2 + 1e-8*diag(dim(cov2)[1]) )

    est.r <- eigen(inv_cov2.s,  symmetric=TRUE)
    inv_cov2.s <- est.r$vectors %*% diag(sqrt(est.r$values)) %*% t(est.r$vectors)
    
    
    m <- (inv_cov1.s %*% cov12) %*% inv_cov2.s
    m <- m %*% t(m)
    
    eigen.m <- eigen(m, symmetric=TRUE)
    print("canonical correlation values:")
    print(eigen.m$values)
    inv_corr <- diag(1. / sqrt(eigen.m$values[1:ncan]+1e-8))
 
    A1 <- inv_cov1.s  %*% eigen.m$vectors[,1:ncan] %*% inv_corr
    
    A2 <- inv_cov2.s %*% inv_cov2.s %*% t(cov12) %*% A1 

    return(list("A1"=t(A1), "A2"=t(A2), "cca"=eigen.m$values))

}


estimate.cca.matrix <- function(obj1,obj2, basis1, basis2, lambda1, lambda2, ncan=NULL){
    ###################
    # obj1: functional object1
    # obj2: functional object1
    # basis1: funcational basis object 1
    # basis2: functional basis object 2
    ###################
    

    score1 <- inprod(obj1, basis1)
    n <- dim(score1)[1]
    inv_score1 <-  diag(1 / sqrt(lambda1))
    score1 <- score1 %*% inv_score1

    score2 <- inprod(obj2, basis2)
    inv_score2 <- diag(1 / sqrt(lambda2))
    score2 <- score2 %*% inv_score2

    D.matrix <- (t(score1) %*% score2)/n
    return(D.matrix)

}

estimate.multivariate.cca <- function(X,Y,ncan){
    # data X: n*p1
    # data Y: n*p2
    stopifnot(dim(X)[1] == dim(Y)[1])
    n <- dim(X)[1]
    XY.cov <- (t(X) %*% Y) / n
    YX.cov <- (t(Y) %*% X) / n
    X.cov <- (t(X) %*% X) / n
    Y.cov <- (t(Y) %*% Y) / n
    X.inv.cov <- solve(X.cov)

    est.r <- eigen(X.inv.cov, symmetric=TRUE)
    X.inv.cov.s <- est.r$vectors %*% diag(sqrt(est.r$values)) %*% t(est.r$vectors)



    Y.inv.cov <- solve(Y.cov)

    est.r <- eigen(Y.inv.cov,  symmetric=TRUE)
    Y.inv.cov.s <- est.r$vectors %*% diag(sqrt(est.r$values)) %*% t(est.r$vectors)
    
    
    m <- (X.inv.cov.s %*% XY.cov) %*% Y.inv.cov.s
    m <- m %*% t(m)
    
    eigen.m <- eigen(m, symmetric=TRUE)
    print("canonical correlation values:")
    print(eigen.m$values)
    inv_corr <- diag(1. / sqrt(eigen.m$values[1:ncan]))
 
    A1 <- X.inv.cov.s  %*% eigen.m$vectors[,1:ncan] %*% inv_corr
    
    A2 <- Y.inv.cov.s %*% Y.inv.cov.s %*% YX.cov %*% A1 

    return(list("A1"=t(A1), "A2"=t(A2), "cca"=eigen.m$values))

}