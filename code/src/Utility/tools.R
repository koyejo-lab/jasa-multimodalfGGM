library(bazar)

#check that the eigenvalues are distinct
tools.unique_engv <- function(cov){
    sym <- isSymmetric(cov)
    values <- eigen(cov, symmetric = sym)$values
    valid_idx <- which(abs(values) > 1e-5)
    num <- length(valid_idx)
    values <- values[valid_idx]
    unique.num <- length(almost.unique(values, tol=1e-5))
    return(num == unique.num)
}

tools.sample_cov <-function(data){
    #n*d
    return((t(data) %*% data) / n)
}
#