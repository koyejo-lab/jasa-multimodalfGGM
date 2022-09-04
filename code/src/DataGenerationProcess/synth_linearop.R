#Description:
#Synthesize the linear operator A

library(pracma)

#produce idendity operator, for sanity check
synth.linear_op.identity <- function(k){
    return(diag(k))
}
#produce orthogonal matrix
synth.linear_op.orthogonal <- function(k, km, scale=1){  

    A <- matrix(rnorm(k*km), nrow=k, ncol=km)
    if(k <= km){
        A <- t(gramSchmidt(t(A))$Q)
    } else {
        A <- gramSchmidt(A)$Q
    }
    return(A*scale)
}

#produce block sparse orthogonal matrix
synth.linear_op.sparse_orthogonal <- function(k, km, num_block, scale=1, permute=FALSE){  
    
    w <- floor(km / num_block)
    h <- floor(k / num_block)
    re_w <- km - (num_block - 1)*w
    re_h <- k - (num_block - 1)*h
    A <- matrix(0, nrow=k, ncol=km)
    for( i in 1:num_block){
        if(i== num_block){
            sub_m <- matrix(rnorm(re_w*re_h), nrow=re_h, ncol=re_w)
            w_idx <- ((i-1)*w + 1):km
            h_idx <- ((i-1)*h + 1):k
        }
        else{
            sub_m <- matrix(rnorm(w*h), nrow=h, ncol=w)    
            w_idx <- ((i-1)*w + 1) : (i*w)
            h_idx <- ((i-1)*h + 1) : (i*h)
        }
        if(k <= km){
            sub_m <- t(gramSchmidt(t(sub_m))$Q)
        } else {
            sub_m <- gramSchmidt(sub_m)$Q
        }
        A[h_idx, w_idx] <- sub_m * scale
    }
    if (permute){
        per_idx <- randperm(k,k)
        A <- A[per_idx,]
    }
    return(A)
}




main <- function(){
    k <- 5
    km <- 15
    num_block <- 5
    scale <- 2.5
    A <- synth.linear_op.sparse_orthogonal(k,km,num_block,scale,permute = TRUE)
    B <- A %*% t(A)
    check <- all.equal(diag(B), scale**2*rep(1,k))
    print(paste("Pass the test of row norm: ", check))
    #verify the the off-diagonal term is zero

}

#test functions
#if (!interactive()) {
#  main()
#}