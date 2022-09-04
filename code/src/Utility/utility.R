
### compute mean
utility.mean_data <- function(data_list){
    M <- length(data_list)
    mean_list <- list()
    for(m in 1:M){
        p <- dim(data_list[[m]])[2]
        #n,p.tau
        mean_data <- data_list[[m]][,1,]*(1/p)
        for(i in 2:p){
            mean_data <-  mean_data + data_list[[m]][,i,]*(1/p)
        }
        mean_list[[m]] <- mean_data
    }
    return(mean_list)
}

### graph to coefficient
utility.graph2B <- function(omega, p){
    #omega pk time pk
    k <- as.integer(dim(omega)[1] / p)
    B <- list()
    for(i in 1:p){
        idx <- ((i-1)*k + 1) : (i*k)
        omega_i <- 1./diag(omega[idx,idx])
        yidx <- ((i-1)*k + 1) : (i*k)

        if(i != 1 & i != p){
            xidx <- c(1 : ((i-1)*k), ((i*k) + 1) : (k*p))
        } else if (i ==1){
            xidx <-  ((i*k) + 1) : (k*p)
        } else{
            xidx <-  1 : ((i-1) * k)
        }
        
    
        Bi <- -diag(omega_i) %*% omega[yidx, xidx]  
        B[[i]] <- Bi
    }
    return(B)
}
### store data
utilty.store <- function(){

}

### select number of basis

utility.select_k <- function(corr, threshold=.5){
    return (which(abs(corr) > threshold))
}

utility.select_km <- function(){
    #TBD
}


utility.remove.errormsg <- function(){
    
}

utility.print.graph.info <- function(B,p, thre){
    Bij_norm <- matrix(0, length(B), p-1)
    k <- as.integer(dim(B[[1]])[2] / (p-1))
    alpha <- 0
    s <- 0
    for(i in 1:length(B)){
        s0 <- 0
        for(j in 1:(p-1)){

            Bij <- B[[i]][, ((j-1)*k+1):(j*k)]
            Bij_norm[i,j] <- norm(Bij,"f")
            if(norm(Bij, "f") > thre){
                s0 <- s0 + 1
                max_col <- max(colSums(abs(Bij) > thre))
                max_row <- max(rowSums(abs(Bij) > thre))
                alpha <- max(alpha, max_col, max_row)
            }
        }
        s <- max(s, s0)
    }
    Bij_norm_arr <- as.vector(Bij_norm)
    Bij_norm_arr <- sort(Bij_norm_arr, decreasing=TRUE)
    print(table(findInterval(Bij_norm_arr, c(1e-5,1e-4,1e-3,1e-2,1e-1))))


    print(paste("sparsity:", s))
    print(paste("alpha:", alpha))

}