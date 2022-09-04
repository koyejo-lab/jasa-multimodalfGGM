library(R.matlab)
utility.save_data <- function(A,D, X_score,raw_data, path, filename){
    # X_score is a list and each with dimension (n,p, km)
    # A is a list of estimated A matrix, each with dimension (k,km)
    source <- list()
    source$M <- length(X_score)
    source$p <- dim(X_score[[1]])[2]
    source$data <- X_score
    source$A <- A
    source$D <- D 
    source$raw <- raw_data
    save(source, file=paste(path,"/",filename,".RData",sep=""))
    print("Finish saving data.")

}

utlity.save_graphs <- function(A,B, G.true, thre=1e-3, path,  filename){
    param <- list()
    param$M <- length(A)
    param$p <- length(B)
    param$G <- G.true
    s <- 0
    k <- as.integer(dim(B[[1]])[2] / (param$p-1))
    print(k)
    alpha <- 0
    Bij_norm <- matrix(0, length(B), p-1)
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


    param$alpha <- alpha
    param$sparsity <- s
    param$Am <- A
    param$B <- B
    print(paste("sparsity:", s))
    print(paste("alpha:", alpha))
    save(param, file=paste(path,"/",filename,".RData",sep=""))
    print("Finish saving graphs")
}

utlity.load_estimation <- function(path, filename){
    obj.list = readMat(paste(path,filename,sep="/"))
    print("Finish loading data")
    print("Data names:")
    print(names(obj.list))
    
    return(obj.list)
}