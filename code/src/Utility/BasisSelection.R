
utility.select.k <- function(corr_list, thre=1e-3){
    selectk$num <- length(which(corr_list >= thre))
    selectk$idx <- which(corr_list >= thre)
    return(selectk)
}

utility.select.km <- function(basis_score, thre=1e-3){
    ####
    #input dimension N, p, km
    ####
    bs_s <- sqrt(basis_score**2)
    bs_mean <- colMeans(colMeans(bs_s))

    selectkm$num <- length(which(bs_mean >= thre))
    selectkm$idx <- which(bs_mean >= thre)

    return(selectkm)

}