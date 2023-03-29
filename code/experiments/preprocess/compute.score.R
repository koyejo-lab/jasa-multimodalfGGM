library(fda)
library(R.matlab)
library(wavethresh)

#########################################
# step 1 project fMRI data to the basis #
#########################################

f <- readMat('../data/fmri/fmri_concate_all_reduced_movie.mat')
#remove head and tail where signals are noisy

f_t <- f$ts[,,8:150]
g = array(f_t, dim = c(dim(f_t)[1]*dim(f_t)[2], 143))
obs.time <- seq(15,300,2)

#Uncomment the following lines to select optimal number of basis
#implement basis selection
sum_list <- list()
for( b in c(5,10,20,30,40,50,60,70,80,90,100,150,200,300)){
    basis <- create.fourier.basis(rangeval=c(1, 300), nbasis=b)
    sum.a <- 0
    g.fdobj<-Data2fd(argvals=obs.time,t(g), basis)
    for(i in 1:(dim(g)[1])){
        a <- sum(sqrt((eval.fd(obs.time, g.fdobj[i])- g[i,])**2))/sum(sqrt(g[i,]**2))
        sum.a <- sum.a+a
    }
    sum_list <- append(sum_list,sum.a/(68*26*286))
}
#print the result and select the optimal based on elbo

print(unlist(sum_list))

opt_nbasis <- 100
basis <- create.fourier.basis(rangeval=c(1, 300), nbasis=opt_nbasis)
g.fdobj<-Data2fd(argvals=obs.time,t(g), basis)
f_s <- rowMeans(sqrt(g.fdobj$coefs**2))
cos_list = matrix(0,50)
sin_list = matrix(0,50)

for( i in 1:50){
    cos_list[i] = f_s[paste('cos', toString(i), sep="")]
    sin_list[i] = f_s[paste('sin', toString(i), sep="")]
}
#png("fmri_score.png")     #Create a png file
#plot(seq(1:50),rev(sort(sqrt(cos_list**2+sin_list**2)))) 
#dev.off()
#select k_m
trun_thre = 2. #select the threshold based on elbo
idx = intersect(which(cos_list > trun_thre), which(sin_list > trun_thre))
#idx = which(sqrt(cos_list**2+sin_list**2) > trun_thre)

row_idx = 'const'
for(i in idx){
    row_idx <- c(row_idx, paste('sin', i, sep=""), paste('cos', i, sep=""))
}
select_coef <- matrix(0, length(row_idx), dim(g.fdobj$coefs)[2])
#print(paste('select number of basis',length(row_idx)))

fmri_data <- array(t(select_coef), dim=c(dim(f_t)[1], dim(f_t)[2],length(row_idx)), dimnames=list(NULL,NULL,row_idx)) 
print('dimension of fmri data:')
print(dim(fmri_data))


f <- readMat('../data/eeg/eeg_concate_all_down70_reduced_movie.mat')

f_arr <- array(f$ts, dim=c(dim(f$ts)[1]*dim(f$ts)[2], 2048))


waveletFamily <- list(c("DaubExPhase",1),
                      c("DaubExPhase",2),
                      c("DaubExPhase",3),
                      c("DaubExPhase",4),
                      c("DaubExPhase",5),
                      c("DaubExPhase",6),
                      c("DaubExPhase",7),
                      c("DaubExPhase",8),
                      c("DaubExPhase",9),
                      c("DaubExPhase",10),
                      c("Coiflets",1),
                      c("Coiflets",2),
                      c("Coiflets",3),
                      c("Coiflets",4),
                      c("Coiflets",5),
                      c("DaubLeAsymm",4),
                      c("DaubLeAsymm",5),
                      c("DaubLeAsymm",6),
                      c("DaubLeAsymm",7),
                      c("DaubLeAsymm",8),
                      c("DaubLeAsymm",9),
                      c("DaubLeAsymm",10),                                                                                                                                    
                      c("Lawton",3),
                      c("LinaMayrand",5.4))

#Uncomment the following lines to select optimal number of basis and wavelet family

baseSelect <- data.frame(WaveletFam = character(),Entropy=numeric(),stringsAsFactors = FALSE)
for(i in waveletFamily){
  coefs <- array(0, dim=2048)
  for(j in 1:dim(f_arr)[1]){
    waveletDecomp = wd(family = i[1],data = f_arr[j,], filter.number = i[2])
    nthresh = waveletDecomp$nlevel-1
    s_id = 1

    for (k in 0:nthresh) {
        coefd <- abs(accessD(waveletDecomp, level = k))
        end_id = length(coefd) + s_id
        coefs[s_id:(end_id-1)] <- coefs[s_id:(end_id-1)] + coefd
        #coefs <- list(list(sqrt(accessD(waveletDecomp, level = k)**2)),coefs)
        s_id <- end_id
    }
    coefc <- abs(accessC(waveletDecomp, level = 0))
    end_id <- length(coefc) + s_id
    coefs[s_id:(end_id-1)] <- coefs[s_id:(end_id-1)] + coefc
  }
  coefs <- coefs/(26*68) #take average of the coefficient
  EntropyB <- Shannon.entropy(abs(coefs)/max(abs(coefs)))/2048
  baseSelect<-rbind(baseSelect,data.frame(waveletFam = paste(i[1],i[2]),Entropy = EntropyB,stringsAsFactors = FALSE))
}
print('Optimal basis family')
print(strsplit(min(baseSelect[,1]), " ")[[1]][1])
print(as.integer(strsplit(min(baseSelect[,1]), " ")[[1]][2]))


coefs <- array(0, dim=2048)
for(j in 1:dim(f_arr)[1]){
    waveletDecomp = wd(family = 'Coiflets',data = f_arr[j,], filter.number = 1)
    nthresh = waveletDecomp$nlevel-1
    s_id = 1
    for (k in 0:nthresh) {
        coefd <- sqrt(accessD(waveletDecomp, level = k)**2)
        end_id = length(coefd) + s_id
        coefs[s_id:(end_id-1)] <- coefs[s_id:(end_id-1)] + coefd
        #coefs <- list(list(sqrt(accessD(waveletDecomp, level = k)**2)),coefs)
        s_id <- end_id
        
    }
    coefc <- sqrt(accessC(waveletDecomp, level = 0)**2)
    end_id <- length(coefc) + s_id
    coefs[s_id:(end_id-1)] <- coefs[s_id:(end_id-1)] + coefc
}

coefs <- coefs/(dim(f_arr)[1])
#print(length(which(coefs>0.75)))
#print(length(which(coefs>0.8)))
#select k_m
idx <- which(coefs>0.78)

#extract score 
score_arr = array(0, dim = c(dim(f_arr)[1], length(idx)))

for(j in 1:dim(f_arr)[1]){
    coef <- array(0, dim=2048)
    waveletDecomp = wd(family = 'Coiflets',data = f_arr[j,], filter.number = 1)
    nthresh = waveletDecomp$nlevel-1
    s_id = 1
    for (k in 0:nthresh) {
        coefd <- sqrt(accessD(waveletDecomp, level = k)**2)
        end_id = length(coefd) + s_id
        coef[s_id:(end_id-1)] <-  coefd
        #coefs <- list(list(sqrt(accessD(waveletDecomp, level = k)**2)),coefs)
        s_id <- end_id
        
    }
    coefc <- sqrt(accessC(waveletDecomp, level = 0)**2)
    end_id <- length(coefc) + s_id
    coef[s_id:(end_id-1)] <-  coefc
    score_arr[j,] <- coef[idx]
}
#select k_m

print(dim(score_arr))
