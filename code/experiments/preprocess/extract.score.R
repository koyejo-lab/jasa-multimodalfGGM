library(fda)
library(R.matlab)
library(wavethresh)


############################
# extract fMRI score       #
#                          #
############################


#estimate score from fmri data
fmri <- readMat('../data/fmri/fmri_concate_run3_reduced_movie.mat')

#remove first 7 samples and bottom 7 samples
fmri.t <- fmri$ts[,,8:150]
print(dim)
#flatten the array
fmri.arr = array(fmri.t, dim = c(dim(fmri.t)[1]*dim(fmri.t)[2], 143))
obs.time <- seq(15,300,2)

#apply fourier basis
basis <- create.fourier.basis(rangeval=c(1, 300), nbasis=100)

#select basis id based on spectral power
fmri.fdobj<-Data2fd(argvals=obs.time, t(fmri.arr), basis)
fmri.rmeans <- rowMeans(sqrt(fmri.fdobj$coefs**2))
par(bg="white")
cos.m = matrix(0,50)
sin.m = matrix(0,50)

for( i in 1:50){
    cos.m[i] = fmri.rmeans[paste('cos', toString(i), sep="")]
    sin.m[i] = fmri.rmeans[paste('sin', toString(i), sep="")]
}
#set the threshold for truncation
trun.thre = 2.
idx = intersect(which(cos.m > trun.thre), which(sin.m > trun.thre))

row.idx = 'const'
for(i in idx){
    row.idx <- c(row.idx, paste('sin', i, sep=""), paste('cos', i, sep=""))
}

#extract data

#select coefficients based idx

select.coef <- matrix(0, length(row.idx), dim(fmri.fdobj$coefs)[2])
rownames(select.coef) <- row.idx
m <- 1
for(id in row.idx){
    select.coef[m,] <- fmri.fdobj$coefs[id,]
    m <- m+1
}
print(paste('number of selected basis', length(row.idx)))

#reshape the array back into original data
fmri.score <- array(t(select.coef), dim=c(dim(fmri.t)[1], dim(fmri.t)[2], length(row.idx)), dimnames=list(NULL,NULL,row.idx))



############################
#  extract eeg score       #
#                          #
############################


#estimate score from eeg data

eeg <- readMat('../data/eeg/eeg_concate_run3_down70_reduced_movie.mat')
t.frame <- dim(eeg$ts)[3]
eeg.arr <- array(eeg$ts, dim=c(dim(eeg$ts)[1]*dim(eeg$ts)[2], t.frame))


#perform basis selection from wavelet basis family

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
                      c("Coiflets",5))

baseSelect <- data.frame(WaveletFam = character(),Entropy=numeric(),stringsAsFactors = FALSE)

for(i in waveletFamily){
  coefs <- array(0, dim=t.frame)
  for(j in 1:dim(eeg.arr)[1]){
    waveletDecomp <- wd(family = i[1],data = eeg.arr[j,], filter.number = i[2])
    nthresh <- waveletDecomp$nlevel-1
    s_id <- 1

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
  EntropyB <- Shannon.entropy(abs(coefs)/max(abs(coefs))) / t.frame
  baseSelect<-rbind(baseSelect,data.frame(waveletFam = paste(i[1],i[2]),Entropy = EntropyB,stringsAsFactors = FALSE))
}

#select the basis family that has smallest entropy
select.family <- strsplit(min(baseSelect[,1]), " ")[[1]]
print(paste("select basis family:", select.family[1], "filter number:", as.integer(select.family[2])))


coefs <- array(0, dim=t.frame)
for(j in 1:dim(eeg.arr)[1]){
    waveletDecomp <- wd(family = "Coiflets", data = eeg.arr[j,], filter.number = 1)
    nthresh <- waveletDecomp$nlevel-1
    s_id <- 1
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
coefs <- coefs / dim(eeg.arr)[1]
trun.thre.eeg <- 0.7
select.idx <- which(coefs > trun.thre.eeg)
print(paste("number of selected basis:",length(select.idx)))

#extract score 
eeg.score <- array(0, dim = c(dim(eeg.arr)[1], length(select.idx)))


for(j in 1:dim(eeg.arr)[1]){
    coef <- array(0, dim=t.frame)
    waveletDecomp <- wd(family = "Coiflets", data = eeg.arr[j,], filter.number = 1)
    nthresh <- waveletDecomp$nlevel-1
    s_id <- 1
    for (k in 0:nthresh) {
        coefd <- sqrt(accessD(waveletDecomp, level = k)**2)
        end_id = length(coefd) + s_id
        coef[s_id:(end_id-1)] <-  coefd
        s_id <- end_id
        
    }
    coefc <- sqrt(accessC(waveletDecomp, level = 0)**2)
    end_id <- length(coefc) + s_id
    coef[s_id:(end_id-1)] <-  coefc
    eeg.score[j,] <- coef[select.idx]
}

#reshape the data back to original shape
eeg.score2 <- array(eeg.score, dim=c(dim(eeg$ts)[1], dim(eeg$ts)[2],length(select.idx)))

data <- list(fmri=fmri.score, eeg=eeg.score2)
save(data, file='../data/resting_data/fmri_eeg_score_movie_run3.Rdata')
print('score saved!')
