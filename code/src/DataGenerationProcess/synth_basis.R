#Description:
#Synthesize bases



library(fda)
### generate fourier basis
# 0.3. Observed basis function values given a vector of time points
synth.fourier.bases.m <- function(obs.time, k, rangeval=c(0, 1)){
  # Input:
  #   obs.time, the observation grid of length tau
  #   k, number of basis functions
  #   rangeval, time interval, a vector of length two
  # Output:
  #   tau * k matrix of observed function values
  fourier <- create.fourier.basis(rangeval=rangeval, nbasis=k)

  fourier.m <- eval.basis(obs.time, fourier) 
  return(fourier.m)
}

### generate bspline basis
synth.bspline.bases.m <- function(obs.time, k, rangeval=c(0, 1)){
    # Input:
    #   obs.time, the observation grid of length tau
    #   k, number of basis functions
    #   rangeval, time interval, a vector of length two
    # Output:
    #   tau * k matrix of observed function values
    bspline <- create.bspline.basis(rangeval=rangeval, nbasis=k, norder=4)
    bspline.m <- eval.basis(obs.time, bspline)
    return(bspline.m)
}

### gnerate linear baplise bases
synth.linear.bases.m <- function(obs.time, k, rangeval=c(0, 1)){
    # Input:
    #   obs.time, the observation grid of length tau
    #   k, number of basis functions
    #   rangeval, time interval, a vector of length two
    # Output:
    #   tau * k matrix of observed function values
    linear <- create.bspline.basis(rangeval=rangeval, nbasis=k, norder=2)
    linear.m <- eval.basis(obs.time, linear)
    return(linear.m)
}

### generate exponential basis
synth.exponential.bases.m <- function(obs.time, k, rangeval=c(0, 1)){
    # Input:
    #   obs.time, the observation grid of length tau
    #   k, number of basis functions
    #   rangeval, time interval, a vector of length two
    # Output:
    #   tau * k matrix of observed function values
    exponential <-create.exponential.basis(rangeval=rangeval, nbasis=k)
    exponential.m <- eval.basis(obs.time, exponential)
    return(exponential.m)
}

