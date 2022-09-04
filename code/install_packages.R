ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    
    if (length(new.pkg)){
        #chooseCRANmirror()
        install.packages(new.pkg, dependencies = TRUE,repos='http://cran.us.r-project.org')
    } 
        
    sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("doParallel", 
              "RSpectra",
              "fda", 
              "mvtnorm", 
              "fgm", 
              "poweRlaw", 
              "matrixcalc",
              "far",
              "R.matlab",
              "pracma",
              "plotly",
              "fields",
              "viridis",
              "wordspace"
              )
ipak(packages)

#the "fgm packages" is for fgm: Partial Separability and Functional Gaussian Graphical Models