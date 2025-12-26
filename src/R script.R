# Title        : Non- parametric method to find MLE method for 1 marker 
# Objective    : Contains implementation of the model (EM-algorithm) and supporting functions
# Created by   : Loyce Kayanula and Kristan. A. Schneider
# Created on   : 31.01.23
# Last modified: 31.12.23


# load or install if needed required pakages 
package_names <- c("caret", "mltools", "data.table", "dplyr",  "purrr", "tidyr", "stringr", "openxlsx")

# Check if each package is installed, install if missing, then load
for (package in package_names) {
    if (!requireNamespace(package, quietly = TRUE)) {
      install.packages(package, dependencies = TRUE)
    }
    library(package, character.only = TRUE)
}

#***********************************************************************************************
#Function for importing data
#***********************************************************************************************


DatImp <- function(path){
  if(substring(path,nchar(path)-3,nchar(path))==".xls"){
    dat <- openxlsx::read.xlsx(path,1)
  }
  else{
    if(substring(path,nchar(path)-4,nchar(path))==".xlsx"){
      dat <- openxlsx::read.xlsx(path,1)
    }
    else{
      if(substring(path,nchar(path)-3,nchar(path))==".txt"){
        dat <- read.table(path,header=TRUE, sep="\t")
      }
      else{
        if(substring(path,nchar(path)-3,nchar(path))==".csv"){
          dat <- read.csv(path,header=TRUE,sep=";")
        }
      }
    }  
  }
  dat
}  


#***********************************************************************************************
#This function calculates the ML estimates for MOI and the haplotype frequency estimates
#***********************************************************************************************

Obsv <- function(n) { 
  if(n == 1) {
    return (matrix(c(1), nrow = 1)) 
  }
  h <- array(0, c(2^n, n))
  h[1:2, 1] <- c(0, 1)
  for (i in 2:n){
    h[(2^(i - 1) + 1):2^i, 1:(i - 1)] <- h[1:2^(i - 1), 1:(i - 1)] 
    h[(2^(i-1)+1):2^i, i] <- 1
  }
  h <- h[2:2^n, 1:n]
  h
}

ax <- function(sampli,n){
  # return a listed list
  # 1st sublist ... returns sets Ax for all x with y Ax represented as vectors of numerics
  # 2nd sublist ... matches first list and returns signs (-1)^(|x|-|y|) 
  # 3rd list is 0-1 matrix of samples y that will be needed
  # 4th ... sampli 
  # smapli .... vector of samples as numerics
  lst <-list(NULL)
  lst1 <-  list(NULL)

  ms <- do.call(rbind,(lapply(sampli,function(x) as.integer(unlist(strsplit(x,""))))))
  ms <- data.frame(ms)
  
  rownames(ms) <- paste(sampli)
  for (i in 1:nrow(ms)) { 
    z <- c(1:n)[ms[i,]==1]
    absx <- sum(ms[i,])
    Obsv1 <- Obsv(absx)
    H <- Obsv1 %*% (2^(z-1))
    
    obs <- array(0,c(nrow(H),n))
    obs[,z] <- Obsv1
    nam <- apply(obs,1,function(x) paste(x,collapse=""))
    rownames(H) <- nam
    lst[[i]] <- H
    lst1[[i]] <- (-1)^(absx - rowSums(Obsv1))

  } 
  ob1 <- unique(unlist(lapply(lst,rownames)))
  ob1 <- do.call(rbind,(lapply(ob1,function(x) as.integer(unlist(strsplit(x,""))))))
  rownames(ob1) <-  paste(ob1 %*% (2^(0:(n-1))))

  list(lst,lst1,ob1,sampli)  
}

#powers for the multinomial

Mpws <- function(x, M) {
  cumprod (rep(x, M))               
}

Mpws1 <- function(x, M) {
  c(1, cumprod(rep(x, M - 1)))        
}

## PGF of MOI 
Gkap <-  function(kap, t) {         
  # kap .... prob mass function of m 
  # t ... argument
  M <- length(kap)
  t(kap%*%sapply(t,Mpws,M))
}

## 1st derivative of PGF of m
dGkap <-  function(kap,t){  
  M <- length(kap)
  t(((1:M)*kap)%*%sapply(t,Mpws1,M))
}


## This is the main function
NonparametricEst <- function (dat, M, p=NULL, km = NULL, maxit = 1000,eps = 10^(-8) ){
  # dat...    data in the following form: it is a list with the fist element containing a 0-1 matrix with the moluecular marker in one-hot encoding
  #          with different samples corresponding to rows and colums to aleles. The second list element is a vector indicating how often the respective observations occur in the data
  # M ...    maximum MOI
  # p.....   initial haplotype frequencies with default set to NULL
  # km ...   initial MOI distrbution with default set to NULL
  # maxit... maximum iterations for the algorithm
  # eps ...  numerical threshold for convergence of algorithm
  
  
  if(is.null(km)){
    km <- dpois(1:M, 1)  
    km <- km / sum(km)
  }
  
  nx <- dat[[2]]              # number of times sample is Obsverved
  
  N <- sum(nx)
  Niorig <- nx %*% dat[[1]]    # the number of times each allele occurs in the samples
  
  pick.1 <- Niorig > 0          ## selects only alleles actually present
  
  n <- sum(pick.1)
  
  if(is.null(p)){
    p <- rep(1/n,n)
  }else{
    p <- c(p[pick.1])
    p <- p/sum(p)
  }   
  
  norig <- ncol(dat[[1]])
  
  Ni <- Niorig[pick.1]
  
  datnew <- dat[[1]][,pick.1]
  
  rownames(datnew) <-   apply(datnew,1,function(x) paste(x,collapse=""))
  names(nx) <- rownames(datnew)
  sampli <- rownames(datnew) # c(datnew %*% (2)^(0:(n-1)))
  

  Axlist <- ax(sampli,n)                     
  nsam <- length(sampli)                      
  
  Ptx <- array(,c(nsam,1))                    
  Rmx <- array(,c(nsam,M))                   
  Tx <-  array(,c(nsam,n))                   
  
  if(n==1){  # if just 1 lineage 
    warning("Only one lineage observed, all MOI distributions are equally likely, estimation not meaningful")
    p <- 1
    km <- rep(NA,M)
    
  }else{    
    ## no degenerate data for EM algorithm
    
    delta <- 1                             
    t <- 0      
    
    while (delta > eps && t < maxit) {
      t <- t +1
      #print("test1")
      Pvec <- Axlist[[3]]%*% p  
      #print("test")
      rownames(Pvec) <- rownames(Axlist[[3]])      
      Gsump <- Gkap(km,Pvec)                        
      rownames(Gsump) <- rownames(Pvec)           
      
      dGsump <- dGkap(km,Pvec)                    
      rownames(dGsump) <- rownames(Pvec)          
      
      sump.mpw <- sapply(Pvec,Mpws,M)
      
      colnames(sump.mpw ) <- rownames(Pvec) 
      
      for(i in 1:length(Axlist[[1]])){                
        pick <- paste(Axlist[[1]][[i]])                              
        Ptx[i] <- Gsump[pick,]%*%Axlist[[2]][[i]]                     
        Tx[i,] <- Axlist[[2]][[i]]%*%(Axlist[[3]][pick,]*dGsump[pick,])      
        Rmx[i,] <- sump.mpw[,pick]%*%as.matrix(Axlist[[2]][[i]])            
        
      }
      
      cPtx <- c(Ptx)                                                
      
      p1 <- rowSums(t(Tx*nx/cPtx)*p)                           
      p1 <- p1/sum(p1)                             
      
      Rmt <- colSums(Rmx/cPtx*nx)
      Rmt <- Rmt*km
      km1 <- Rmt/sum(Rmt)                            
      
      delta <- sqrt(sum((p-p1)^2)+(sum((km-km1)^2))) 
      
      p <- p1
      km <- km1
      #print(km)
      #print(delta)
    }
    if (t == maxit) {
      warning(paste("EM algorithm failed to converge within", maxit, "iterations "))
    }  
    
    
  }
  
  pout <- array(0,norig)
  pout[pick.1] <- p
  names(pout) <- colnames(Niorig) 
  
  Niout <- array(0,norig)
  Niout[pick.1] <- Ni
  meanMOI <- km %*% (1:M) 
  names(km) <- paste("P[MOI=",1:M,"]",sep="")
  
  obs_prev <- Niout/N
  names(obs_prev) <- colnames(Niorig) 
  est_prev <- c(1-Gkap(km,1-pout))
  est_prev[!pick.1] <- 0
  names(est_prev) <- colnames(Niorig) 
  
  out <- list(pout, km, c(meanMOI), M, Niorig, N, obs_prev, est_prev, 1-sum(pout^2))
  names(out) <- c( "MLE of lineage freqs.", "MLE of MOI distribution", "average MOI", "maximum MOI for algorithm", "lineage counts", "sample size", "observed prevalence","estimated prevalences", "heterozygosity")
  
  out
  
}


## User friendly wrap around main function
MOI.NP <- function(data, markername, M,  CI=FALSE, B=1000, alpha_level = 0.95, p = NULL, km=NULL, maxit =1000, eps=10^-8){
  # data ..........  input data
  # markername ....  name of the markers present in the data
  # M .............  maximum MOI for the algorithm
  # CI ............  if TRUE, Bootstrap confidence intervals will be perfomed 
  # B .............  number of bootsrap repeats used for CIs - ignored if CI = FALSE
  # alpha_level ...  coverage of bootsrap CIs - ignored if CI = FALSE
  # p .............  initial frequency distribution with default set to NULL
  # km ............  initial  MOI distribution with default set to NULL
  # maxit .........  maximum iterations for the algorithm
  # eps ...........  numerical threshold for convergence of algorithm

  # first transform data into one-hot encoding, i.e., each sample
  # corresponds to a row of 0 or 1 indicating the lineages present
  data <- data.frame(SampleID = factor(data[,1]), marker = factor(data[,markername]))
  colnames(data) <- c("SampleID", markername)
  Dt <- as.data.frame(one_hot(data.table(data), col = markername, dropUnusedLevels = TRUE, sparsifyNAs = TRUE))    # Perform one-hot encoding using one_hot
  Dt <- aggregate(. ~ SampleID , data = Dt, FUN = max)
  
  pattern <- paste0("^", markername, "_")  # the ^tells that the string has to start with the given characters
    
  sel_markr <- colnames(Dt)[str_detect(colnames(Dt), pattern)]
  marker <- as.data.table (Dt %>% select(SampleID, all_of(sel_markr)))

  #all unique infections per marker  and their counts 
  unique_samps <-  marker[, .(COUNT = .N), by = setdiff(names(marker), c("SampleID"))]
    
  # excluding all samples with no infection i.e, rows with zero
  unique_samps <- unique_samps %>% filter(rowSums(select(., -COUNT)) > 0)
    
  lst1 <- unique_samps[, setdiff(names(unique_samps), "COUNT"), with = FALSE]
  lst2 <- unique_samps[,"COUNT"]
  N <- sum(lst2) 
  newdt <- list(as.matrix(lst1), lst2[[1]]) 
  
  estimates.NonP <- NonparametricEst(newdt, M, p, km,  maxit=maxit, eps=eps )
 
  out <- list()
  out$estimates <- estimates.NonP

  
  if(CI){
    
    p <- estimates.NonP[[1]]
    km <- estimates.NonP[[2]]
    dataB <- data.table(t(rmultinom(B,N, unlist(lst2/N))))
    
    unique_samps <-  dataB[, .(COUNT = .N), by = c(colnames(dataB))]
    Bweight <- c(unique_samps[,"COUNT"])
    unique_samps1 <- unique_samps[,-"COUNT"]
    B1 <- nrow(unique_samps1)
    newdtB <- newdt
    boot.rep <- list()
    for(b in 1:B1){
      
      newdtB[[2]] <- unlist(unique_samps1[b,])
            
      sel <- newdtB[[2]] > 0
      newdtB[[1]] <- newdt[[1]][sel,]
      newdtB[[2]] <- newdtB[[2]][sel]
      
      boot.rep[[b]] <- NonparametricEst(newdtB, M, p, km,  maxit=maxit, eps=eps )
      
    }
    
    ### Compile bootstrap CIs

    CI_lev <- c((1-alpha_level)/2,1-(1-alpha_level)/2)
    
    BCIs_NP <- list()
    
    BCIs_NP$"CI average MOI" <- quantile(rep(sapply(boot.rep,function(x) x$`average`),each = Bweight),CI_lev,na.rm=TRUE)
    
    BCIs_NP$"CI heterozygosity" <- quantile(rep(sapply(boot.rep,function(x) x$`heterozygosity`),each = Bweight),CI_lev,na.rm=TRUE)
    
    BCIs_NP$"CIs MOI distribution" <-  t(apply(sapply(boot.rep,function(x) x$`MLE of MOI distribution`),1,function(y) quantile(rep(y, each =Bweight),CI_lev ,na.rm=TRUE)))
    
    BCIs_NP$"CIs lineage frequencies" <-   t(apply(sapply(boot.rep,function(x) x$`MLE of lineage freqs.`),1,function(y) quantile(rep(y, each =Bweight),CI_lev,na.rm=TRUE)))
    
    BCIs_NP$"prevalences" <-  t(apply(sapply(boot.rep,function(x) x$`estimated prevalences`),1,function(y) quantile(rep(y, each =Bweight),CI_lev,na.rm=TRUE)))
    
    out$CIs <- BCIs_NP  
  }
  
  out
}
