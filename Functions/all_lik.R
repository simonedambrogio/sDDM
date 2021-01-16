dsDDM <- function(data, theta, qu, dp, s, B_intercept, B_slope, h, nDT, N) {
  
  #RT ddm and RT ddm position
  xpos = seq(-5,5,length.out=1024); steps = xpos[2] - xpos[1]
  data$RTddm<-ifelse(data$choice>0, data$rts, -data$rts)
  data$RTddm_pos<-sapply(1:nrow(data), function(i){ which.min(abs(xpos - data$RTddm[i]))})
  ust=unique(data[, c('social', 'temptation')])
  
  probs = NULL
  for (i in nrow(ust)) {
    
    rts <- sddm_parallel(dp=dp, B_intercept=B_intercept, B_slope=B_slope, theta=theta, qu=qu, s=s,
                         temptation=ust$temptation[i], social_strength=ust$social[i], 
                         dt=h, nDT=nDT, N=N)
    rts = rts[rts!=0]
    xdens = density(rts, from=-5, to=5, n=1024, bw=0.11) #matches the prob. space of RTs from above
    idx = which(data$social==ust$social[i] &  data$temptation==ust$temptation[i])
    probs = steps*xdens$y[data$RTddm_pos[idx]]
  }
  
  probs[probs==0] = 1e-100
  return(probs)
}

all_lik <- function(par_chains, data, svalues, nDT=0.3, h=0.005, N=3000, cores){
  
  Subjects=unique(data$Subject)
  ret <- mclapply(Subjects, function(j){
    input_data <- na.omit(data[data$Subject==j, ]) #j is the 
    dummy <- list()
    for(cp in 1:par_chains){#for each chain
      
      dummy[[cp]] <- sum(log(dsDDM(data = input_data, 
                                   dp = svalues[[j]][[cp]][1], 
                                   theta = svalues[[j]][[cp]][2], 
                                   s = svalues[[j]][[cp]][3], 
                                   qu = svalues[[j]][[cp]][4],
                                   B_intercept = svalues[[j]][[cp]][5], 
                                   B_slope = svalues[[j]][[cp]][6], 
                                   h = h, N = N, nDT = nDT
        
      )))
      
    }
    return(dummy)
  }, mc.cores = cores)
  
  return(ret)
}
