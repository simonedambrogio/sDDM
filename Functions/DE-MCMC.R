hdemcmc<-function(data, nDT, z, N, dt, fixations,
                  chain_length, burnin, thining, par_chains, cores=4 #The more the better
){
  
  # get number of individuals
  Subjects<-length(unique(data$Subject))
  USubjects<-unique(data$Subject) #uniques in the other script
  
  #Set boundaries and minimal sd
  dummy=rbind(bound_dp, bound_theta, bound_s, bound_qu, bound_B_intercept, bound_B_slope)
  lb<-dummy[,1]
  ub<-dummy[,2]
  
  #define tuning parameters
  tuning_size<-2.38/sqrt(2*length(lb))
  tuning_size_parent=2.38/sqrt(2*length(2))
  
  #create some empty lists which we fill later
  lik_cur<-chain_sets<-lik_prop<-cur<-propl<-svalues<-rep(list(list()),Subjects)
  current_parent_average<-current_parent_sigma<-proposed_parent_average<-proposed_parent_sigma<-parent_average_chain<-parent_sigma_chain<-list()
  
  #initialize parent chains
  set_range<-(ub-lb)*0.25# prevent starting values to be too close to upper or lower bound
  for(cp in 1:par_chains){
    parent_average_chain[[cp]]<-matrix(NA,chain_length,length(lb)) # chain_length is 2000
    parent_sigma_chain[[cp]]<-matrix(NA,chain_length,length(lb))
    parent_average_chain[[cp]][1,]<-runif(length(lb),min=c(lb)+set_range,max=c(ub)-set_range)
    parent_sigma_chain[[cp]][1,]<-runif(length(lb),max=ub,min=min_sd)
  }
  
  #initialize individual chains depending on model version
  for (j in 1:Subjects){
    for(cp in 1:par_chains){
      svalues[[j]][[cp]]<-runif(length(lb),min=c(lb)+set_range,max=c(ub)-set_range)#starting values independent of prior
      chain_sets[[j]][[cp]]<-matrix(NA,chain_length,length(lb))
      chain_sets[[j]][[cp]][1,]<-svalues[[j]][[cp]]
    }}
  
  #calculate likelihood of individual parameters
  lik_cur=all_lik(par_chains, data, svalues, cores = cores)
  
  keep=matrix(TRUE,Subjects,par_chains)#increases speed: if proposals are rejected we keep old likelihood calculation  
  sumlogl<-rep(0,Subjects) #for the DIC calculation
  
  for (chain_i in c(2:chain_length)){
    
    #get the currect markov states 
    #individuals
    for (j in 1:Subjects){
      for(cp in 1:par_chains){ #for all parallel chains
        cur[[j]][[cp]]<- chain_sets[[j]][[cp]][chain_i-1,]
      }} 
    #parents
    for(cp in 1:par_chains){ #for all parallel chains
      current_parent_average[[cp]]<- parent_average_chain[[cp]][chain_i-1,]
      current_parent_sigma[[cp]]<- parent_sigma_chain[[cp]][chain_i-1,]
    }
    
    #get new proposals for the parent destribution by using the differential evolution algorithm 
    for(cp in 1:par_chains){ #for all parallel chains
      doitagain<-TRUE
      ##########################################
      count=0
      while(doitagain){ #while new proposal are outside of bounds: try again
        count=count+1
        index=sample(c(1:par_chains)[c(1:par_chains)!=cp],2)#sample 2 alternative chains
        proposed_parent_average[[cp]]<-current_parent_average[[cp]]+(current_parent_average[[index[1]]]-current_parent_average[[index[2]]])*runif(1, 0, tuning_size_parent)+runif(length(ub),min=-0.001,max=0.001)
        proposed_parent_sigma[[cp]]<-current_parent_sigma[[cp]]+ (current_parent_sigma[[index[1]]]-current_parent_sigma[[index[2]]])*runif(1, 0, tuning_size_parent)+runif(length(ub),min=-0.001,max=0.001)
        
        sigma_full_log=(proposed_parent_sigma[[cp]]>=10 | proposed_parent_sigma[[cp]]<=min_sd)
        avg_full_log=(proposed_parent_average[[cp]]<=lb | proposed_parent_sigma[[cp]]<=0.01 | proposed_parent_average[[cp]]>=ub)
        
        sigma_log=any(sigma_full_log)
        avg_log=any(avg_full_log)
        doitagain<-any(avg_log | sigma_log)
        if(count>1000){
          if(sigma_log){
            for(i in 1:sum(sigma_full_log)) proposed_parent_sigma[[cp]][sigma_full_log][i] = runif(1,
                                                                                                 max=ub[sigma_full_log][i],
                                                                                                 min=min_sd)} 
          if(avg_log){
            for( in in 1:sum(avg_full_log)){proposed_parent_average[[cp]][avg_full_log][i] = runif(1,
                                                                                                 min=lb[avg_full_log][i]+set_range[avg_full_log][i],
                                                                                                 max=ub[avg_full_log][i]-set_range[avg_full_log][i])} 
          
          #proposed_parent_sigma[[cp]][proposed_parent_sigma[[cp]]<=min_sd+0.001]=1;
          doitagain=F
          } #if after 100 attempts not working reset sigma (doesn't happen anymore but while loops are dangerous without stop condition )
      }
    }
    
    
    # parent crossover (swap chain stains with small likelihood)
    for(cp in 1:par_chains){ #for all parallel chains
      if (runif(1)<0.05){
        index=sample(c(1:par_chains)[c(1:par_chains)!=cp],1)#crossover step between 2 chains
        dummy=proposed_parent_average[[cp]]
        proposed_parent_average[[cp]]=proposed_parent_average[[index]]
        proposed_parent_average[[index]]=dummy
      }
      if (runif(1)<0.05){
        index=sample(c(1:par_chains)[c(1:par_chains)!=cp],1)#crossover step between 2 chains
        dummy=proposed_parent_sigma[[cp]]
        proposed_parent_sigma[[cp]]=proposed_parent_sigma[[index]]
        proposed_parent_sigma[[index]]=dummy
      } 
    }
    
    
    #accept or withdraw proposals for the parent distribution
    #independently for each parameter using the accepted states of t-1
    for(cp in 1:par_chains){ #for all parallel chains
      for (p in 1:length(ub)){
        
        likelihood=cur_dummy=prop_dummy=0
        
        cur_dummy<-cur_dummy+(sum(dnorm(cur[[j]][[cp]][p],parent_average_chain[[cp]][chain_i-1,p],parent_sigma_chain[[cp]][chain_i-1,p],log=T)))
        prop_dummy<-prop_dummy+(sum(dnorm(cur[[j]][[cp]][p],proposed_parent_average[[cp]][p],proposed_parent_sigma[[cp]][p],log=T)))
        
        parent_propl=prop_dummy
        parent_cur=cur_dummy
        parent_ratio <- exp(parent_propl-parent_cur)
        if(parent_ratio>runif(1)){parent_average_chain[[cp]][chain_i,p]<-proposed_parent_average[[cp]][p];parent_sigma_chain[[cp]][chain_i,p]<-proposed_parent_sigma[[cp]][p]
        }else {  parent_average_chain[[cp]][chain_i,p]<-current_parent_average[[cp]][p];parent_sigma_chain[[cp]][chain_i,p]<-current_parent_sigma[[cp]][p]  }
      }}
    
    
    #get new propsoals for individual states
    for (j in 1:Subjects){
      for(cp in 1:par_chains){ #for all parallel chains
        doitagain<-TRUE
        while(doitagain){
          index=sample(c(1:par_chains)[c(1:par_chains)!=cp],2) #sample 2 alternative chains
          propl[[j]][[cp]]<-cur[[j]][[cp]]+(cur[[j]][[index[1]]]-cur[[j]][[index[2]]])*tuning_size+runif(length(lb),min=-0.001,max=0.001)
          doitagain<- any(propl[[j]][[cp]]<lb)|any(propl[[j]][[cp]]>ub)#any paramter out of bounds?
        }}}
    
    
    #crossover individual level  (swap chain stains with small likelihood)
    for (j in 1:Subjects){
      for(cp in 1:par_chains){
        if (runif(1)<0.05){
          chain_index<-sample(c(1:par_chains)[c(1:par_chains)!=cp],1)
          paramter_index<-sample(length(lb),1)
          dummy=propl[[j]][[cp]][paramter_index]
          propl[[j]][[cp]][paramter_index]<-propl[[j]][[chain_index]][paramter_index]
          propl[[j]][[chain_index]][paramter_index]=dummy
        }}}
    
    
    
    
    
    #calculate likelihood of new proposels
    lik_prop=all_lik(par_chains, data, propl, cores = cores)
    
    #################################
    #Accept or reject; if reject keep indcate to keep
    for (j in 1:Subjects){
      for(cp in 1:par_chains){
        
        lik_prior_cur = (sum(dnorm(cur[[j]][[cp]], parent_average_chain[[cp]][chain_i, ], parent_sigma_chain[[cp]][chain_i,], log=T)))
        lik_prior_prop = (sum(dnorm(propl[[j]][[cp]], parent_average_chain[[cp]][chain_i, ], parent_sigma_chain[[cp]][chain_i,], log=T)))
        
        ratio<-exp((lik_prop[[j]][[cp]]+lik_prior_prop)-(lik_cur[[j]][[cp]]+lik_prior_cur))
        if(runif(1)<ratio){
          chain_sets[[j]][[cp]][chain_i,] <- propl[[j]][[cp]];
          keep[j,cp] <- FALSE;
          sumlogl[j] <- sumlogl[j]+lik_prop[[j]][[cp]]
        }else{  
          
          chain_sets[[j]][[cp]][chain_i,] <- cur[[j]][[cp]];
          keep[j,cp] <- TRUE;
          sumlogl[j] <- sumlogl[j]+lik_cur[[j]][[cp]]
          
        }
      }
    }
    
    
    #if keep=true than the likelihood stays else the lik from the proposel
    for (j in 1:Subjects){
      for(cp in 1:par_chains){
        lik_cur[[j]][[cp]]<-ifelse(keep[j,cp]==TRUE,lik_cur[[j]][[cp]],lik_prop[[j]][[cp]])#increases speed: if proposals are rejected we keep old likelihood calculation  
      }}
    
    
    # parent migration steps where chain-states are swapped across parallel chains
    if (runif(1)<0.05){
      chain_index <- sample(par_chains,sample(2:par_chains,1))
      
      dummy_average <- parent_average_chain[[chain_index[1]]][chain_i,]
      dummy_sigma <- parent_sigma_chain[[chain_index[1]]][chain_i,]
      for(index in 1:(length(chain_index)-1)){
        parent_average_chain[[chain_index[index]]][chain_i,] <-  parent_average_chain[[chain_index[index+1]]][chain_i,]
        parent_sigma_chain[[chain_index[index]]][chain_i,] <-  parent_sigma_chain[[chain_index[index+1]]][chain_i,]
      }
      parent_average_chain[[chain_index[length(chain_index)]]][chain_i,] <- dummy_average 
      parent_sigma_chain[[chain_index[length(chain_index)]]][chain_i,] <- dummy_sigma 
    }
    
    
    # circular migration steps where chain-states are swapped across parallel chains
    for (j in 1:Subjects){
      if (runif(1)<0.05){
        chain_index <- sample(par_chains,sample(2:par_chains,1))
        dummy <- chain_sets[[j]][[chain_index[1]]][chain_i,]
        lik_dummy <- lik_cur[[j]][[chain_index[1]]] 
        for(index in 1:(length(chain_index)-1)){
          chain_sets[[j]][[chain_index[index]]][chain_i,] <-  chain_sets[[j]][[chain_index[index+1]]][chain_i,]
          lik_cur[[j]][[chain_index[index]]] <- lik_cur[[j]][[chain_index[index+1]]] 
        }
        
        chain_sets[[j]][[chain_index[length(chain_index)]]][chain_i,] <- dummy
        lik_cur[[j]][[chain_index[length(chain_index)]]] <- lik_dummy
      }}
    
    
    ##################################
    cat("Iteration: ", chain_i, " / ", chain_length,"    Count = ", count,  "\r", sep = "")#show progress 
  }
  
  #stopCluster(cl)#stop cluster
  
  #####
  #Claculating DIC
  sumlogl_hat<-rep(0,Subjects)
  #get mean posterior
  for (j in 1:Subjects){
    parmater_hat<-apply(do.call(rbind,chain_sets[[j]]),2,mean)
    
    sumlogl_hat[j]=sumlogl_hat[j]+ sum(log(dsDDM(data = data_sim, 
                                                 theta = parmater_hat[1], 
                                                 dp = parmater_hat[2], 
                                                 s = parmater_hat[3], 
                                                 qu = parmater_hat[4], 
                                                 B_intercept = parmater_hat[5], 
                                                 B_slope = parmater_hat[6], 
                                                 h = dt, nDT = nDT,
                                                 N = N)))
  }
  
  #calculate DIC
  P=2*(sumlogl_hat - (sumlogl/((chain_length-1)*par_chains)))
  DIC=sum(-2*(sumlogl_hat-P))
  
  
  ####
  #rearrange and return
  for(cp in 1:par_chains) for (j in 1:Subjects) chain_sets[[j]][[cp]]<-chain_sets[[j]][[cp]][seq(burnin+1,chain_length,thining),]
  for(cp in 1:par_chains) parent_sigma_chain[[cp]]<-parent_sigma_chain[[cp]][seq(burnin+1,chain_length,thining),]
  for(cp in 1:par_chains) parent_average_chain[[cp]]<-parent_average_chain[[cp]][seq(burnin+1,chain_length,thining),]
  
  return(list(chain_sets,parent_sigma_chain,parent_average_chain,DIC))
  
}








##############
# 
# # results=list(chain_sets,parent_sigma_chain,parent_average_chain,DIC)
# # save(results, file="results.RData")
# load("results.RData")
# 
# results[[1]][[1]]
# 
# boundary=c(); a=c(); b=c()
# for(cp in 1:par_chains) boundary=c(boundary, results[[1]][[1]][[cp]][, 1])
# for(cp in 1:par_chains) a=c(a, results[[1]][[1]][[cp]][, 2])
# for(cp in 1:par_chains) b=c(b, results[[1]][[1]][[cp]][, 3])
# hist(b)
# 
# #Plot Boundary
# iter=rep(1:100, 24); chain=rep(1:24, each=100)
# data.frame(boundary=boundary, iter=iter, chain=chain) %>% 
#   mutate(chain=factor(chain)) %>% 
#   ggplot() +
#   geom_density(aes(boundary, color=chain)) +
#   geom_density(aes(boundary), size=2) +
#   ggpubr::theme_pubr() + theme(legend.position = "none")
# 
# data.frame(boundary=boundary, iter=iter, chain=chain) %>% 
#   mutate(chain=factor(chain)) %>% 
#   ggplot(aes(iter, boundary, color=chain)) +
#   geom_line() +
#   facet_wrap(~chain) +
#   ggpubr::theme_pubr() + theme(legend.position = "none")
# 
# 
# #Plot a
# iter=rep(1:100, 24); chain=rep(1:24, each=100)
# data.frame(a=a, iter=iter, chain=chain) %>% 
#   mutate(chain=factor(chain)) %>% 
#   ggplot() +
#   geom_density(aes(a, color=chain)) +
#   geom_density(aes(a), size=2) +
#   ggpubr::theme_pubr() + theme(legend.position = "none")
# 
# data.frame(a=a, iter=iter, chain=chain) %>% 
#   mutate(chain=factor(chain)) %>% 
#   ggplot(aes(iter, a, color=chain)) +
#   geom_line() +
#   facet_wrap(~chain) +
#   ggpubr::theme_pubr() + theme(legend.position = "none")
# 
# 
# #Plot b
# iter=rep(1:100, 24); chain=rep(1:24, each=100)
# data.frame(b=b, iter=iter, chain=chain) %>% 
#   mutate(chain=factor(chain)) %>% 
#   ggplot() +
#   geom_density(aes(b, color=chain)) +
#   geom_density(aes(b), size=2) +
#   ggpubr::theme_pubr() + theme(legend.position = "none")
# 
# data.frame(b=b, iter=iter, chain=chain) %>% 
#   mutate(chain=factor(chain)) %>% 
#   ggplot(aes(iter, b, color=chain)) +
#   geom_line() +
#   facet_wrap(~chain) +
#   ggpubr::theme_pubr() + theme(legend.position = "none")
