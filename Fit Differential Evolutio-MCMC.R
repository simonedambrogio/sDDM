# Social Informations: {-7, -5, -3, -1, 1, 3, 5, 7}
# Temptation (Closest Value): -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5
# We have a table of 8 x 11 = 88 cells.

RcppParallel::setThreadOptions(numThreads = 1) #this is critical for running on Mac OS for some reason.
social_info=seq(-7, 7, 2)
temptation=-5:5

#Load Libraries
library(dplyr)
library(ggplot2); library(ggpubr); library(RColorBrewer)
theme=theme(text = element_text(size=15), strip.background = element_rect(fill="white"))
library(DEoptim)
df_lapply=function(X, FUN) do.call(rbind.data.frame, lapply(X, FUN)) %>% tibble::as_tibble()

############################ Trials per cell ############################
n_trial=3; print(paste('Each participant has to complite', n_trial*length(social_info)*length(temptation) ,'trials'))
n_subject=50

#Create variable with n_trial per cell
social=rep(social_info, each=n_trial)

#Replicate social for each temptation level
social_temptation=df_lapply(1:length(temptation), function(t_i) data.frame(social=social, temptation=temptation[t_i] ))

#Replicate social_temptation for each subject
data_sim=df_lapply(1:n_subject, function(s_i) rbind(cbind(data.frame(Subject=s_i), social_temptation)))

##### Simulate choice and rts #####
Rcpp::sourceCpp('Functions/rsddm.cpp')

nDT=0.3; dt=0.005

dp=.4; B_intercept=0.15; B_slope=1; theta=1; qu=0.9; s=0.1

data.frame(temptation=seq(-5, 5, 0.1), 
           B = 1/(1+exp(B_intercept*(seq(-5, 5, 0.1)-B_slope)))
) %>% ggplot(aes(temptation, B)) +
  geom_line(color='firebrick4', size=2)+
  theme_pubr() + theme +
  geom_hline(yintercept = 0.5, linetype=2, color="gray") +
  geom_vline(xintercept = 0, linetype=2, color="gray") +
  scale_y_continuous(limits = c(0,1))

dp_s=rnorm(n_subject,dp,0.1); B_intercept_s=rnorm(n_subject,B_intercept,0.03); B_slope_s=rnorm(n_subject,B_slope,0.1)
theta_s=rnorm(n_subject,theta,0.1); qu_s=runif(n_subject,qu-qu/3,qu+qu/3); s_s=runif(n_subject,s-s/3,s+s/3)

rts=lapply(1:n_subject, function(s_i){
  data=data_sim[data_sim$Subject==s_i, ]
  dp=dp_s[s_i]; B_intercept=B_intercept_s[s_i]; B_slope=B_slope_s[s_i]; theta=theta_s[s_i]; qu=qu_s[s_i]; s=s_s[s_i]
  sapply(1:nrow(data), function(row_i){
    sddm_parallel(dp=dp, theta=theta, nDT=nDT, B_intercept = B_intercept, B_slope = B_slope, 
                  temptation = data$temptation[row_i], social_strength=data$social[row_i], 
                  dt=dt, qu=qu, s=s, N=1)
  })
}) %>% unlist

data_sim=data_sim %>% mutate(rts=rts, choice=ifelse(rts<0, 0, 1), rts=abs(rts))

# Plots 
data_sim %>% 
  mutate(social=factor(social, labels = seq(-7, 7, 2)), Subject=factor(Subject),
         choice=ifelse(choice==0, 1, 0)) %>%  #Change choice if we wnat to plot the probability of cheating
  ggplot() +
  #stat_smooth(aes(temptation, choice, group=Subject), geom='line', alpha=0.5, se=FALSE, method = 'glm', color='gray', alpha=0.2, method.args=list('binomial'))+
  geom_smooth(aes(temptation, choice, color=social), method = 'glm', method.args=list('binomial'), se=F, size=2) +
  theme_pubr() + theme +
  scale_color_brewer(palette = 'Set1') +
  scale_x_continuous(breaks = -5:5) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x="Temptation", y="P(Cheating)") +
  facet_wrap(~social)


ggplot() +
  # geom_line(data = df2_1, aes(RV, fit), size = 1, color="navy") +
  #geom_line(aes(RV, rts, color=Natt)) +
  geom_point(data=data_sim %>% 
               Rmisc::summarySEwithin(., measurevar = "rts", withinvars = c("social", "temptation", "Subject")) %>% 
               mutate(temptation=as.numeric(as.character(temptation))),
             aes(temptation, rts, group=Subject), size=1, color="gray", alpha=0.5
             ) +
  geom_pointrange(data=data_sim %>% 
                    Rmisc::summarySEwithin(., measurevar = "rts", withinvars = c("social", "temptation")) %>% 
                    mutate(temptation=as.numeric(as.character(temptation))),
                  aes(temptation, rts, ymin=rts-se, ymax=rts+se, color=social), size=1) +
  labs(x="Temptation", y="Mean Reaction Time") +
  scale_x_continuous(breaks = -5:5) +
  theme_pubr() + theme +
  scale_color_brewer(palette = 'Set1') +
  facet_wrap(~social)

data_sim %>% 
  mutate(rts=ifelse(choice==1, rts, -rts), social=factor(social)) %>% 
  ggplot(aes(rts, color=social, fill=social)) +
  geom_density(size=2, alpha=1/6) +
  scale_x_continuous(limits = c(-5,5)) +
  theme_pubr() + theme +
  scale_color_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Set1') +
  facet_wrap(~social)

############# Recovery Parameters #############
#Hierarchical Differential Evolution MCMC
source('Functions/DE-MCMC.R')
source('Functions/all_lik.R')

# Parameters used: dp=.4; theta=1; qu=0.5; s=0.5; B_intercept=0.5; B_slope=0.5;
bound_dp=c(dp-dp/1.1, dp+dp/0.2)
bound_theta=c(theta-theta/1.1, theta+theta/0.4)
bound_s=c(s-s/1.1, s+s/0.2)
bound_qu=c(qu-qu/2, qu+qu/0.5)
bound_B_intercept=c(B_intercept-B_intercept/1.1, B_intercept+B_intercept/0.2)
bound_B_slope=c(B_slope-B_slope/1.1, B_slope+B_slope/0.5)
#Minimal simga value to reduce likelyhood of chains getting stucked (not sure if necessary)
min_sd=0.01; h=0.005; N=3000

chain_length=6000; burnin=2000; thining=10; par_chains=24; cores=4

fit_s= hdemcmc(data_sim, nDT=0.3, N=3000, dt=0.005,
               chain_length=chain_length,
               burnin=burnin,
               thining=thining,
               par_chains=par_chains, #parallel chains
               cores=cores); save(fit_s, file = "Result/fit_s.RData")

load(file = "Result/fit_s.RData")

#Create Parent Posterior data_plot
avg_dp=c(); avg_theta=c(); avg_qu=c(); avg_s=c(); avg_B_intercept=c(); avg_B_slope=c();
for(cp in 1:par_chains) avg_dp=c(avg_dp, fit_s[[3]][[cp]][, 1])
for(cp in 1:par_chains) avg_theta=c(avg_theta, fit_s[[3]][[cp]][, 2])
for(cp in 1:par_chains) avg_s=c(avg_s, fit_s[[3]][[cp]][, 3])
for(cp in 1:par_chains) avg_qu=c(avg_qu, fit_s[[3]][[cp]][, 4])
for(cp in 1:par_chains) avg_B_intercept=c(avg_B_intercept, fit_s[[3]][[cp]][, 5])
for(cp in 1:par_chains) avg_B_slope=c(avg_B_slope, fit_s[[3]][[cp]][, 6])

iteration = rep(1:((chain_length-burnin)/thining), par_chains)
iteration_s = rep(iteration, n_subject)

chain = rep(1:par_chains, each=((chain_length-burnin)/thining))
chain_s = rep(chain, n_subject)

data_plot_avg = data.frame(dp=avg_dp,
                           theta=avg_theta,
                           s=avg_s,
                           qu=avg_qu,
                           B_intercept=avg_B_intercept,
                           B_slope=avg_B_slope,
                           iteration=iteration, chain=factor(chain))

apply(data_plot_avg[,1:6], 2, mean)

#Plot dp
data_plot_avg %>% 
  ggplot() +
  #geom_density(aes(boundary, color=chain)) +
  geom_density(aes(dp), size=2, alpha=1/4,
               color=brewer.pal(n = 8, name = "Blues")[8], 
               fill=brewer.pal(n = 8, name = "Blues")[8]) +
  ggpubr::theme_pubr() + 
  theme(legend.position = "none")
#density(data_plot_avg$dp, n = 9600)$x[which.max(density(data_plot_avg$dp, n = 9600)$y)]

data_plot_avg %>%
  ggplot(aes(iteration, dp, color=chain)) +
  geom_line() +
  ggpubr::theme_pubr() + 
  scale_color_manual(values = sample(brewer.pal(n = 8, name = "Blues"), 24, replace = T)) +
  theme(legend.position = "none")

#Plot Intercept
data_plot_avg %>% 
  mutate(chain=factor(chain)) %>% 
  ggplot() +
  geom_density(aes(theta), size=2, alpha=1/4,
               color=brewer.pal(n = 8, name = "RdPu")[8],
               fill=brewer.pal(n = 8, name = "RdPu")[8]) +
  ggpubr::theme_pubr() + 
  theme(legend.position = "none")
#density(data_plot_avg$theta, n = 9600)$x[which.max(density(data_plot_avg$theta, n = 9600)$y)]

data_plot_avg %>%
  ggplot(aes(iteration, theta, color=chain)) +
  geom_line() +
  ggpubr::theme_pubr() + 
  scale_color_manual(values = sample(brewer.pal(n = 8, name = "RdPu"), 24, replace = T)) +
  theme(legend.position = "none")

#Plot s
data_plot_avg %>% 
  mutate(chain=factor(chain)) %>% 
  ggplot() +
  geom_density(aes(s), size=2, alpha=1/4,
               color=brewer.pal(n = 8, name = "BuGn")[8],
               fill=brewer.pal(n = 8, name = "BuGn")[8]) +
  ggpubr::theme_pubr() + 
  theme(legend.position = "none")
#density(data_plot_avg$s, n = 9600)$x[which.max(density(data_plot_avg$s, n = 9600)$y)]

data_plot_avg %>%
  ggplot(aes(iteration, s, color=chain)) +
  geom_line() +
  ggpubr::theme_pubr() + 
  scale_color_manual(values = sample(brewer.pal(n = 8, name = "BuGn"), 24, replace = T)) +
  theme(legend.position = "none")

#Plot 
data_plot_avg %>% 
  mutate(chain=factor(chain)) %>% 
  ggplot() +
  geom_density(aes(B_intercept), size=2, alpha=1/4,
               color=brewer.pal(n = 8, name = "BuGn")[8],
               fill=brewer.pal(n = 8, name = "BuGn")[8]) +
  ggpubr::theme_pubr() + 
  theme(legend.position = "none")
#density(data_plot_avg$B_intercept, n = 9600)$x[which.max(density(data_plot_avg$B_intercept, n = 9600)$y)]

data_plot_avg %>%
  ggplot(aes(iteration, B_intercept, color=chain)) +
  geom_line() +
  ggpubr::theme_pubr() + 
  scale_color_manual(values = sample(brewer.pal(n = 8, name = "BuGn"), 24, replace = T)) +
  theme(legend.position = "none")

#Plot Slope
data_plot_avg %>% 
  mutate(chain=factor(chain)) %>% 
  ggplot() +
  geom_density(aes(B_slope), size=2, alpha=1/4,
               color=brewer.pal(n = 8, name = "BuGn")[8],
               fill=brewer.pal(n = 8, name = "BuGn")[8]) +
  ggpubr::theme_pubr() + 
  theme(legend.position = "none")
#density(data_plot_avg$B_slope, n = 9600)$x[which.max(density(data_plot_avg$B_slope, n = 9600)$y)]

data_plot_avg %>%
  ggplot(aes(iteration, B_slope, color=chain)) +
  geom_line() +
  ggpubr::theme_pubr() + 
  scale_color_manual(values = sample(brewer.pal(n = 8, name = "BuGn"), 24, replace = T)) +
  theme(legend.position = "none") +
  facet_wrap(~chain)

#Subject Posteriors
#Create Subject Posterior data_plot
dp=c(); theta=c(); qu=c(); s=c(); B_intercept=c(); B_slope=c(); 
for(j in Subjects) for(cp in 1:par_chains) dp=c(dp, fit_s[[1]][[j]][[cp]][, 1])
for(j in Subjects) for(cp in 1:par_chains) theta=c(theta, fit_s[[1]][[j]][[cp]][, 2])
for(j in Subjects) for(cp in 1:par_chains) s=c(s, fit_s[[1]][[j]][[cp]][, 3])
for(j in Subjects) for(cp in 1:par_chains) qu=c(qu, fit_s[[1]][[j]][[cp]][, 4])
for(j in Subjects) for(cp in 1:par_chains) B_intercept=c(B_intercept, fit_s[[1]][[j]][[cp]][, 5])
for(j in Subjects) for(cp in 1:par_chains) B_slope=c(B_slope, fit_s[[1]][[j]][[cp]][, 6])

data_plot = data.frame(Subject = factor(rep(1:n_subject, each=length(iteration))),
                       dp=dp,
                       theta=theta,
                       qu=qu,
                       s=s,
                       B_intercept=B_intercept,
                       B_slope=B_slope,
                       iteration=iteration_s, chain=factor(chain_s))

#Personal drift Posterior
data_plot %>% 
  ggplot() +
  geom_density(aes(dp, color=Subject)) +
  #geom_density(aes(dp), size=2) +
  ggpubr::theme_pubr() + 
  theme(legend.position = "none")

data_plot %>%
  ggplot(aes(iteration, dp, color=chain)) +
  geom_line() +
  facet_wrap(~Subject) +
  ggpubr::theme_pubr() + 
  scale_color_manual(values = sample(brewer.pal(n = 8, name = "Blues"), 24, replace = T)) +
  theme(legend.position = "none")

#theta Posterior
data_plot %>% 
  ggplot() +
  geom_density(aes(theta, color=Subject)) +
  ggpubr::theme_pubr() + 
  theme(legend.position = "none")

data_plot %>%
  ggplot(aes(iteration, theta, color=chain)) +
  geom_line() +
  facet_wrap(~Subject) +
  ggpubr::theme_pubr() + 
  scale_color_manual(values = sample(brewer.pal(n = 8, name = "Blues"), 24, replace = T)) +
  theme(legend.position = "none")

#s Posterior
data_plot %>% 
  ggplot() +
  geom_density(aes(s, color=Subject)) +
  ggpubr::theme_pubr() + 
  theme(legend.position = "none")

data_plot %>%
  ggplot(aes(iteration, s, color=chain)) +
  geom_line() +
  facet_wrap(~Subject) +
  ggpubr::theme_pubr() + 
  scale_color_manual(values = sample(brewer.pal(n = 8, name = "Blues"), 24, replace = T)) +
  theme(legend.position = "none")

#qu Posterior
data_plot %>% 
  ggplot() +
  geom_density(aes(qu, color=Subject)) +
  ggpubr::theme_pubr() + 
  theme(legend.position = "none")

data_plot %>%
  ggplot(aes(iteration, qu, color=chain)) +
  geom_line() +
  facet_wrap(~Subject) +
  ggpubr::theme_pubr() + 
  scale_color_manual(values = sample(brewer.pal(n = 8, name = "Blues"), 24, replace = T)) +
  theme(legend.position = "none")

#B_intercept Posterior
data_plot %>% 
  ggplot() +
  geom_density(aes(B_intercept, color=Subject)) +
  ggpubr::theme_pubr() + 
  theme(legend.position = "none")

data_plot %>%
  ggplot(aes(iteration, B_intercept, color=chain)) +
  geom_line() +
  facet_wrap(~Subject) +
  ggpubr::theme_pubr() + 
  scale_color_manual(values = sample(brewer.pal(n = 8, name = "Blues"), 24, replace = T)) +
  theme(legend.position = "none")

#B_slope Posterior
data_plot %>% 
  ggplot() +
  geom_density(aes(B_slope, color=Subject)) +
  ggpubr::theme_pubr() + 
  theme(legend.position = "none")

data_plot %>%
  ggplot(aes(iteration, B_slope, color=chain)) +
  geom_line() +
  facet_wrap(~Subject) +
  ggpubr::theme_pubr() + 
  scale_color_manual(values = sample(brewer.pal(n = 8, name = "Blues"), 24, replace = T)) +
  theme(legend.position = "none")
