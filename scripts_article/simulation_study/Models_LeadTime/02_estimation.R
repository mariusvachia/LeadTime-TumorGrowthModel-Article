library(LeadTime)
#-------------------------------------------------------------------------------
# Scenario 1
#-------------------------------------------------------------------------------
for(i in 1:200){
  load(paste0("datasets/scenario1/sim_",i,".RData"))
  mod <- LeadTime::TumorGrowthModel(data_simu,
                                    d0 = 0.5,
                                    distribution.t0 = "normal",
                                    initpar = c(-1,-1,-9,-5,0.5,4,2),
                                    n_cores = 18)
  # save(mod, file = paste0("est/scenario1/mod_",i,".RData"))
}

results <- data.frame(matrix(nrow = 200, ncol = 16))
for(i in 1:200){
  load(paste0("est/scenario1/mod_",i,".RData"))
  results[i,] <- c(i,mod$convergence,mod$par,sqrt(diag(mod$vcov)))
}

colnames(results) <- c("i", "istop", "est_tau1", "est_tau2", "est_eta", "est_beta0", "est_beta1", "est_mu", "est_sigma", "sd_tau1", "sd_tau2", "sd_eta", "sd_beta0", "sd_beta1", "sd_mu", "sd_sigma")
df <- as.data.frame(results)
table(df$istop)
TruePar <- round(c(-log(2.36),log(2.36) - log(3.00),-8.75,-5.04,0.56,log(55),log(11.8)),3)
tab <- data.frame(matrix(NA, nrow = 7, ncol = 5))
for(j in 1:7){
  tab[j,] <- c(TruePar[j],
               round(mean(df[df$istop == 1,2+j]),3),
               round(mean(df[df$istop == 1,9+j]),3),
               round(sd(df[df$istop == 1,2+j]),3),
               round(sum(TruePar[j] >= df[df$istop == 1,2+j]-1.96*df[df$istop == 1,9+j] & TruePar[j] <= df[df$istop == 1,2+j]+1.96*df[df$istop == 1,9+j])/nrow(df[df$istop == 1,]),3))
}
tab[,6] <- round((tab[,2] - tab[,1]),3)
tab[,7] <- round((tab[,2] - tab[,1])/tab[,4],3)
colnames(tab) <- c("True Value","Est","model SD","empirical SD","coverage","abs bias","std bias")
rownames(tab) <- c("-log(tau1)","log(tau1) - log(tau2)","log(eta)","beta0","beta1","log(mu)","log(sigma)")
tab

#-------------------------------------------------------------------------------
# Scenario 2
#-------------------------------------------------------------------------------
for(i in 1:200){
  load(paste0("datasets/scenario2/sim_",i,".RData"))
  mod <- LeadTime::TumorGrowthModel(data_simu,
                                    d0 = 0.5,
                                    distribution.t0 = "normal",
                                    initpar = c(-1,-1,-9,-5,0.5,4,2),
                                    n_cores = 18)
  # save(mod, file = paste0("est/scenario2/mod_",i,".RData"))
}

results <- data.frame(matrix(nrow = 200, ncol = 16))
for(i in 1:200){
  load(paste0("est/scenario2/mod_",i,".RData"))
  results[i,] <- c(i,mod$convergence,mod$par,sqrt(diag(mod$vcov)))
}

colnames(results) <- c("i", "istop", "est_tau1", "est_tau2", "est_eta", "est_beta0", "est_beta1", "est_mu", "est_sigma", "sd_tau1", "sd_tau2", "sd_eta", "sd_beta0", "sd_beta1", "sd_mu", "sd_sigma")
df <- as.data.frame(results)
table(df$istop)
TruePar <- round(c(-log(2.36),log(2.36) - log(3.00),-8.75,-5.04,0.56,log(55),log(11.8)),3)
tab <- data.frame(matrix(NA, nrow = 7, ncol = 5))
for(j in 1:7){
  tab[j,] <- c(TruePar[j],
               round(mean(df[df$istop == 1,2+j]),3),
               round(mean(df[df$istop == 1,9+j]),3),
               round(sd(df[df$istop == 1,2+j]),3),
               round(sum(TruePar[j] >= df[df$istop == 1,2+j]-1.96*df[df$istop == 1,9+j] & TruePar[j] <= df[df$istop == 1,2+j]+1.96*df[df$istop == 1,9+j])/nrow(df[df$istop == 1,]),3))
}
tab[,6] <- round((tab[,2] - tab[,1]),3)
tab[,7] <- round((tab[,2] - tab[,1])/tab[,4],3)
colnames(tab) <- c("True Value","Est","model SD","empirical SD","coverage","abs bias","std bias")
rownames(tab) <- c("-log(tau1)","log(tau1) - log(tau2)","log(eta)","beta0","beta1","log(mu)","log(sigma)")
tab


#-------------------------------------------------------------------------------
# Scenario 3
#-------------------------------------------------------------------------------
for(i in 1:200){
  load(paste0("datasets/scenario3/sim_",i,".RData"))
  mod <- LeadTime::TumorGrowthModel(data_simu,
                                    d0 = 0.5,
                                    distribution.t0 = "normal",
                                    initpar = c(-1,-1,-9,-5,0.5,4,2),
                                    n_cores = 18)
  # save(mod, file = paste0("est/scenario3/mod_",i,".RData"))
}

results <- data.frame(matrix(nrow = 200, ncol = 16))
for(i in 1:200){
  load(paste0("est/scenario3/mod_",i,".RData"))
  results[i,] <- c(i,mod$convergence,mod$par,sqrt(diag(mod$vcov)))
}

colnames(results) <- c("i", "istop", "est_tau1", "est_tau2", "est_eta", "est_beta0", "est_beta1", "est_mu", "est_sigma", "sd_tau1", "sd_tau2", "sd_eta", "sd_beta0", "sd_beta1", "sd_mu", "sd_sigma")
df <- as.data.frame(results)
table(df$istop)
TruePar <- round(c(-log(2.36),log(2.36) - log(3.00),-8.75,-5.04,0.56,log(55),log(11.8)),3)
tab <- data.frame(matrix(NA, nrow = 7, ncol = 5))
for(j in 1:7){
  tab[j,] <- c(TruePar[j],
               round(mean(df[df$istop == 1,2+j]),3),
               round(mean(df[df$istop == 1,9+j]),3),
               round(sd(df[df$istop == 1,2+j]),3),
               round(sum(TruePar[j] >= df[df$istop == 1,2+j]-1.96*df[df$istop == 1,9+j] & TruePar[j] <= df[df$istop == 1,2+j]+1.96*df[df$istop == 1,9+j])/nrow(df[df$istop == 1,]),3))
}
tab[,6] <- round((tab[,2] - tab[,1]),3)
tab[,7] <- round((tab[,2] - tab[,1])/tab[,4],3)
colnames(tab) <- c("True Value","Est","model SD","empirical SD","coverage","abs bias","std bias")
rownames(tab) <- c("-log(tau1)","log(tau1) - log(tau2)","log(eta)","beta0","beta1","log(mu)","log(sigma)")
tab
