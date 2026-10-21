library(LeadTime)
#-------------------------------------------------------------------------------
# Simulation, T0 normal
#-------------------------------------------------------------------------------
n <- 2000
p = list(p1 = 2.36, p2 = 3.00, p3 = exp(-8.75), p4 = -5.04, p5 = 0.56, p6 = 55, p7 = 11.8)
tau1 <- p$p1
tau2 <- p$p2
eta <- p$p3
beta0 <- p$p4
beta1 <- p$p5
mu <- p$p6
sigma <- p$p7
par <- c(-log(tau1),log(tau1)-log(tau2),log(eta),beta0,beta1,log(mu),log(sigma))
#-------------------------------------------------------------------------------
# Scenario 1, full screening, 100% participation & 100% visits
#-------------------------------------------------------------------------------
load(paste0("seeds/seed_200.RData"))
for(i in 1:200){
  set.seed(seed_200[i])
  data_simu <- SimuTGM(n,
                       p,
                       d0 = 0.5,
                       distribution.t0 = "normal",
                       screening_start = 50,
                       screening_stop = 75,
                       screening_interval = 2,
                       prob_program = 1,
                       prob_visit = 1,
                       nor1 = 0.2,
                       nor2 = 0.15,
                       surv_size = 6e-4,
                       surv_growth = 4e-3)
  #save(data_simu, file = paste0(path,"simu200/sim_",i,".RData"))
}


#-------------------------------------------------------------------------------
# Scenario 2, realistic, 75% participation & 75% visits
#-------------------------------------------------------------------------------
load(paste0("seeds/seed_rea75.RData"))
for(i in 1:200){
  set.seed(seed_rea75[i])
  data_simu <- SimuTGM(n,
                       p,
                       d0 = 0.5,
                       distribution.t0 = "normal",
                       screening_start = 50,
                       screening_stop = 75,
                       screening_interval = 2,
                       prob_program = 0.75,
                       prob_visit = 0.75,
                       nor1 = 0.2,
                       nor2 = 0.15,
                       surv_size = 6e-4,
                       surv_growth = 4e-3)
  # save(data_simu, file = paste0(path,"simu_rea75/sim_",i,".RData"))
}

#-------------------------------------------------------------------------------
# Scenario 3, extremely low screening, 25% participation & 25% visits
#-------------------------------------------------------------------------------
load(paste0("seeds/seed_ext25.RData"))
for(i in 1:200){
  set.seed(seed_ext25[i])
  data_simu <- SimuTGM(n,
                       p,
                       d0 = 0.5,
                       distribution.t0 = "normal",
                       screening_start = 50,
                       screening_stop = 75,
                       screening_interval = 2,
                       prob_program = 0.25,
                       prob_visit = 0.25,
                       nor1 = 0.2,
                       nor2 = 0.15,
                       surv_size = 6e-4,
                       surv_growth = 4e-3)
  # save(data_simu, file = paste0(path,"simu_ext25/sim_",i,".RData"))
}


#-------------------------------------------------------------------------------
# Description
#-------------------------------------------------------------------------------
# scenario 1 :
scenario1 <- as.data.frame(matrix(numeric(0), nrow = 200, ncol = 14))
colnames(scenario1) <- c("rep","nsym","nscr","nprog","total_screening","tumor_size_med","tumor_size_Q1","tumor_size_Q3",
                         "dsym_med","dsym_Q1","dsym_Q3","dscr_med","dscr_Q1","dscr_Q3")
for(i in 1:200){
  load(paste0("datasets/scenario1/sim_",i,".RData"))
  total_screening <- 0
  nprog <- 0
  for(j in 1:nrow(data_simu)){
    total_screening <- total_screening + length(data_simu$tscr[[j]])
    nprog <- nprog + ifelse(length(data_simu$tscr[[j]]) > 0, 1,0)
  }
  scenario1[i,] <- c(i,sum(data_simu$det == 1),sum(data_simu$det == 2), nprog, total_screening,
                     quantile(data_simu$ddet,c(0.5,0.25,0.75)),
                     quantile(data_simu$ddet[data_simu$det==1],c(0.5,0.25,0.75)),
                     quantile(data_simu$ddet[data_simu$det==2],c(0.5,0.25,0.75)))
}
colMeans(scenario1[, c(2:5)], na.rm = TRUE)
apply(scenario1[,c(6:14)], 2, median, na.rm = TRUE)


# scenario 2 :
scenario2 <- as.data.frame(matrix(numeric(0), nrow = 200, ncol = 14))
colnames(scenario2) <- c("rep","nsym","nscr","nprog","total_screening","tumor_size_med","tumor_size_Q1","tumor_size_Q3",
                         "dsym_med","dsym_Q1","dsym_Q3","dscr_med","dscr_Q1","dscr_Q3")
for(i in 1:200){
  load(paste0("datasets/scenario2/sim_",i,".RData"))
  total_screening <- 0
  nprog <- 0
  for(j in 1:nrow(data_simu)){
    total_screening <- total_screening + length(data_simu$tscr[[j]])
    nprog <- nprog + ifelse(length(data_simu$tscr[[j]]) > 0, 1,0)
  }
  scenario2[i,] <- c(i,sum(data_simu$det == 1),sum(data_simu$det == 2),nprog , total_screening,
                     quantile(data_simu$ddet,c(0.5,0.25,0.75)),
                     quantile(data_simu$ddet[data_simu$det==1],c(0.5,0.25,0.75)),
                     quantile(data_simu$ddet[data_simu$det==2],c(0.5,0.25,0.75)))
}
colMeans(scenario2[, c(2:5)], na.rm = TRUE)
apply(scenario2[,c(6:14)], 2, median, na.rm = TRUE)




# scenario 3 :
scenario3 <- as.data.frame(matrix(numeric(0), nrow = 200, ncol = 14))
colnames(scenario3) <- c("rep","nsym","nscr","nprog","total_screening","tumor_size_med","tumor_size_Q1","tumor_size_Q3",
                         "dsym_med","dsym_Q1","dsym_Q3","dscr_med","dscr_Q1","dscr_Q3")
for(i in 1:200){
  load(paste0("datasets/scenario3/sim_",i,".RData"))
  total_screening <- 0
  nprog <- 0
  for(j in 1:nrow(data_simu)){
    total_screening <- total_screening + length(data_simu$tscr[[j]])
    nprog <- nprog + ifelse(length(data_simu$tscr[[j]]) > 0, 1,0)
  }
  scenario3[i,] <- c(i,sum(data_simu$det == 1),sum(data_simu$det == 2),nprog, total_screening,
                     quantile(data_simu$ddet,c(0.5,0.25,0.75)),
                     quantile(data_simu$ddet[data_simu$det==1],c(0.5,0.25,0.75)),
                     quantile(data_simu$ddet[data_simu$det==2],c(0.5,0.25,0.75)))
}
colMeans(scenario3[, c(2:5)], na.rm = TRUE)
apply(scenario3[,c(6:14)], 2, median, na.rm = TRUE)
