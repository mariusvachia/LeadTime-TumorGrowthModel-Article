library(LeadTime)

for(i in 1:200){
  load(paste0("datasets/scenario2/sim_",i,".RData"))
  load(paste0("est/scenario2/mod_",i,".RData"))
  cor <- LeadTime::LeadTimeCor(data_simu,
                               variables=list(id = "id",
                                              det = "det",
                                              tdet = "tdet",
                                              ddet = "ddet",
                                              scr = "scr",
                                              tscr = "tscr",
                                              last = "death",
                                              tlast = "tdeath_w1"),
                               model = mod,
                               overdiagnosis = FALSE,
                               d0 = 0.1,
                               distribution.t0 = "normal",
                               n_cores = 18)
  # save(cor, file = paste0("cor/cor_",i,".RData"))
}

load(paste0("datasets/scenario2/sim_",1,".RData"))
data_simu$rep <- 1
data_true_value <- subset(data_simu, det == 2)
for(i in 2:200){
  load(paste0("datasets/scenario2/sim_",i,".RData"))
  data_simu$rep <- i
  data_true_value <- rbind(data_true_value, subset(data_simu, det == 2))
}
nrow(data_true_value)
nrow(data_true_value)/200
summary(round(data_true_value$surv_w1,1))
summary(round(data_true_value$lead_time,2))

survival_true_value <- survival::survfit(survival::Surv(data_true_value$tdeath_w1 - data_true_value$tsym, data_true_value$death) ~ 1)
summary(survival_true_value, times = c(2, 5, 10))
3243/200
(4572+3243)/200
(7230+4572+3243)/200

lt_est <- data.frame(matrix(nrow = 0, ncol = 7))
surv_res <- data.frame(matrix(nrow = 200, ncol = 6))
for(i in 1:200){
  load(paste0("datasets/scenario2/sim_",i,".RData"))
  load(paste0("cor/cor_",i,".RData"))
  df <- data.frame(rep = rep(i, sum(cor$detection == 2)),
                   id = cor$id[cor$detection == 2],
                   observed = cor$observed[cor$detection == 2],
                   leadtime = cor$leadtime[cor$detection == 2],
                   corrected = cor$corrected[cor$detection == 2],
                   lower_pred = cor$predint$`2.5%`,
                   upper_pred = cor$predint$`97.5%`)
  lt_est <- rbind(lt_est, df)
  surv_res[i,] <- c(summary(cor$survival$cor, times = c(2, 5, 10))$surv, summary(cor$survival$cor, times = c(2, 5, 10))$std.err)
}

TruePar <- round(summary(survival_true_value, times = c(2, 5, 10))$surv, 3)
tab <- data.frame(matrix(NA, nrow = 3, ncol = 5))
for(j in 1:3){
  tab[j,] <- c(TruePar[j],
               round(mean(surv_res[,j]),3),
               round(mean(surv_res[,3+j]),3),
               round(sd(surv_res[,j]),3),
               round(sum(TruePar[j] >= surv_res[,j]-1.96*surv_res[,3+j] & TruePar[j] <= surv_res[,j]+1.96*surv_res[,3+j])/nrow(surv_res),3))
}
tab[,6] <- round((tab[,2] - tab[,1]),3)
tab[,7] <- round((tab[,2] - tab[,1])/tab[,4],3)
colnames(tab) <- c("True Value","Est","est SE","empirical SE","coverage","abs bias","std bias")
rownames(tab) <- c("2 years", "5 years", "10 years")
tab




summary(data_true_value$lead_time)
summary(lt_est$leadtime)
summary(lt_est$leadtime - data_true_value$lead_time)



agg_lt_true <- data_true_value$lead_time
agg_lt_est <- numeric(0)
for(i in 1:200){
  load(paste0("cor/cor_",i,".RData"))
  agg_lt_est <- c(agg_lt_est,cor$leadtime[cor$detection == 2])
}


plot(agg_lt_true,agg_lt_est, ylim = c(0,35))
abline(a = 0, b = 1, col = "red")
hist(agg_lt_true-agg_lt_est,200,xlim = c(-10,10))
