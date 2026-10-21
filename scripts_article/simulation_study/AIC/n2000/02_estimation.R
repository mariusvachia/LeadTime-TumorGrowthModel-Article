library(LeadTime)

load("env.RData")
distribution <- c("weibull","gompertz","gamma","lnorm","normal")
model.distrib <- c("weibull","gompertz","gamma","log-normal","normal")
init <- list(weibull = c(2,4),
             gompertz = c(-4,-2),
             gamma = c(-2,4),
             lnorm = c(4,-2),
             normal = c(4,2))

for(rep in 1:2500){
  i <- ((rep-1)%/%500)+1
  j <- (((rep-1)%%500)%/%100)+1
  k <- (rep-1)%%100+1
  load(paste0("datasets/",distribution[i],"/sim_",k,".RData"))

  mod <- TumorGrowthModel(data_simu,
                          initpar = c(-0.8,-0.2,-8.75,-5,0.5,init[[j]]),
                          distribution.t0 = model.distrib[j],
                          n_cores = 16)

  # save(mod, file = paste0("est/",distribution[i],"/",distribution[j],"/mod_",k,".RData"))
}
