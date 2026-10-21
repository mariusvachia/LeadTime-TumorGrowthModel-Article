results <- data.frame(matrix(nrow = 2500, ncol = 19))
# dossier_i <- list.dirs(path, recursive = FALSE)
distribution <- c("weibull","gompertz","gamma","lnorm","normal")
for(i in 1:5){
  dossier_i <- paste0("est/",distribution[i])# list.dirs(dossier_i[i], recursive = FALSE)
  for(j in 1:5){
    dossier_j <- paste0(dossier_i,"/",distribution[j])# list.files(dossier_j[j], full.names = TRUE)
    for(k in 1:100){
      fi <- paste0(dossier_j,"/mod_",k,".RData")
      if(file.exists(fi)){
        load(paste0(dossier_j,"/mod_",k,".RData"))
        results[(i-1)*500+(j-1)*100+k,] <- c(distribution[i],distribution[j],k,mod$convergence,mod$par[1:7],sqrt(diag(mod$vcov))[1:7],ifelse(mod$convergence == TRUE,mod$aic,NA))
      }
    }
  }
}
colnames(results) <- c("simulated", "estimated","rep", "convergence", "p1", "p2","p3","p4", "p5", "p6", "p7" ,
                       "sd_p1", "sd_p2", "sd_p3","sd_p4", "sd_p5","sd_p6", "sd_p7", "aic")

criterion <- data.frame(matrix(0, nrow = 5, ncol = 5))
for(i in 1:5){
  base_i <- subset(results, simulated == distribution[i])
  for(j in 1:100){
    base_j <- subset(base_i, rep == j)
    dist <- base_j$estimated[which.min(base_j$aic)]
    criterion[i,which(distribution == dist)] <- criterion[i,which(distribution == dist)] + 1
  }
}
rownames(criterion) <- distribution
colnames(criterion) <- distribution
criterion


nun <- numeric(0)
aic_res <- data.frame(matrix(0, nrow = 5, ncol = 15))
for(i in 1:5){
  base_i <- subset(results, simulated == distribution[i])
  nun <- c(nun,c(distribution[i],"sd","%"))
  for(j in 1:5){
    base_j <- subset(base_i, estimated == distribution[j])
    aic_res[i,3*j-2] <- mean(as.numeric(base_j$aic), na.rm = T)
    aic_res[i,3*j-1] <- sd(as.numeric(base_j$aic), na.rm = T)
    aic_res[i,3*j] <- sum(!is.na(as.numeric(base_j$aic))/100)
  }
}
rownames(aic_res) <- distribution
colnames(aic_res) <- nun
aic_res
