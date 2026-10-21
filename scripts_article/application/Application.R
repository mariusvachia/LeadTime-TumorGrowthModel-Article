#-------------------------------------------------------------------------------
# # please load functions in file ("LeadTime/R/functions.R")
#-------------------------------------------------------------------------------

load("ExampleData.RData")
names(ExampleData)[2] <- "det"
names(ExampleData)[3] <- "ddet"
names(ExampleData)[4] <- "tdet"
names(ExampleData)[5] <- "scr" # Screening test results
names(ExampleData)[6] <- "tscr" # Age at screening test results
names(ExampleData)[7] <- "last" # Vital status at last follow-up
names(ExampleData)[8] <- "tlast"  # Age at last follow-up
mod_normal0.1 <- LeadTime::TumorGrowthModel(ExampleData,
                                            d0 = 0.1,
                                            distribution.t0 = "normal",
                                            initpar = c(-0.8415672, -0.2570451, -8.7500000, -4.7500000, 0.5600000, 4.0073332, 2.4680995),
                                            n_cores = 18)
# save(mod_normal0.1, file = "mod_normal0.1.RData")

mod_gamma0.1 <- LeadTime::TumorGrowthModel(ExampleData,
                                           d0 = 0.1,
                                           distribution.t0 = "gamma",
                                           initpar = c(-0.8415672, -0.2570451, -8.7500000, -4.7500000, 0.5600000, -2.7725887,  3.9765615),
                                           n_cores = 18)
# save(mod_gamma0.1, file = "mod_gamma0.1.RData")

mod_weibull0.1 <- LeadTime::TumorGrowthModel(ExampleData,
                                             d0 = 0.1,
                                             distribution.t0 = "weibull",
                                             initpar = c(-0.8415672, -0.2570451, -8.7500000, -4.7500000, 0.5600000, 1.6094379, 4.0943446),
                                             n_cores = 18)
# save(mod_weibull0.1, file = "mod_weibull0.1.RData")

mod_gompertz0.1 <- LeadTime::TumorGrowthModel(ExampleData,
                                              d0 = 0.1,
                                              distribution.t0 = "gompertz",
                                              initpar = c(-0.8415672, -0.2570451, -8.7500000, -4.7500000,  0.5600000, -5.1159958, -2.4651040),
                                              n_cores = 18)
# save(mod_gompertz0.1, file = "mod_gompertz0.1.RData")

mod_lnorm0.1 <- LeadTime::TumorGrowthModel(ExampleData,
                                           d0 = 0.1,
                                           distribution.t0 = "log-normal",
                                           initpar = c(-0.8415672, -0.2570451, -8.7500000, -4.7500000,  0.5600000, 4.0000000, -1.3862944),
                                           n_cores = 18)
# save(mod_lnorm0.1, file = "mod_lnorm0.1.RData")


load("mod_gompertz0.1.RData")
load("mod_weibull0.1.RData")
load("mod_normal0.1.RData")
load("mod_gamma0.1.RData")
load("mod_lnorm0.1.RData")


aic_vec <- c(mod_gompertz0.1$aic,mod_weibull0.1$aic,mod_normal0.1$aic,mod_gamma0.1$aic,mod_lnorm0.1$aic)
names(aic_vec) <- c("Gompertz", "Weibull", "Normal", "Gamma", "Log-normal")
aic_vec


# IC 95%
c(exp(-mod_normal0.1$par[1]), exp(-mod_normal0.1$par[1]-1.96*sqrt(mod_normal0.1$vcov[1,1])),exp(-mod_normal0.1$par[1]+1.96*sqrt(mod_normal0.1$vcov[1,1])))
c(
  exp(-mod_normal0.1$par[1] - mod_normal0.1$par[2]),
  exp((-mod_normal0.1$par[1] - mod_normal0.1$par[2]) - 1.96 * sqrt(mod_normal0.1$vcov[1,1] + mod_normal0.1$vcov[2,2] + 2 * mod_normal0.1$vcov[1,2])),
  exp((-mod_normal0.1$par[1] - mod_normal0.1$par[2]) + 1.96 * sqrt(mod_normal0.1$vcov[1,1] + mod_normal0.1$vcov[2,2] + 2 * mod_normal0.1$vcov[1,2]))
)
c(exp(mod_normal0.1$par[3]), exp(mod_normal0.1$par[3]-1.96*sqrt(mod_normal0.1$vcov[3,3])),exp(mod_normal0.1$par[3]+1.96*sqrt(mod_normal0.1$vcov[3,3])))
c(mod_normal0.1$par[4], mod_normal0.1$par[4]-1.96*sqrt(mod_normal0.1$vcov[4,4]),mod_normal0.1$par[4]+1.96*sqrt(mod_normal0.1$vcov[4,4]))
c(mod_normal0.1$par[5], mod_normal0.1$par[5]-1.96*sqrt(mod_normal0.1$vcov[5,5]),mod_normal0.1$par[5]+1.96*sqrt(mod_normal0.1$vcov[5,5]))
c(exp(mod_normal0.1$par[6]), exp(mod_normal0.1$par[6]-1.96*sqrt(mod_normal0.1$vcov[6,6])),exp(mod_normal0.1$par[6]+1.96*sqrt(mod_normal0.1$vcov[6,6])))
c(exp(mod_normal0.1$par[7]), exp(mod_normal0.1$par[7]-1.96*sqrt(mod_normal0.1$vcov[7,7])),exp(mod_normal0.1$par[7]+1.96*sqrt(mod_normal0.1$vcov[7,7])))


qgamma(0.5,2.09,4.36)*log(2) * 365.25

curve(dgamma(x,2.08,2.65)*log(2), xlim = c(0,5))
qgamma(0.5,2.08,2.65)*log(2)
qgamma(0.5,2.08,2.65)*log(2) * 365.25
qgamma(0.5,2.36,3.00)*log(2)
qgamma(0.5,2.36,3.00)*log(2) * 365.25

# Median Doubling Time
qgamma(0.5,exp(-mod_normal0.1$par[1]),exp(-mod_normal0.1$par[1])/exp(mod_normal0.1$par[2]))*log(2)*365.25
library(numDeriv)
mu_est <- mod_normal0.1$par[1:2]
V_est  <- mod_normal0.1$vcov[1:2, 1:2]

f_transformation <- function(theta) {
  shape <- exp(-theta[1])
  rate  <- exp(-theta[1]) / exp(theta[2])
  res   <- qgamma(0.5, shape = shape, rate = rate) * log(2) * 365.25
  return(res)
}

gradient <- jacobian(f_transformation, mu_est)

variance_delta <- gradient %*% V_est %*% t(gradient)
se_delta       <- sqrt(as.numeric(variance_delta))

val_ponctuelle <- f_transformation(mu_est)
ic_inf <- val_ponctuelle - 1.96 * se_delta
ic_sup <- val_ponctuelle + 1.96 * se_delta

cat("Estimation :", val_ponctuelle, "\n")
cat("SE Delta   :", se_delta, "\n")
cat("IC 95%     : [", ic_inf, ";", ic_sup, "]\n")



# Sensitivity 5mm & 10mm
Sensitivity(5,mod_normal0.1$par[4],mod_normal0.1$par[5],0.1)
Sensitivity(10,mod_normal0.1$par[4],mod_normal0.1$par[5],0.1)
mu_est <- mod_normal0.1$par[4:5]
V_est  <- mod_normal0.1$vcov[4:5, 4:5]

Sensitivity <- function(d, p1, p2, d0){
  return(ifelse(d >= d0, 1 / (1 + exp(-(p1 + p2 * d))), 0))
}

f_sens_5 <- function(p) {
  Sensitivity(d = 5, p1 = p[1], p2 = p[2], d0 = 0.1)
}
grad_5  <- jacobian(f_sens_5, mu_est)
se_5    <- sqrt(as.numeric(grad_5 %*% V_est %*% t(grad_5)))
val_5   <- f_sens_5(mu_est)

f_sens_10 <- function(p) {
  Sensitivity(d = 10, p1 = p[1], p2 = p[2], d0 = 0.1)
}
grad_10 <- jacobian(f_sens_10, mu_est)
se_10   <- sqrt(as.numeric(grad_10 %*% V_est %*% t(grad_10)))
val_10  <- f_sens_10(mu_est)

resultats <- data.frame(
  Diametre = c("5 mm", "10 mm"),
  Sensibilite = c(val_5, val_10),
  SE = c(se_5, se_10),
  IC_inf = c(val_5 - 1.96 * se_5, val_10 - 1.96 * se_10),
  IC_sup = c(val_5 + 1.96 * se_5, val_10 + 1.96 * se_10)
)
print(resultats)


# Cumulative age at onset at 50 and 75
cdf_t0(50,exp(mod_normal0.1$par[6]),exp(mod_normal0.1$par[7]),"normal")
1 - cdf_t0(75,exp(mod_normal0.1$par[6]),exp(mod_normal0.1$par[7]),"normal")
mu_est <- mod_normal0.1$par[6:7]
V_est  <- mod_normal0.1$vcov[6:7, 6:7]

f_mesure_1 <- function(p) {
  res <- cdf_t0(50, exp(p[1]), exp(p[2]), "normal")
  return(res)
}
grad_1 <- jacobian(f_mesure_1, mu_est)
se_1   <- sqrt(as.numeric(grad_1 %*% V_est %*% t(grad_1)))
val_1  <- f_mesure_1(mu_est)

f_mesure_2 <- function(p) {
  res <- 1 - cdf_t0(75, exp(p[1]), exp(p[2]), "normal")
  return(res)
}
grad_2 <- jacobian(f_mesure_2, mu_est)
se_2   <- sqrt(as.numeric(grad_2 %*% V_est %*% t(grad_2)))
val_2  <- f_mesure_2(mu_est)

resultats_cdf <- data.frame(
  Mesure = c("CDF à 50 ans", "1 - CDF à 75 ans"),
  Estimation = c(val_1, val_2),
  SE = c(se_1, se_2),
  IC_inf = c(val_1 - 1.96 * se_1, val_2 - 1.96 * se_2),
  IC_sup = c(val_1 + 1.96 * se_1, val_2 + 1.96 * se_2)
)
print(resultats_cdf)



# Median Volume at symptoms
root_fun <- function(d){
  return(cdf_dsym(d,exp(-mod_normal0.1$par[1]),exp(-mod_normal0.1$par[1])/exp(mod_normal0.1$par[2]),exp(mod_normal0.1$p[3]),0.1)-0.5)
}
med_vsym <- uniroot(root_fun,c(20,30))$root
mu_est <- mod_normal0.1$par[1:3]
V_est  <- mod_normal0.1$vcov[1:3, 1:3]

get_median_diameter <- function(p) {

  root_fun <- function(d) {
    shape <- exp(-p[1])
    rate  <- exp(-p[1]) / exp(p[2])
    param3 <- exp(p[3])

    val <- cdf_dsym(d, shape, rate, param3, 0.1) - 0.5
    return(val)
  }

  res_uniroot <- uniroot(root_fun, c(20, 30))

  return(res_uniroot$root)
}

gradient <- jacobian(get_median_diameter, mu_est)

variance_delta <- gradient %*% V_est %*% t(gradient)
se_delta       <- sqrt(as.numeric(variance_delta))

val_ponctuelle <- get_median_diameter(mu_est)
ic_inf <- val_ponctuelle - 1.96 * se_delta
ic_sup <- val_ponctuelle + 1.96 * se_delta

cat("Diamètre médian au diagnostic (mm) :", val_ponctuelle, "\n")
cat("SE Delta                            :", se_delta, "\n")
cat("IC 95%                              : [", ic_inf, ";", ic_sup, "]\n")


#-------------------------------------------------------------------------------
# Graph
#-------------------------------------------------------------------------------
load("mod_normal0.1.RData")
library(ggplot2)
library(patchwork)
p1 <- ggplot() +
  xlim(0, 3.5) + # Définir les limites de l'axe X
  stat_function(fun = function(x) dgamma(x/log(2),exp(-mod_normal0.1$par[1]),exp(-mod_normal0.1$par[1])/exp(mod_normal0.1$par[2]))/log(2), size = 1) +
  geom_vline(xintercept = qgamma(0.5,exp(-mod_normal0.1$par[1]),exp(-mod_normal0.1$par[1])/exp(mod_normal0.1$par[2]))*log(2), linetype = "dashed", color = "black", size = 0.8) +
  ylab("Density") +
  xlab("Volume doubling time (years)") +
  theme_minimal()+
  annotate("text", x = qgamma(0.5,exp(-mod_normal0.1$par[1]),exp(-mod_normal0.1$par[1])/exp(mod_normal0.1$par[2]))*log(2) + 1.2, y = 0.6, label = paste("Median =", round(qgamma(0.5,exp(-mod_normal0.1$par[1]),exp(-mod_normal0.1$par[1])/exp(mod_normal0.1$par[2]))*log(2), 2),"years"), color = "black")
p1

source("C:/Users/mariu/Desktop/full_save_10_04_25/pckLTcor/LeadTime/R/ScreeningSensitivity.R", echo=TRUE)
p2 <- ggplot() +
  xlim(0.1, 20) + # Définir les limites de l'axe X
  stat_function(fun = function(x) Sensitivity(x,mod_normal0.1$par[4],mod_normal0.1$par[5],0.1), size = 1) +
  ylab("Screening sensitivity") +
  xlab("Diameter (mm)") +
  theme_minimal()
p2

p3 <- ggplot() +
  xlim(0, 100) + # Définir les limites de l'axe X
  stat_function(fun = function(x) pdf_t0(x,exp(mod_normal0.1$par[6]),exp(mod_normal0.1$par[7]),"normal") , size = 1) +
  ylab("Density") +
  xlab("Age at tumor onset") +
  theme_minimal()
p3


med_vsym <- uniroot(root_fun,c(20,30))$root
p4 <- ggplot() +
  xlim(0, 90) + # Définir les limites de l'axe X
  stat_function(fun = function(x) pdf_dsym(x,exp(-mod_normal0.1$par[1]),exp(-mod_normal0.1$par[1])/exp(mod_normal0.1$par[2]),exp(mod_normal0.1$par[3]),0.5), size = 1) +
  geom_vline(xintercept = med_vsym, linetype = "dashed", color = "black", size = 0.8) +
  ylab("Density") +
  xlab("Diameter at symptom onset (mm)") +
  theme_minimal() +
  annotate("text", x = med_vsym + 26, y = 0.03, label = paste("Median =", round(med_vsym, 2),"mm"), color = "black")

p4

(p1 + p2) / (p3 + p4)  # & theme(legend.position = "bottom")

#-------------------------------------------------------------------------------
# Lead Time estimation
#-------------------------------------------------------------------------------
library(survminer)
load("ExampleData.RData")
load("mod_normal0.1.RData")

cor_normal_0.1_cancer_deaths <- LeadTime::LeadTimeCor(ExampleData,
                                                      model = mod_normal0.1,
                                                      overdiagnosis = FALSE,
                                                      d0 = 0.1,
                                                      distribution.t0 = "normal",
                                                      n_cores = 18)
# save(cor_normal_0.1_cancer_deaths, file = "cor_normal_0.1_cancer_deaths.RData")
load("cor_normal_0.1_cancer_deaths.RData")

summary(cor_normal_2009_0.1_cancer_deaths$survival$scr, times = c(2,5))
summary(cor_normal_2009_0.1_cancer_deaths$survival$cor, times = c(2,5))
summary(cor_normal_2009_0.1_cancer_deaths$survival$sym, times = c(2,5))

sum(cor_normal_0.1_cancer_deaths$overdetected)
sum(cor_normal_0.1_cancer_deaths$overdiagnosis)



df <- data.frame(id = cor_normal_0.1_cancer_deaths$id,
                 detection = cor_normal_0.1_cancer_deaths$detection,
                 observed = cor_normal_0.1_cancer_deaths$observed,
                 leadtime = cor_normal_0.1_cancer_deaths$leadtime,
                 death = cor_normal_0.1_cancer_deaths$death)
sym <- subset(df, detection == 1)
scr <- subset(df, detection == 2)

survival <- list(
  scr_obs    = survival::survfit(survival::Surv(scr$observed, scr$death) ~ 1),
  scr_lea    = survival::survfit(survival::Surv(ifelse(scr$observed - scr$leadtime > 0, scr$observed - scr$leadtime, 0), scr$death) ~ 1),
  sym_obs    = survival::survfit(survival::Surv(sym$observed, sym$death) ~ 1)
)

leg <- c("Observed screen-detected", "Corrected lead time screen-detected", "Symptomatically detected")
pal <- c("violetred4", "dodgerblue3", "aquamarine3")


ggsurv <- survminer::ggsurvplot_combine(
  fit = survival,
  censor = FALSE,
  xlab = "Years from diagnosis",
  ylab = "Survival probability",
  break.time.by = 1,
  conf.int = TRUE,
  xlim = c(0,5),
  ylim = c(0.75,1),
  legend.labs = leg,
  palette = pal,
  ggtheme = ggplot2::theme_bw()
)

ggsurv$plot <- ggsurv$plot +
  ggplot2::labs(color = NULL, fill = NULL)

print(ggsurv)


hist(scr$leadtime, 30)
