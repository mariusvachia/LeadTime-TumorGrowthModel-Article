library(LeadTime)
#-------------------------------------------------------------------------------
# Weibull
#-------------------------------------------------------------------------------
n <- 2000
p = list(p1 = 2.32, p2 = 3.00, p3 = exp(-8.75), p4 = -4.75, p5 = 0.56, p6 = 5, p7 = 60)
tau1 <- p$p1
tau2 <- p$p2
eta <- p$p3
beta0 <- p$p4
beta1 <- p$p5
k <- p$p6
lambda <- p$p7
par <- c(-log(tau1),log(tau1)-log(tau2),log(eta),beta0,beta1,log(k),log(lambda))
#-------------------------------------------------------------------------------
# Data simulation
#-------------------------------------------------------------------------------
load(paste0("seeds/seed_weibull.RData"))
for(i in 1:100){
  set.seed(seed_weibull[i])
  data_simu <- SimuTGM(n,
                       p,
                       d0 = 0.5,
                       distribution.t0 = "weibull",
                       screening_start = 50,
                       screening_stop = 75,
                       screening_interval = 2,
                       prob_program = 1,
                       prob_visit = 1,
                       nor1 = 0.2,
                       nor2 = 0.15,
                       surv_size = 6e-4,
                       surv_growth = 4e-3)

  # save(data_simu, file = paste0(path,"weibull/sim_",i,".RData"))
}
#-------------------------------------------------------------------------------
# Gompertz
#-------------------------------------------------------------------------------
n <- 2000
p = list(p1 = 2.32, p2 = 3.00, p3 = exp(-8.75), p4 = -4.75, p5 = 0.56, p6 = 0.006, p7 = 0.085)
tau1 <- p$p1
tau2 <- p$p2
eta <- p$p3
beta0 <- p$p4
beta1 <- p$p5
a <- p$p6
b <- p$p7
par <- c(-log(tau1),log(tau1)-log(tau2),log(eta),beta0,beta1,log(a),log(b))
#-------------------------------------------------------------------------------
# Data simulation
#-------------------------------------------------------------------------------
load(paste0("seeds/seed_gompertz.RData"))
for(i in 1:100){
  set.seed(seed_gompertz[i])
  data_simu <- SimuTGM(n,
                       p,
                       d0 = 0.5,
                       distribution.t0 = "gompertz",
                       screening_start = 50,
                       screening_stop = 75,
                       screening_interval = 2,
                       prob_program = 1,
                       prob_visit = 1,
                       nor1 = 0.2,
                       nor2 = 0.15,
                       surv_size = 6e-4,
                       surv_growth = 4e-3)

  # save(data_simu, file = paste0(path,"gompertz/sim_",i,".RData"))
}
#-------------------------------------------------------------------------------
# Gamma
#-------------------------------------------------------------------------------
n <- 2000
p = list(p1 = 2.32, p2 = 3.00, p3 = exp(-8.75), p4 = -4.75, p5 = 0.56, p6 = 16, p7 = 0.3)
tau1 <- p$p1
tau2 <- p$p2
eta <- p$p3
beta0 <- p$p4
beta1 <- p$p5
a <- p$p6
b <- p$p7
par <- c(-log(tau1),log(tau1)-log(tau2),log(eta),beta0,beta1,-log(a),log(a)-log(b))
#-------------------------------------------------------------------------------
# Data simulation
#-------------------------------------------------------------------------------
load(paste0("seeds/seed_gamma.RData"))
for(i in 1:100){
  set.seed(seed_gamma[i])
  data_simu <- SimuTGM(n,
                       p,
                       d0 = 0.5,
                       distribution.t0 = "gamma",
                       screening_start = 50,
                       screening_stop = 75,
                       screening_interval = 2,
                       prob_program = 1,
                       prob_visit = 1,
                       nor1 = 0.2,
                       nor2 = 0.15,
                       surv_size = 6e-4,
                       surv_growth = 4e-3)

  # save(data_simu, file = paste0(path,"gamma/sim_",i,".RData"))
}
#-------------------------------------------------------------------------------
# Log-normal
#-------------------------------------------------------------------------------
n <- 2000
p = list(p1 = 2.32, p2 = 3.00, p3 = exp(-8.75), p4 = -4.75, p5 = 0.56, p6 = 4, p7 = 0.25)
tau1 <- p$p1
tau2 <- p$p2
eta <- p$p3
beta0 <- p$p4
beta1 <- p$p5
a <- p$p6
b <- p$p7
par <- c(-log(tau1),log(tau1)-log(tau2),log(eta),beta0,beta1,a,log(b))
#-------------------------------------------------------------------------------
# Data simulation
#-------------------------------------------------------------------------------
load(paste0("seeds/seed_lnorm.RData"))
for(i in 1:100){
  set.seed(seed_lnorm[i])
  data_simu <- SimuTGM(n,
                       p,
                       d0 = 0.5,
                       distribution.t0 = "log-normal",
                       screening_start = 50,
                       screening_stop = 75,
                       screening_interval = 2,
                       prob_program = 1,
                       prob_visit = 1,
                       nor1 = 0.2,
                       nor2 = 0.15,
                       surv_size = 6e-4,
                       surv_growth = 4e-3)

  # save(data_simu, file = paste0(path,"lnorm/sim_",i,".RData"))
}
#-------------------------------------------------------------------------------
# Normal
#-------------------------------------------------------------------------------
n <- 2000
p = list(p1 = 2.32, p2 = 3.00, p3 = exp(-8.75), p4 = -4.75, p5 = 0.56, p6 = 55, p7 = 11.8)
tau1 <- p$p1
tau2 <- p$p2
eta <- p$p3
beta0 <- p$p4
beta1 <- p$p5
mu <- p$p6
sigma <- p$p7
par <- c(-log(tau1),log(tau1)-log(tau2),log(eta),beta0,beta1,log(mu),log(sigma))
#-------------------------------------------------------------------------------
# Data simulation
#-------------------------------------------------------------------------------
load(paste0("seeds/seed_normal.RData"))
for(i in 1:100){
  set.seed(seed_normal[i])
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

  # save(data_simu, file = paste0(path,"normal/sim_",i,".RData"))
}



curve(pdf_t0(x,5,60,"weibull"), xlim = c(0,100), ylim = c(0,0.05), col = "firebrick3", lwd = 1.5)
curve(pdf_t0(x,0.006,0.085,"gompertz"), xlim = c(0,100), col = "springgreen3", add = TRUE, lwd = 1.5)
curve(pdf_t0(x,16,0.30,"gamma"), xlim = c(0,100), col = "navy", add = TRUE, lwd = 1.5)
curve(pdf_t0(x,4,0.25,"log-normal"), xlim = c(0,100), col = "goldenrod3", add = TRUE, lwd = 1.5)
curve(pdf_t0(x,55,11.8,"normal"), xlim = c(0,100), col = "lightblue3", add = TRUE, lwd = 1.5)





# Graph article
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# 1. Data Generation
x_vals <- seq(0, 100, length.out = 200)
df <- data.frame(
  x = x_vals,
  Weibull      = pdf_t0(x_vals, 5, 60, "weibull"),
  Gompertz     = pdf_t0(x_vals, 0.006, 0.085, "gompertz"),
  Gamma        = pdf_t0(x_vals, 16, 0.30, "gamma"),
  `Log-normal` = pdf_t0(x_vals, 4, 0.25, "log-normal"),
  `Truncated normal` = pdf_t0(x_vals, 55, 11.8, "normal"),
  check.names  = FALSE
) %>%
  pivot_longer(-x, names_to = "Distribution", values_to = "Density")

# Order
target_order <- c("Gompertz", "Weibull", "Truncated normal", "Gamma", "Log-normal")
df$Distribution <- factor(df$Distribution, levels = target_order)

# 2. Background color
df_background <- df %>%
  rename(Distribution_bg = Distribution)

# 3. Labels
df_labels <- data.frame(
  Distribution = factor(target_order, levels = target_order),
  x_pos = 4,
  y_pos = 0.035,
  param_text = c(
    "shape = 0.006\nrate = 0.085",
    "shape = 5\nscale = 60",
    "mean = 55\nsd = 11.8",
    "shape = 16\nrate = 0.30",
    "meanlog = 4\nsdlog = 0.25"
  )
)

# --- Split in two groups ---
target_top <- c("Gompertz", "Weibull", "Truncated normal")
target_bottom <- c("Gamma", "Log-normal")

# ==============================================================================
# GRAPH UPPER (3 panels)
# ==============================================================================
p_top <- ggplot(df %>% filter(Distribution %in% target_top), aes(x = x, y = Density)) +
  geom_line(data = df_background,
            aes(x = x, y = Density, group = Distribution_bg),
            color = "grey93", linewidth = 0.5) +
  geom_line(aes(color = Distribution), linewidth = 1) +
  geom_text(data = df_labels %>% filter(Distribution %in% target_top),
            aes(x = x_pos, y = y_pos, label = param_text),
            hjust = 0, vjust = 1,
            size = 3.2, color = "grey30", fontface = "italic",
            inherit.aes = FALSE) +
  facet_wrap(~Distribution, nrow = 1) +
  labs(x = NULL, y = "Density") +       # No title
  ylim(0, 0.035) +
  scale_color_brewer(palette = "Set1", drop = FALSE) + # drop = FALSE Keep color
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 11)
  )

# ==============================================================================
# GRAPHIQUE LOWER (2 panels)
# ==============================================================================
p_bottom <- ggplot(df %>% filter(Distribution %in% target_bottom), aes(x = x, y = Density)) +
  geom_line(data = df_background,
            aes(x = x, y = Density, group = Distribution_bg),
            color = "grey93", linewidth = 0.5) +
  geom_line(aes(color = Distribution), linewidth = 1) +
  geom_text(data = df_labels %>% filter(Distribution %in% target_bottom),
            aes(x = x_pos, y = y_pos, label = param_text),
            hjust = 0, vjust = 1,
            size = 3.2, color = "grey30", fontface = "italic",
            inherit.aes = FALSE) +
  facet_wrap(~Distribution, nrow = 1) +
  labs(x = "Age at tumor onset", y = "Density") +
  ylim(0, 0.035) +
  scale_color_brewer(palette = "Set1", drop = FALSE) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 11)
  )

# ==============================================================================
# Merge with patchwork
# ==============================================================================
layout_design <- "
  AAAAAA
  #BBBB#
"
p_top + p_bottom + plot_layout(design = layout_design)










