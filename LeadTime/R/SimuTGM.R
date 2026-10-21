#' Data simulation for the tumor growth model
#'
#' @param n The number of individuals in the simulated sample.
#' @param p The list of parameters used for the simulation. The argument is in
#' list form so that a single value can be used for a parameter or a vector of
#' size n with linear predictors.
#' @param d0 The smallest diameter of a tumor detectable by screening defining
#' the starting point of the tumor growth function. Default is 0.5
#' @param distribution.t0 A str term to define the distribution of age at onset.
#' @param screening_start The minimum age for screening.
#' @param screening_stop The maximum age for screening
#' @param screening_interval The time interval between two screenings.
#' @param prob_program The proportion of individual attending a screening program
#' @param prob_visit The probability for a patient to attend each visit.
#' @param nor1 The mean parameter (in years) of the Gaussian noise added to the screening interval.
#' @param nor2 The standard deviation parameter of the Gaussian noise to the screening interval.
#' @param surv_size Survival regression coefficient associated with tumour size,
#' modulating the rate of the exponential distribution to simulate survival time following diagnosis.
#' @param surv_growth Survival regression coefficient associated with the inverse tumour growth rate,
#' modulating the rate of the exponential distribution to simulate survival time following diagnosis.
#'
#' @importFrom stats qgamma rbinom runif rweibull rnorm rlnorm rexp
#'
#' @returns A dataframe with n observations and 11 variables.
#' \item{r}{The tumor growth rate}
#' \item{vsym}{The tumor volume at symptoms}
#' \item{dsym}{The tumor diameter at symptoms}
#' \item{t0}{The age at tumor onset}
#' \item{tsym}{The age at symptoms}
#' \item{det}{The detection mode}
#' \item{vdet}{The tumor volume at detection}
#' \item{ddet}{The tumor diameter at detection}
#' \item{tdet}{The age at tumor detection}
#' \item{scr}{The screening results}
#' \item{tscr}{The patient age at screening results}
#' @export
SimuTGM <- function(n,
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
                    surv_growth = 4e-3){
  data_simu <- data.frame(id = c(1:n))
  #-----------------------------------------------------------------------------
  # natural history of cancer
  #-----------------------------------------------------------------------------
  v0 <- Volume(d0)
  data_simu$r <- qgamma(runif(n,0,1),p$p1,rate = p$p2)

  data_simu$vsym <- v0-log(runif(n,0,1))/(p$p3*data_simu$r)
  data_simu$dsym <- Diameter(data_simu$vsym)
  if(distribution.t0 == "weibull"){
    data_simu$t0 <- rweibull(n,p$p6,p$p7)
  }else if(distribution.t0 == "gompertz"){
    data_simu$t0 <- (1/p$p7) * log(1 - (1/p$p6) * log(runif(n,0,1)))
  }else if(distribution.t0 == "gamma"){
    data_simu$t0 <- qgamma(runif(n,0,1),p$p6,rate = p$p7)
  }else if(distribution.t0 == "log-normal"){
    data_simu$t0 <- rlnorm(n,p$p6,p$p7)
  }else if(distribution.t0 == "normal"){
    data_simu$t0 <- qnorm(runif(n, pnorm(0,p$p6,p$p7),1),p$p6,p$p7)
  }
  data_simu$tsym <- data_simu$t0 + data_simu$r*log(data_simu$vsym/v0)

  data_simu$program <- rbinom(n,1,prob_program)
  start <- screening_start + runif(n,0,screening_interval)

  data_simu$det <- numeric(n)
  data_simu$vdet <- data_simu$vsym
  data_simu$ddet <- data_simu$dsym
  data_simu$tdet <- data_simu$tsym

  #-----------------------------------------------------------------------------
  # screening and cancer detection
  #-----------------------------------------------------------------------------
  tscr <- list()
  scr <- list()
  for(i in 1:n){
    # print(i)
    if(data_simu$tsym[i] > start[i] & data_simu$program[i] & start[i] <= screening_stop){

      tscr[[i]] <- start[i] + cumsum(c(0,rnorm(max(round((min(data_simu$tsym[i],screening_stop) - start[i])/screening_interval),0),screening_interval+nor1,nor2)))

      vis <- rbinom(length(tscr[[i]]),1,prob_visit)
      tscr[[i]] <- tscr[[i]][vis == 1 & data_simu$tsym[i] > tscr[[i]] & tscr[[i]] <= screening_stop]
      dscr <- Diameter(v0*exp((tscr[[i]]-data_simu$t0[i])/data_simu$r[i]))
      pscr <- Sensitivity(dscr,p$p4,p$p5,d0)
      scr[[i]] <- rbinom(length(pscr),1,pscr)
    }else{
      tscr[[i]] <- numeric(0)
      scr[[i]] <- numeric(0)
    }
    if(sum(scr[[i]][tscr[[i]]>data_simu$t0[i]]) > 0){
      det <- min(which(tscr[[i]]>data_simu$t0[i] & scr[[i]]==1))
      scr[[i]] <- scr[[i]][1:det]
      tscr[[i]] <- tscr[[i]][1:det]

      data_simu$det[i] <- 1
      data_simu$tdet[i] <- tscr[[i]][det]
      data_simu$vdet[i] <- v0*exp((data_simu$tdet[i]-data_simu$t0[i])/data_simu$r[i])
      data_simu$ddet[i] <- Diameter(data_simu$vdet[i])
    }
  }
  data_simu$det <- data_simu$det + 1
  data_simu$scr <- scr
  data_simu$tscr <- tscr

  #-----------------------------------------------------------------------------
  # Lead time
  #-----------------------------------------------------------------------------
  data_simu$undet_time <- data_simu$tdet - data_simu$t0
  data_simu$lead_time <- data_simu$tsym - data_simu$tdet
  data_simu$length_time <- data_simu$tsym - data_simu$t0

  #-----------------------------------------------------------------------------
  # Survival
  #-----------------------------------------------------------------------------
  nscr <- sum(data_simu$det==2)
  data_simu$death <- 1
  data_simu$tdeath_w0 <- data_simu$tsym + rexp(n, data_simu$dsym * surv_size + surv_growth / data_simu$r)
  data_simu$tdeath_w1 <- data_simu$tdeath_w0
  data_simu$tdeath_w1[data_simu$det == 2] <- data_simu$tsym[data_simu$det == 2] + rexp(nscr, data_simu$ddet[data_simu$det == 2] * surv_size + surv_growth / data_simu$r[data_simu$det == 2])

  data_simu$surv_w0 <- data_simu$tdeath_w0 - data_simu$tsym
  data_simu$surv_w1 <- data_simu$tdeath_w1 - data_simu$tdet
  data_simu$benefit <- data_simu$surv_w1 - data_simu$surv_w0


  #----------------------------------------------------------------------------
  # To avoid bug in Tumor growth model
  #----------------------------------------------------------------------------
  data_simu$vdet <- ifelse(data_simu$ddet < d0+0.001,Volume(d0+0.001),data_simu$vdet)
  data_simu$ddet <- ifelse(data_simu$ddet < d0+0.001,d0+0.001,data_simu$ddet)
  return(data_simu)
}


# remove.packages("LeadTime")
# unlink("C:/Users/mariu/AppData/Local/R/win-library/4.5/LeadTime", recursive = TRUE, force = TRUE)


# set.seed(1)
# n <- 2000
# p = list(p1 = 2.32, p2 = 3.00, p3 = exp(-8.75), p4 = -4.75, p5 = 0.56, p6 = 50, p7 = 12, p8 = 0.8)
# p <- list(p1 = 2.08, p2 = 2.65, p3 = exp(-9.30), p4 = -4.36, p5 = 0.54, p6 = 50.97, p7 = 12.52)
# data_simu <- SimuTGM(n,
#                        p,
#                        d0 = 0.1,
#                        distribution.t0 = "normal",
#                        screening_start = 50,
#                        screening_stop = 75,
#                        screening_interval = 2,
#                        prob_program = 0.75,
#                        prob_visit = 0.75,
#                        nor1 = 0.2,
#                        nor2 = 0.15,
#                        surv_size = 6e-4,
#                        surv_growth = 4e-3)
