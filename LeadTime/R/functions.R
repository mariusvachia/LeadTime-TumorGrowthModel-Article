#-------------------------------------------------------------------------------
# Density
#-------------------------------------------------------------------------------
pdf_t0 <- function(t0, p1, p2, distribution.t0 = "normal"){
  if(distribution.t0 == "weibull"){
    return(dweibull(t0,p1,p2)) # shape = 6.5, scale = 70
  }else if(distribution.t0 == "gompertz"){
    return(p2*p1*exp(p1+p2*t0-p1*exp(p2*t0))) # shape = 0.0015, scale = 0.95
  }else if(distribution.t0 == "gamma"){
    return(dgamma(t0,p1,p2)) # shape = 30, rate = 0.55
  }else if(distribution.t0 == "log-normal"){
    return(dlnorm(t0,p1,p2)) # log(mean) = 4, log(sd) = 0.18 #
  }else if(distribution.t0 == "normal"){
    return(dnorm(t0,p1,p2)/(1-pnorm(0,p1,p2))) # Truncated before 0
  }
}

pdf_r <- function(r,p1,p2){
  return(dgamma(r,p1,p2))
}

cpdf_vsym <- function(v,r,eta,v0){
  return(dexp(v - v0,r*eta))
}

pdf_vsym <- function(v, p1, p2, p3, v0){
    return(p3*p1*p2^p1/(p2 + p3*(v - v0))^(p1+1))
}

pdf_dsym <- function(d, p1, p2, p3, d0){
  return(pi/2*d^2*pdf_vsym(Volume(d),p1,p2,p3,Volume(d0)))
}

# cpdf_tsym <- function(t,r,t0,eta,v0){
#   return(eta*v0*exp(-eta*r*(v0*exp((t-t0)/r)-v0)+(t-t0)/r))
# }


#-------------------------------------------------------------------------------
# Cumulative
#-------------------------------------------------------------------------------
cdf_t0 <- function(t0, p1, p2, distribution.t0 = "normal"){
  if(distribution.t0 == "weibull"){
    return(pweibull(t0,p1,p2)) # shape = 6.5, scale = 70
  }else if(distribution.t0 == "gompertz"){
    return(1 - exp(-p1*(exp(p2*t0)-1))) # shape = 0.0015, scale = 0.95
  }else if(distribution.t0 == "gamma"){
    return(pgamma(t0,p1,p2)) # shape = 30, rate = 0.55
  }else if(distribution.t0 == "log-normal"){
    return(plnorm(t0,p1,p2)) # log(mean) = 4, log(sd) = 0.18 #
  }else if(distribution.t0 == "normal"){
    return(pnorm(t0,p1,p2)/(1-pnorm(0,p1,p2))) # Truncated before 0
  }
}


cdf_r <- function(r,p1,p2){
  return(pgamma(r,p1,p2))
}

ccdf_vsym <- function(v,r,eta,v0){
  return(pexp(v - v0,r*eta))
}

cdf_vsym <- function(v, p1, p2, p3, v0){
  return(1 - (p2/(p2+p3*(v - v0)))^p1)
}
cdf_vsym <- Vectorize(cdf_vsym,"v")

cdf_dsym <- function(d, p1, p2, p3, d0){
  return(cdf_vsym(Volume(d),p1,p2,p3,Volume(d0)))
}
