#' Print of tumor growth model
#'
#' @param x A TumorGrowth model class object.
#' @param ... Other parameters to be passed through to printing functions. (Ignored in this function)
#'
#' @importFrom stats pnorm qnorm
#'
#' @export
print.TumorGrowthModel <- function(x,...){
  if(x$convergence == FALSE){stop("Convergence = FALSE")}
  par <- x$par
  sd <- rep(NA, length(par))
  sd[!x$fix] <- sqrt(diag(x$vcov))
  lower <- x$par - stats::qnorm(0.975)*sd
  upper <- x$par + stats::qnorm(0.975)*sd
  p <- 2 * (1 - stats::pnorm(abs(as.numeric(x$par) / sd)))
  sig <- ifelse(!is.na(p),paste0(ifelse(p <= 0.1 & p > 0.05, '.', ''), ifelse(p <= 0.05, '*', ''),
                ifelse(p <= 0.01, '*', ''), ifelse(p <= 0.001, '*', '')),NA)
  p <- ifelse(p < 0.001,"<0.001",round(p,3))
  summ <- data.frame(cbind(round(par,3),round(lower,3),round(upper,3),p,sig))

  summary <- rbind(c("","","","",""),summ[1:length(x$growth),],c("","","","",""),summ[(length(x$growth)+1):(length(x$growth)+length(x$symptoms)),],c("","","","",""),summ[(length(x$growth)+length(x$symptoms)+1):((length(x$growth)+length(x$symptoms))+length(x$sensitivity)),],c("","","","",""),summ[((length(x$growth)+length(x$symptoms))+length(x$sensitivity)+1):((length(x$growth)+length(x$symptoms))+length(x$sensitivity)+length(x$onset)),])
  parnames <- c("Inverse growth rate  ",paste(" ",x$growth),"Volume at symptoms  ", paste(" ",x$symptoms),"Screening sensitivity  ", paste(" ",x$sensitivity),"Age at tumor onset  ", paste(" ",x$onset))

  names <- parnames
  for(i in 1:length(parnames)){
    for(j in nchar(parnames[i]):max(nchar(parnames))+1){names[i] <- paste(names[i],"")}
  }

  output <- cbind(names, summary)
  colnames(output) <- c("parameters     ","estimate", "lower", "upper", "p-value","")
  cat("---------------Continuous tumor growth model---------------\n\n")
  print(output, row.names = FALSE, quote = FALSE)
  cat("\nModel successfuly estimated with log-likelihood :", x$lld, "\n")
}



#' Print of a lead time estimation
#'
#' @param x A LeadTime class object.
#' @param ... Other parameters to be passed through to printing functions. (Ignored in this function)
#'
#' @export
print.LeadTime <- function(x,...){
  plot(x)
  if(!is.null(x$model)){print(x$model)}
  cat("\n\n---------------Lead time correction---------------")
  if(length(x$detection[x$detection == 1]) > 1){
    cat("\n\nSymptomatically detected :", sum(ifelse(x$detection == 1,1,0)))
    cat("\n   Deaths :", sum(ifelse(x$detection == 1 & x$death == 1,1,0)))
    cat("\n   Average follow-up time after diagnosis :", round(sum(x$observed[x$detection == 1])/sum(ifelse(x$detection == 1,1,0)),2))
  }
  cat("\n\nScreen detected :", sum(ifelse(x$detection == 2,1,0)))
  cat("\n   Deaths :", sum(ifelse(x$detection == 2 & x$death == 1,1,0)))
  cat("\n   Average follow-up time after diagnosis :", round(sum(x$observed[x$detection == 2])/sum(ifelse(x$detection == 2,1,0)),2))
  cat("\n   Number of full corrected patients (lead time >= post-diagnosis follow-up):", sum(ifelse(x$corrected == 0 & x$detection == 2,1,0)))

  cat("\n   Mean of the lead time (correction) :", round(sum(x$leadtime[x$detection == 2])/sum(ifelse(x$detection == 2,1,0)),2))
  # cat("\n   Relative part of lead time over the follow-up time:", round(sum(100*x$leadtime[x$detection == 2])/sum(x$observed[x$detection == 2]),2),"%")
}
