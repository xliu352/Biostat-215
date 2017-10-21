bc.median.ci <- function(fit, type = "linear", alpha = 0.05, table = FALSE) {

  if((min(fit$surv) < 0.5) == FALSE) {
    stop("Estimate for median survival time does not exist. Can not compute
         confidence interval.")
  }
  
  median.est <- min(fit$time[fit$surv == max(fit$surv[fit$surv - 0.5 <= 0])])
  
  surv.est <- fit$surv
  se.surv  <- fit$std.err * surv.est
  
  if (type != "linear" & type != "log" & type != "asin"){
     stop("type must be one of the three options: linear, log, asin")
  }

  if (type == "linear") {
    zScore <- (surv.est - (1 - 0.5)) / se.surv
  } else if (type == "log") {
    zScore <- (log(-log(surv.est)) - log(-log(1 - 0.5))) * (surv.est * log(surv.est)) /
      se.surv
  } else if (type == "asin") {
    zScore <- 2 * (asin(sqrt(surv.est)) - asin(sqrt(1 - 0.5))) * sqrt(surv.est * (1 - surv.est)) /
      se.surv
  }
  res <- data.frame(time = fit$time, z = zScore)
 
  #- Finding limits
  LL <- ifelse(max(zScore) < qnorm(1 - alpha / 2),
         NA, max(min(fit$time[zScore <  qnorm(1 - alpha / 2)])))
  
  UL <- ifelse(min(zScore) > -qnorm(1 - alpha / 2),
               NA, min(max(fit$time[zScore > -qnorm(1 - alpha / 2)])))
  
  results <- list()
  results$median <- round(median.est, 3)
  results$lower  <- ifelse(is.na(LL), NA, round(LL, 3))
  results$upper  <- ifelse(is.na(UL), NA, round(UL, 3))
  results$type   <- type
  results$alpha  <- alpha
  if(table == TRUE) results$table  <- res
  
  return(results)
}

