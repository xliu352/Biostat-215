arcsin.ci <- function(fit, alpha = 0.05, t = NULL) {
  
  #- Extract necessary values from survfit object
  surv.est  <- fit$surv
  sigma     <- fit$std.err #This is sigma, not the SE of S(t)
  surv.time <- fit$time
  
  #- Get limits element wise
  lims <- data.frame(time = fit$time, surv = surv.est,  sigma = fit$std.err)
  
  for(i in 1:length(surv.est)) {
  lims$LL[i] <- sin(max(0, asin(sqrt(surv.est[i])) - 0.5 * qnorm(1 - alpha / 2) * sigma[i] * 
                  sqrt(surv.est[i] / (1 - surv.est[i]))))^2
  
  lims$UL[i] <- sin(min(pi / 2, asin(sqrt(surv.est[i])) +  0.5 * qnorm(1 - alpha / 2) * sigma[i] * 
                  sqrt(surv.est[i] / (1 - surv.est[i]))))^2
  }
  
  lims <- round(lims, 4)
  
  return(lims)
}