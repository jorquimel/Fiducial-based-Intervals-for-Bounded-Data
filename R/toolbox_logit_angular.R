#------------------------------------------------------------------------------#
# Program Name   : toolbox for logit and angular transformation
# Author         : Jorge Quiroz
# Usage          : Bounded data 0-1
#                : Logit-normal
#                : Asin-sqrt-normal
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                 Normal distribution   Prediction interval
#-------------------------------------------------------------------------------

normal_prediction_interval <- function(x, conf.level = 0.95) {
  n <- length(x)
  mean_x <- mean(x)
  sd_x <- sd(x)
  alpha <- 1 - conf.level
  t_crit <- qt(1 - alpha/2, df = n - 1)
  margin <- t_crit * sd_x * sqrt(1 + 1/n)
  lower <- mean_x - margin
  upper <- mean_x + margin
  c(lower = lower, upper = upper)
}

#------------------------------------------------------------------------------#
#                           Angular transformation
#------------------------------------------------------------------------------#
asin_tranform <- function(x){
  if (any(x <= 0 | x >= 1)) stop ("x must be in (0, 1)")
  y <- 2*asin(sqrt(x))
  return(as.numeric(y))
}

asin_back <- function(x){
  y <- (sin(x/2))^2
  return(as.numeric(y))
}

#------------------------------------------------------------------------------#
#                           End of file
#------------------------------------------------------------------------------#
