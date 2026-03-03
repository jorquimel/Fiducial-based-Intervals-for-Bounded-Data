#------------------------------------------------------------------------------#
# Program Name   : toolbox_kumaraswamy GPQ
# Author         : Jorge Quiroz
# Usage          : Pivots for alpha and GPQ for beta of
#                : Kumaraswamy distribution
#------------------------------------------------------------------------------#
library(extraDistr)
library(rootSolve)
library(numbers)
library(greybox)

#------------------------------------------------------------------------------#
#             Exact pivotal quantity for parameter alpha
#-------------------------------------------------------------------------------
pivot_a <- function(a, x_data, target,
                    type = c("greenwood", "fdist", "wang"),
                    zero = c("linear", "quadratic")){
  type <- match.arg(type)
  zero <- match.arg(zero)
  n <- length(x_data)

  x_data <- sort(x_data)
  logF <- log(1 - x_data^a)


  if(type == "greenwood"){
    pivot <- sum(logF^2)/(sum(logF))^2
  }
  if(type == "fdist"){
    r1 <- floor(n/2);  r2 <- n - r1
    num <- sum(logF[1:r1]) + r2*logF[r1]
    den <- sum(logF[(r1 + 1):n]) - r2*logF[r1]

    pivot <- (r2*num)/(r1*den)
  }
  if(type == "wang"){
    cf <- n - 1:n
    si <- cumsum(logF) + cf*logF
    si <- si[-n]
    stotal <- sum(logF)

    pivot <- 2*sum(log(stotal/si))
  }

  if(zero == "linear"){
    return(pivot - target)
  } else {
    return((pivot - target)^2)
  }

}

#------------------------------------------------------------------------------#
#        Generalized pivotal quantity for parameters alpha and beta
#-------------------------------------------------------------------------------
gpq_ab <- function(x_data, nrep = 2e3, bound = c(1e-3, 100),
                   type = c("greenwood", "fdist", "wang"),
                   zero = c("linear", "quadratic")){
  type <- match.arg(type)
  zero <- match.arg(zero)
  lbnd <- bound[1];  ubnd <- bound[2]

  n <- length(x_data)    # number of obs

  # greenwood approach
  if(type == "greenwood"){
    dfa <- n - 1
    rgreenwood <- replicate(nrep, rGreenwood(df = dfa))

    if(zero == "linear"){
      Ra <- sapply(rgreenwood, function(x){uniroot(f = pivot_a,
                                                   interval = c(lbnd, ubnd),
                                                   extendInt = "yes",
                                                   x_data = x_data,
                                                   target = x,
                                                   type = "greenwood", zero = "linear")$root
      })
    }
    if(zero == "quadratic"){
      Ra <- sapply(rgreenwood, function(x){optimize(f = pivot_a,
                                                    interval = c(lbnd, ubnd),
                                                    x_data = x_data,
                                                    target = x,
                                                    type = "greenwood", zero = "quadratic")$minimum
      })
    }
  }

  # fdist approach
  if(type == "fdist"){
    r1 <- floor(n/2);   r2 <- n - r1
    df1 <- 2*r1;  df2 <- 2*r2
    rf_a <- rf(nrep, df1 = df1, df2 = df2)
    if(zero == "linear"){
      Ra <- sapply(rf_a, function(x){uniroot(f = pivot_a,
                                             interval = c(lbnd, ubnd),
                                             extendInt = "yes",
                                             x_data = x_data,
                                             target = x,
                                             type = "fdist", zero = "linear")$root
      })
    }
    if(zero == "quadratic"){
      Ra <- sapply(rf_a, function(x){optimize(f = pivot_a,
                                              interval = c(lbnd, ubnd),
                                              x_data = x_data,
                                              target = x,
                                              type = "fdist", zero = "quadratic")$minimum
      })
    }
  }

  # Wang approach
  if(type == "wang"){
    dfa <- 2*(n - 1)
    rchisq_a <- rchisq(nrep, df = dfa)

    if(zero == "linear"){
      Ra <- sapply(rchisq_a, function(x){uniroot(f = pivot_a,
                                                 interval = c(lbnd, ubnd),
                                                 extendInt = "yes",
                                                 x_data = x_data,
                                                 target = x,
                                                 type = "wang", zero = "linear")$root
      })
    }
    if(zero == "quadratic"){
      Ra <- sapply(rchisq_a, function(x){optimize(f = pivot_a,
                                                  interval = c(lbnd, ubnd),
                                                  x_data = x_data,
                                                  target = x,
                                                  type = "wang", zero = "quadratic")$minimum
      })
    }
  }

  #  Generate values for Rb
  dfb <- 2*n
  r_chisq <- rchisq(nrep, df = dfb)

  Rb <- mapply(FUN = function(chsq, a){
    b <- chsq/(-2*sum(log(1 - x_data^a)))
  }, r_chisq, Ra)

  # return gpq
  return(data.frame(Ra = Ra, Rb = Rb))
}

#------------------------------------------------------------------------------#
#              Kumaraswamy parameters:  beta, mean, quantile
#-------------------------------------------------------------------------------
kumar_mean <- function(a, b){
  mk <- b*beta(1 + 1/a, b)
  return(mk)
}

kumar_quantile <- function(a, b, P = 0.95, side = c("lower", "upper")){
  side <- match.arg(side)
  Pu <- (1 + P)/2
  Pl <- 1 - Pu

  if(side == "lower"){
    Qp <- (1 - (1 - Pl)^(1/b))^(1/a)
  } else {
    Qp <- (1 - (1 - Pu)^(1/b))^(1/a)
  }
  return(Qp)
}

kumar_parameter <- function(df, P = 0.95, conf = c(0.025, 0.05, 0.50, 0.95, 0.975)){
  # CI for parameter b
  conf_int_b <- quantile(df$Rb, probs = conf)

  # CI for mean
  gpq_mean <- mapply(FUN = kumar_mean, a = df$Ra, b = df$Rb)
  conf_int_mean <- quantile(gpq_mean, probs = conf)

  # two-sided CI for Quantile
  tilow <- mapply(FUN = kumar_quantile,
                  a = df$Ra, b = df$Rb, side = "lower", P = P)

  tiup <- mapply(FUN = kumar_quantile,
                 a = df$Ra, b = df$Rb, side = "upper", P = P)

  tolint <- c(quantile(tilow, probs = c(0.025, 0.05)),
              quantile(tiup, probs = c(0.95, 0.975)))

  return(list(conf_int_b = conf_int_b,
              conf_int_mean = conf_int_mean,
              tolint = tolint))
}

#--------------------------------------------------------------------------
#                 Calculate Fiducial prediction interval
#--------------------------------------------------------------------------
predict_interval <- function(df, clevel = 0.95){
  # Calculate fiducial prediction intervals using the
  # GPQ distribution provided as dataframe df

  nrun <- nrow(df)
  Ru <- runif(nrun)

  FQ <- (1 - Ru^(1/df$Rb))^(1/df$Ra)
  pred_int <- quantile(FQ, probs = clevel)
  #print(pred_int)
  return(pred_int)
}

predict_interval_two <- function(df, clevel = 0.95){
  # Calculate fiducial prediction intervals using the
  # GPQ distribution provided as dataframe df
  uplevel <- (1 + clevel)/2
  lowlevel <- 1 - uplevel
  nrun <- nrow(df)
  Ru <- runif(nrun)

  FQ <- (1 - Ru^(1/df$Rb))^(1/df$Ra)
  pred_int <- quantile(FQ, probs = c(lowlevel, uplevel))
  #print(pred_int)
  return(pred_int)
}

#------------------------------------------------------------------------------#
#               Functions to generate greenwood statistics
#-------------------------------------------------------------------------------
greenwood.dist <- function(df = 5){
  #
  # generate random values from greenwood distribution
  # df: number of ordered uniform variable (usually: sample size - 1)
  #
  ordered_unif <- c(0, sort(runif(df)), 1)
  space_unif <- diff(ordered_unif)
  rgw <- sum( space_unif^2 )/(sum( space_unif ))^2
  return(rgw)
}

rGreenwood <- function(df = 5){
  #
  # generate random values from greenwood statistic
  # df: number of ordered uniform variable (Use: sample size - 1)
  #
  ordered_unif <- c(0, sort(runif(df)), 1)
  space_unif <- diff(ordered_unif)
  rgw <- sum( space_unif^2 )
  return(rgw)
}

#------------------------------------------------------------------------------#
#                            End of file
#------------------------------------------------------------------------------#
