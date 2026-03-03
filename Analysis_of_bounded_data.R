#------------------------------------------------------------------------------#
# Program Name   : Analysis of Bounded Data
# Author         : Jorge Quiroz
# Usage          : Analysis of Continuous data between 0 and 1
#                : Logit-normal
#                : Asin-sqrt-normal
#                : Kumaraswamy: Three pivotals for alpha + GPQ for beta
#------------------------------------------------------------------------------#
options(digits=7)
start.time <- proc.time()[3]

#------------------------------------------------------------------------------#
#                         Set working directory
#-------------------------------------------------------------------------------
# Set working directory to the location of this script.
# If running interactively, set it manually to the code_repo_new folder.
this_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) getwd()
)
setwd(this_dir)

#------------------------------------------------------------------------------#
#                                Sources and library
#-------------------------------------------------------------------------------
library(tolerance)
library(readxl)
library(MASS)
library(reshape2)

source("R/toolbox_kumaraswamy_GPQ.R")
source("R/toolbox_logit_angular.R")

#------------------------------------------------------------------------------#
#                            Enter Data  Here
#-------------------------------------------------------------------------------
set.seed(123)  # set seed for reproducibility

monomer_ds <-read_excel("data/Example_data.xlsx")

print(monomer_ds, digits = 5)

summary(monomer_ds)
sd(monomer_ds$`Percent monomer`)

# Express response of interest as percentage (0, 1)
resp <- monomer_ds$`Percent monomer`/100

#------------------------------------------------------------------------------#
#                            Logit transformation
#-------------------------------------------------------------------------------

# Mean
mean_logit0 <- t.test(qlogis(resp))$conf.int
mean_logit <- plogis(mean_logit0)

# Tolerance interval
tolint_logit0 <- normtol.int(qlogis(resp),
                             alpha = 0.05, P = 0.95,
                             side = 2)[c(4, 5)]

tolint_logit <- plogis(c(tolint_logit0$`2-sided.lower`, tolint_logit0$`2-sided.upper`))

# Prediction interval
predint_logit0 <- normal_prediction_interval(qlogis(resp),
                                             conf.level = 0.90)
predint_logit <- plogis(predint_logit0)

#tabulate results
out_logit <- rbind(mean_logit, predint_logit, tolint_logit)

out_logit <- as.data.frame(round(100*out_logit, 2))
colnames(out_logit) <- c("95%LCL", "95%UCL")
out_logit$Method <- "Logistic"
out_logit$Parameter <- c("Mean", "PredInt", "TolInt")

out_logit

#-------------------------------------------------------------------------------
#                       arcsine-square-root transformation
#-------------------------------------------------------------------------------
# Mean
mean_asin0 <- t.test(asin_tranform(resp))$conf.int
mean_asin <- asin_back(mean_asin0)

# Tolerance interval
tolint_asin0 <- normtol.int(asin_tranform(resp),
                            alpha = 0.05, P = 0.95,
                            side = 2)[c(4, 5)]

tolint_asin <- asin_back(tolint_asin0)

# Prediction interval
predint_asin0 <- normal_prediction_interval(asin_tranform(resp),
                                            conf.level = 0.90)
predint_asin <- asin_back(predint_asin0)

# Tabulate results
out_asin <- rbind(mean_asin, predint_asin, tolint_asin)

out_asin <- as.data.frame(round(100*out_asin, 2))
colnames(out_asin) <- c("95%LCL", "95%UCL")
out_asin$Method <- "Arcsine"
out_asin$Parameter <- c("Mean", "PredInt", "TolInt")

out_asin

#-------------------------------------------------------------------------------
#             Kumaraswamy Approach 1 : Greenwood statistics
#-------------------------------------------------------------------------------
gpq_green <- gpq_ab(x_data = resp,
                    type = "greenwood", nrep = 1e4,
                    bound = c(0.01, 15))
result_green <- kumar_parameter(gpq_green)

predint_green <- predict_interval_two(gpq_green, clevel = 0.90)

out_green <- rbind(result_green$conf_int_mean[c(1, 5)],
                   predint_green[c(1, 2)], result_green$tolint[c(1, 4)])

out_green <- as.data.frame(round(100*out_green, 2))
colnames(out_green) <- c("95%LCL", "95%UCL")
out_green$Method <- "Approach1"
out_green$Parameter <- c("Mean", "PredInt", "TolInt")

out_green

#-------------------------------------------------------------------------------
#                Kumaraswamy Approach 2 : F-distribution
#-------------------------------------------------------------------------------
gpq_fdist <- gpq_ab(x_data = resp,
                    type = "fdist", nrep = 1e4,
                    bound = c(0.01, 15))

result_fdist <- kumar_parameter(gpq_fdist)

predint_fdist <- predict_interval_two(gpq_fdist, clevel = 0.90)

out_fdist <- rbind(result_fdist$conf_int_mean[c(1, 5)],
                   predint_fdist[c(1, 2)],
                   result_fdist$tolint[c(1, 4)])

out_fdist <- as.data.frame(round(100*out_fdist, 2))
colnames(out_fdist) <- c("95%LCL", "95%UCL")
out_fdist$Method <- "Approach2"
out_fdist$Parameter <- c("Mean", "PredInt", "TolInt")

out_fdist

#-------------------------------------------------------------------------------
#             Kumaraswamy Approach 3 : Wang
#-------------------------------------------------------------------------------
gpq_wang <- gpq_ab(x_data = resp,
                   type = "wang", nrep = 1e4,
                   bound = c(0.01, 2))

result_wang <- kumar_parameter(gpq_wang)

predint_wang <- predict_interval_two(gpq_wang, clevel = 0.90)

out_wang <- rbind(result_wang$conf_int_mean[c(1, 5)], predint_wang[c(1, 2)],
                  result_wang$tolint[c(1, 4)])

out_wang <- as.data.frame(round(100*out_wang, 2))
colnames(out_wang) <- c("95%LCL", "95%UCL")
out_wang$Method <- "Wang"
out_wang$Parameter <- c("Mean", "PredInt", "TolInt")

out_wang

#-------------------------------------------------------------------------------
#                       tabulate data
#-------------------------------------------------------------------------------

combined_df <- rbind(out_asin, out_logit, out_green, out_fdist, out_wang)

long_df <- melt(combined_df, id.vars = c("Method", "Parameter"),
                variable.name = "CL",
                value.name = "Bound")

tab <- xtabs(Bound ~ Parameter + Method + CL, data = long_df)
ftable(tab)

#------------------------------------------------------------------------------#
#                                Closing
#------------------------------------------------------------------------------#
cat("\n"); cat("Time Elapsed (hour): ", (proc.time()[3]-start.time)/3600, "\n")
rm(list=ls(all=TRUE))
#------------------------------------------------------------------------------#
#                             End of script
#------------------------------------------------------------------------------#
