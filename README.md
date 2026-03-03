[README.md](https://github.com/user-attachments/files/25726012/README.md)
# Fiducial-Based Statistical Intervals for Bounded Bioanalytical Data Using the Kumaraswamy Distribution

**Jorge Quiroz, Jingwei Xiong, Satrajit Roychoudhury, Heliang Shi**

## Overview

This repository provides supplementary R code for the paper:

> Quiroz, J., Xiong, J., Roychoudhury, S., & Shi, H. (2026). Fiducial-Based Statistical Intervals for Bounded Bioanalytical Data Using the Kumaraswamy Distribution. *Pharmaceutical Statistics*.

The code implements fiducial inference methods for the Kumaraswamy distribution applied to continuous data bounded between 0 and 1. It constructs confidence intervals for the mean, prediction intervals, and tolerance intervals using generalized pivotal quantities (GPQ), and compares them with classical transformation-based approaches (logit-normal and arcsine-square-root-normal).

## Methods Implemented

| Method | Description |
|--------|-------------|
| **Approach 1 (Greenwood)** | GPQ for parameter *a* based on the Greenwood statistic of uniform spacings |
| **Approach 2 (F-distribution)** | GPQ for parameter *a* based on an F-distribution criterion from split samples |
| **Approach 3 (Wang)** | GPQ for parameter *a* based on Wang's chi-square statistic |
| **Logit-normal** | Logistic transformation followed by normal-theory intervals |
| **Arcsine-sqrt-normal** | Arcsine-square-root transformation followed by normal-theory intervals |

## Repository Structure

```
code_repo_new/
├── R/
│   ├── toolbox_kumaraswamy_GPQ.R   # GPQ functions for Kumaraswamy parameters
│   └── toolbox_logit_angular.R     # Logit and angular transformation utilities
├── data/
│   └── Example_data.xlsx           # Example dataset (percent monomer)
├── Analysis_of_bounded_data.R      # Main analysis script
├── README.md
└── LICENSE
```

## Requirements

- **R** >= 4.0.0
- Required packages:
  - `tolerance`, `readxl`, `MASS`, `reshape2`
  - `extraDistr`, `rootSolve`, `numbers`, `greybox`

Install missing packages in R:

```r
install.packages(c("tolerance", "readxl", "MASS", "reshape2",
                   "extraDistr", "rootSolve", "numbers", "greybox"))
```

## Usage

### Quick Start

1. Open R and set the working directory to the repository folder:

```r
setwd("path/to/code_repo_new")
```

2. Run the analysis:

```r
source("Analysis_of_bounded_data.R")
```

### Using Your Own Data

To apply the methods to your own bounded data:

1. Source the toolbox files:

```r
source("R/toolbox_kumaraswamy_GPQ.R")
source("R/toolbox_logit_angular.R")
```

2. Prepare your data as a numeric vector in the (0, 1) interval:

```r
resp <- your_data / 100  # if data are percentages
```

3. Compute statistical intervals using the Kumaraswamy GPQ approach:

```r
set.seed(123)

# Approach 1: Greenwood
gpq_result <- gpq_ab(x_data = resp, type = "greenwood", nrep = 1e4,
                     bound = c(0.01, 15))

# Confidence interval for the mean
params <- kumar_parameter(gpq_result)
params$conf_int_mean

# Prediction interval (90%)
predict_interval_two(gpq_result, clevel = 0.90)

# Tolerance interval (95% content, 95% confidence)
params$tolint
```

4. For transformation-based intervals:

```r
# Logit transformation
library(tolerance)
mean_logit <- plogis(t.test(qlogis(resp))$conf.int)
tolint_logit <- normtol.int(qlogis(resp), alpha = 0.05, P = 0.95, side = 2)
predint_logit <- plogis(normal_prediction_interval(qlogis(resp), conf.level = 0.90))

# Arcsine-square-root transformation
mean_asin <- asin_back(t.test(asin_tranform(resp))$conf.int)
tolint_asin <- normtol.int(asin_tranform(resp), alpha = 0.05, P = 0.95, side = 2)
predint_asin <- asin_back(normal_prediction_interval(asin_tranform(resp), conf.level = 0.90))
```

## Description of Functions

### `toolbox_kumaraswamy_GPQ.R`

| Function | Description |
|----------|-------------|
| `gpq_ab()` | Generate GPQ samples for Kumaraswamy parameters *a* and *b* |
| `kumar_parameter()` | Compute confidence intervals for *b*, the mean, and tolerance intervals from GPQ samples |
| `kumar_mean()` | Mean of the Kumaraswamy distribution |
| `kumar_quantile()` | Quantile function of the Kumaraswamy distribution |
| `predict_interval()` | One-sided fiducial prediction interval |
| `predict_interval_two()` | Two-sided fiducial prediction interval |
| `pivot_a()` | Exact pivotal quantity for parameter *a* (Greenwood, F-distribution, or Wang) |
| `rGreenwood()` | Generate random Greenwood statistic |
| `greenwood.dist()` | Generate random normalized Greenwood statistic |

### `toolbox_logit_angular.R`

| Function | Description |
|----------|-------------|
| `normal_prediction_interval()` | Prediction interval based on normal theory |
| `asin_tranform()` | Arcsine-square-root transformation for data in (0, 1) |
| `asin_back()` | Back-transformation from arcsine-square-root scale |

## Citation

```bibtex
@article{quiroz2026fiducial,
  title   = {Fiducial-Based Statistical Intervals for Bounded Bioanalytical
             Data Using the {K}umaraswamy Distribution},
  author  = {Quiroz, Jorge and Xiong, Jingwei and Roychoudhury, Satrajit
             and Shi, Heliang},
  journal = {Pharmaceutical Statistics},
  year    = {2026},
  note    = {In review}
}
```

## License

MIT License. See [LICENSE](LICENSE) for details.
