## -----------------------------------------------------------------------------
library(BayesianMCPMod)
library(RBesT)
library(clinDR)
library(dplyr)
library(tibble)
library(reactable)

set.seed(7015)

#' Display Parameters Table
#'
#' This function generates a markdown table displaying the names and values of parameters
#' from a named list.
#'
#' @param named_list A named list where each name represents a parameter name and the list
#'   element represents the parameter value. Date values in the list are automatically
#'   converted to character strings for display purposes.
#'
#' @return Prints a markdown table with two columns: "Parameter Name" and "Parameter Values".
#'   The function does not return a value but displays the table directly to the output.
#'
#' @importFrom knitr kable
#' @examples
#' params <- list("Start Date" = as.Date("2020-01-01"),
#'                "End Date" = as.Date("2020-12-31"),
#'                "Threshold" = 10)
#' display_params_table(params)
#'
#' @export
display_params_table <- function(named_list) {
  display_table <- data.frame()
  value_names <- data.frame()
  for (i in 1:length(named_list)) {
    # dates will display as numeric by default, so convert to char first
    if (class(named_list[[i]]) == "Date") {
      named_list[[i]] = as.character(named_list[[i]])
    }
    if (!is.null(names(named_list[[i]]))) {
      value_names <- rbind(value_names, paste(names(named_list[[i]]), collapse = ', '))
    }
    values <- data.frame(I(list(named_list[[i]])))
    display_table <- rbind(display_table, values)
  }
  
  round_numeric <- function(x, digits = 3) {
    if (is.numeric(x)) {
      return(round(x, digits))
    } else {
      return(x)
    }
  }
  
  display_table[1] <- lapply(display_table[1], function(sublist) {
    lapply(sublist, round_numeric)
  })
  
  class(display_table[[1]]) <- "list"
  
  if (nrow(value_names) == 0) {
    knitr::kable(
      cbind(names(named_list), display_table),
      col.names = c("Name", "Value")
    )
  } else {
    knitr::kable(
      cbind(names(named_list), value_names, display_table),
      col.names = c("Name", "Value Labels", "Value")
    )
  }
}

## ----Historical Data for Control Arm------------------------------------------
data("metaData")
dataset     <- filter(as.data.frame(metaData), bname == "BRINTELLIX")
histcontrol <- filter(
  dataset,
  dose       == 0,
  primtime   == 8,
  indication == "MAJOR DEPRESSIVE DISORDER",
  protid     != 5)

hist_data   <- data.frame(
  trial = histcontrol$nctno,
  est   = histcontrol$rslt,
  se    = histcontrol$se,
  sd    = histcontrol$sd,
  n     = histcontrol$sampsize)

## ----Defining MAP prior function----------------------------------------------
getPriorList <- function (
  
  hist_data,
  dose_levels,
  dose_names    = NULL,
  robust_weight = 0.5
  
) {
  
  sd_tot <- with(hist_data, sum(sd * n) / sum(n))
  
  gmap <- RBesT::gMAP(
    formula    = cbind(est, se) ~ 1 | trial,
    weights    = hist_data$n,
    data       = hist_data,
    family     = gaussian,
    beta.prior = cbind(0, 100 * sd_tot),
    tau.dist   = "HalfNormal",
    tau.prior  = cbind(0, sd_tot / 4))
  
  prior_ctr <- RBesT::automixfit(gmap)
  
  if (!is.null(robust_weight)) {
    
    prior_ctr <- suppressMessages(RBesT::robustify(
      priormix = prior_ctr,
      weight   = robust_weight,
      sigma    = sd_tot))
    
  }
  
  prior_trt <- RBesT::mixnorm(
    comp1 = c(w = 1, m = summary(prior_ctr)[1], n = 1),
    sigma = sd_tot,
    param = "mn")
  
  prior_list <- c(list(prior_ctr),
                  rep(x     = list(prior_trt),
                      times = length(dose_levels[-1])))
  
  if (is.null(dose_names)) {
    
    dose_names <- c("Ctr", paste0("DG_", seq_along(dose_levels[-1])))
    
  }
  
  names(prior_list) <- dose_names
  
  return (prior_list)
  
}

## ----Getting the MAP prior----------------------------------------------------
dose_levels <- c(0, 2.5, 5, 10)

prior_list  <- getPriorList(
  hist_data     = hist_data,
  dose_levels   = dose_levels,
  robust_weight = 0.3)

getESS(prior_list)

## -----------------------------------------------------------------------------
# Guesstimate estimation
exp_guesst  <- DoseFinding::guesst(
  model = "exponential", 
  d = 5, p = 0.2, Maxd = max(dose_levels)
)
emax_guesst <- DoseFinding::guesst(
  model = "emax",
  d = 2.5, p = 0.9
)
sigEmax_guesst <- DoseFinding::guesst(
  model = "sigEmax",
  d = c(2.5, 5), p = c(0.5, 0.95)
)
logistic_guesst <- DoseFinding::guesst(
  model = "logistic",
  d = c(5, 10), p = c(0.1, 0.85)
)

## -----------------------------------------------------------------------------
betaMod_params <- c(delta1 = 1, delta2 = 1)
quadratic_params <- c(delta2 = -0.1)

## -----------------------------------------------------------------------------
mods <- DoseFinding::Mods(
  linear      = NULL,
  # guesstimate scale
  exponential = exp_guesst,
  emax        = emax_guesst,
  sigEmax     = sigEmax_guesst,
  logistic    = logistic_guesst,
  # parameter scale
  betaMod     = betaMod_params,
  quadratic   = quadratic_params,
  # Options for all models
  doses       = dose_levels,
  maxEff      = -1,
  placEff     = -12.8
)

plot(mods)

## -----------------------------------------------------------------------------
display_params_table(mods)

## -----------------------------------------------------------------------------
knitr::kable(DoseFinding::getResp(mods, doses = dose_levels))

## -----------------------------------------------------------------------------
data("metaData")

trial_data <- dplyr::filter(
  dplyr::filter(tibble::tibble(metaData), bname == "BRINTELLIX"),
  primtime == 8,
  indication == "MAJOR DEPRESSIVE DISORDER",
  protid == 5
)

n_patients <- c(128, 124, 129, 122)

## -----------------------------------------------------------------------------
posterior <- getPosterior(
  prior_list = prior_list,
  mu_hat = trial_data$rslt,
  S_hat = trial_data$se,
  calc_ess = TRUE
)

knitr::kable(summary(posterior))

## -----------------------------------------------------------------------------
crit_pval <- getCritProb(
  mods           = mods,
  dose_levels    = dose_levels,
  se_new_trial   = trial_data$se,
  alpha_crit_val = 0.05
)

contr_mat <- getContr(
  mods         = mods,
  dose_levels  = dose_levels,
  sd_posterior = summary(posterior)[, 2]
)

## -----------------------------------------------------------------------------
# # i) the frequentist contrast
# contr_mat_prior <- getContr(
#   mods           = mods,
#   dose_levels    = dose_levels,
#   dose_weights   = n_patients,
#   prior_list     = prior_list)
# # ii) re-estimated frequentist contrasts
# contr_mat_prior <- getContr(
#   mods           = mods,
#   dose_levels    = dose_levels,
#   se_new_trial   = trial_data$se)
# # iii)  Bayesian approach using number of patients for new trial and prior distribution
# contr_mat_prior <- getContr(
#   mods           = mods,
#   dose_levels    = dose_levels,
#   dose_weights   = n_patients,
#   prior_list     = prior_list)

## -----------------------------------------------------------------------------
BMCP_result <- performBayesianMCP(
  posterior_list = posterior,
  contr          = contr_mat, 
  crit_prob_adj  = crit_pval)

## -----------------------------------------------------------------------------
BMCP_result

## -----------------------------------------------------------------------------
# If simple = TRUE, uses approx posterior
# Here we use complete posterior distribution
fit <- getModelFits(
  models      = mods,
  dose_levels = dose_levels,
  posterior   = posterior,
  simple      = FALSE)

## -----------------------------------------------------------------------------
display_params_table(stats::predict(fit, doses = c(0, 2.5, 4, 5, 7, 10)))

## -----------------------------------------------------------------------------
plot(fit)

## -----------------------------------------------------------------------------
plot(fit, cr_bands = TRUE)

## -----------------------------------------------------------------------------
bootstrap_quantiles <- getBootstrapQuantiles(
  model_fits = fit,
  quantiles  = c(0.025, 0.5, 0.975),
  doses      = c(0, 2.5, 4, 5, 7, 10),
  n_samples  = 6
)

## -----------------------------------------------------------------------------
reactable::reactable(
  data = bootstrap_quantiles,
  groupBy = "models",
  columns = list(
    doses = colDef(aggregate = "count", format = list(aggregated = colFormat(suffix = " doses"))),
    "2.5%" = colDef(aggregate = "mean", format = list(aggregated = colFormat(prefix = "mean = ", digits = 2), cell = colFormat(digits = 4))),
    "50%" = colDef(aggregate = "mean", format = list(aggregated = colFormat(prefix = "mean = ", digits = 2), cell = colFormat(digits = 4))),
    "97.5%" = colDef(aggregate = "mean", format = list(aggregated = colFormat(prefix = "mean = ", digits = 2), cell = colFormat(digits = 4)))
  )
)

## -----------------------------------------------------------------------------
# performBayesianMCPMod(
#   posterior_list   = posterior,
#   contr            = contr_mat,
#   crit_prob_adj    = crit_pval,
#   simple           = FALSE)

