#' @description
#' This function takes a regression estimated using fixest with the sunab option
#' and extracts the aggregated event-study coefficients and their variance-covariance matrix
#' @param sunab_fixest The result of a fixest call using the sunab option
#' @returns A list containing beta (the event-study coefficients),
#'          sigma (the variance-covariance matrix), and
#'          cohorts (the relative times corresponding to beta, sigma)

sunab_beta_vcv <-
function(sunab_fixest){

  ## The following code block extracts the weights on individual coefs used in
  # the fixest aggregation ##
  sunab_agg   <- sunab_fixest$model_matrix_info$sunab$agg_period
  sunab_names <- names(sunab_fixest$coefficients)
  sunab_sel   <- grepl(sunab_agg, sunab_names, perl=TRUE)
  sunab_names <- sunab_names[sunab_sel]
  if(!is.null(sunab_fixest$weights)){
    sunab_wgt <- colSums(sunab_fixest$weights * sign(model.matrix(sunab_fixest)[, sunab_names, drop=FALSE]))
  } else {
    sunab_wgt <- colSums(sign(model.matrix(sunab_fixest)[, sunab_names, drop=FALSE]))
  }

  #Construct matrix sunab_trans such that sunab_trans %*% non-aggregated coefs = aggregated coefs,
  sunab_cohorts <- as.numeric(gsub(paste0(".*", sunab_agg, ".*"), "\\2", sunab_names, perl=TRUE))
  sunab_mat     <- model.matrix(~ 0 + factor(sunab_cohorts))
  sunab_trans   <- solve(t(sunab_mat) %*% (sunab_wgt * sunab_mat)) %*% t(sunab_wgt * sunab_mat)

  #Get the coefs and vcv
  sunab_coefs   <- sunab_trans %*% cbind(sunab_fixest$coefficients[sunab_sel])
  sunab_vcov    <- sunab_trans %*% sunab_fixest$cov.scaled[sunab_sel, sunab_sel] %*% t(sunab_trans)

  return(list(beta    = sunab_coefs,
              sigma   = sunab_vcov,
              cohorts = sort(unique(sunab_cohorts))))
}
