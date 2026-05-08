#' Judge data without covariates
#'
#' A simulated dataset for the judge example without covariates.
#' Generated using \code{GenData_nocov()}.
#'
#' @format A data frame with observable variables:
#' \describe{
#'   \item{group}{Group identifier}
#'   \item{X}{Treatment variable}
#'   \item{Y}{Outcome variable}
#'   \item{e}{Error term}
#'   \item{MX}{Residual of X on Z}
#'   \item{Me}{Residual of e on Z}
#'   \item{MY}{Residual of Y on Z}
#' }
"dnc"

#' Interacted data with covariates
#'
#' A simulated dataset for the QOB example with covariates.
#' Generated using \code{GenData_cov()}.
#'
#' @format A data frame with observable variables:
#' \describe{
#'   \item{group}{Group identifier}
#'   \item{groupW}{Covariate group identifier (e.g., State)}
#'   \item{X}{Treatment variable}
#'   \item{Y}{Outcome variable}
#'   \item{e}{Error term}
#'   \item{MX}{Residuals}
#'   \item{Me}{Residuals}
#'   \item{MY}{Residuals}
#' }
"dc"

#' Nonviolent Misdemeanor Prosecution in Suffolk County
#'
#' A dataset containing administrative records from the Suffolk County District
#' Attorney's Office (SCDAO) used to estimate the causal effects of
#' misdemeanor prosecution on recidivism
#' @format A data frame with observable variables
#' @source Agan, A., Doleac, J. L., & Harvey, A. (2023). Misdemeanor Prosecution.
#' The Quarterly Journal of Economics, 138(3), 1453-1505.
"suffolk"
