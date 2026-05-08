#' Leave-Three-Out Variance Estimator for LM Statistic
#'
#' @description
#' Calculates the consistent variance estimator (\eqn{\hat{V}_{LM}}) for the Lagrange Multiplier (LM)
#' test statistic using the "Leave-Three-Out" (L3O) adjustment.
#'
#' @param X Numeric vector of length n. The endogenous variable.
#' @param e Numeric vector of length n. Residuals under the null hypothesis
#'   (\eqn{e = Y - X\beta_0}).
#' @param P Matrix of dimension n x n. Projection matrix of the instruments (and potentially covariates).
#'   Corresponds to matrix \eqn{P} or \eqn{H_Q} in the paper.
#' @param G Matrix of dimension n x n. Weighting matrix used in the JIVE/UJIVE estimator.
#'   For standard JIVE, G is equal to P. For UJIVE with covariates, G is the adjusted
#'   matrix defined in Section 3.1.
#' @param noisy Logical. If \code{TRUE}, print progress dots during the loop.
#'   Defaults to \code{FALSE}.
#'
#' @details
#' This variance estimator is robust to both many weak instruments and heterogeneous treatment effects.
#' It corrects for biases in variance estimation that arise when reduced-form coefficients
#' are not consistently estimable.
#'
#' The function computes the variance estimator defined in Equation (9) of Yap (2025):
#' \deqn{\hat{V}_{LM} = A_1 + A_2 + A_3 + A_4 + A_5}
#'
#' It iterates through each observation \eqn{i} to compute the necessary adjustments
#' (Leave-Three-Out determinants \eqn{D_{ijk}}) and aggregates the components using
#' optimized matrix operations to handle the double sums over \eqn{j} and \eqn{k}.
#'
#' Specifically:
#' \itemize{
#'   \item \strong{A1, A2, A3} capture the core variance components involving interactions between
#'   the instruments, endogenous variable, and residuals.
#'   \item \strong{A4, A5} are correction terms that account for the variability from estimating
#'   the reduced-form coefficients (which cannot be treated as fixed in the many-instrument setting).
#' }
#'
#' The calculation relies on determinants \eqn{D_{ij}} and \eqn{D_{ijk}} derived from the annihilator
#' matrix \eqn{M = I - P} to ensure the estimator is unbiased.
#'
#' @return Scalar. The estimated variance \eqn{\hat{V}_{LM}}.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity".
#'
#' @export
L3Ovar_iloop_cov <- function(X, e, P, G, noisy = FALSE) {
  n <- length(X)
  M <- diag(n) - P

  # Stuff that can be calculated once
  dM <- matrix(diag(M), ncol = 1) # force column vector
  D2 <- dM %*% t(dM) - M * M
  onesN <- matrix(rep(1, n), ncol = 1)
  recD2 <- 1 / D2
  diag(recD2) <- 0
  Me <- M %*% e
  MX <- M %*% X
  Poff <- P - diag(diag(P))
  Goff <- G - diag(diag(G))

  A1vec <- A2vec <- A3vec <- A4vec <- A5vec <- rep(0, n)
  for (i in 1:n) {
    # Calculation conditioned on i
    Mi <- matrix(M[, i], ncol = 1) # force column vector
    Pi <- matrix(P[, i], ncol = 1) # force column vector
    Gicol <- matrix(G[, i], ncol = 1) # force column vector
    Girow <- matrix(G[i, ], ncol = 1) # force column vector
    D3i <- M[i, i] * D2 - (dM %*% t(Mi)^2 + Mi^2 %*% t(dM) - 2 * M * (Mi %*% t(Mi)))
    Di <- matrix(D2[, i], ncol = 1)
    D2D3i <- D2 / D3i
    D2D3i[i, ] <- 0
    D2D3i[, i] <- 0
    diag(D2D3i) <- 0
    recD3i <- 1 / D3i
    recD3i[i, ] <- 0
    recD3i[, i] <- 0
    diag(recD3i) <- 0
    recD2i <- matrix(1 / D2[, i], ncol = 1)
    recD2i[i] <- 0
    Poffi <- Poff[, i]
    Gofficol <- Goff[, i]
    Goffirow <- matrix(Goff[i, ], ncol = 1)

    A11i <- t(X * Girow) %*% D2D3i %*% (Girow * X) * (t(Mi) %*% e)
    A12i <- t((Me) * X * Girow * Mi) %*% recD3i %*% (Girow * X * dM) -
      t(X * Girow * Mi) %*% (recD3i * M) %*% (Girow * X * (Me))
    A13i <- t(dM * X * Girow) %*% ((onesN %x% t(Me)) * recD3i) %*% (Girow * Mi * X)
    A14i <- t((Me) * X * Girow) %*% (M * recD3i) %*% (Gicol * Mi * X)
    A15i <- (t(Mi) %*% e) * (t(Goffirow^2 * recD2i) %*% (dM * X^2)) -
      (t(Goffirow^2 * Mi * recD2i) %*% (Me * X^2))

    A21i <- t(X * Girow) %*% D2D3i %*% (Gicol * e) * (t(Mi) %*% X)
    A22i <- t((MX) * X * Girow * Mi) %*% recD3i %*% (Gicol * e * dM) -
      t(X * Girow * Mi) %*% (recD3i * M) %*% (Gicol * e * (MX))
    A23i <- t(dM * X * Girow) %*% ((onesN %x% t(MX)) * recD3i) %*% (Gicol * Mi * e)
    A24i <- t((MX) * X * Girow) %*% (M * recD3i) %*% (Gicol * Mi * e)
    A25i <- (t(Mi) %*% X) * (t(Goffirow * Gofficol * recD2i) %*% (dM * X * e)) -
      (t(Goffirow * Gofficol * Mi * recD2i) %*% (MX * X * e))

    A31i <- t(e * Gicol) %*% D2D3i %*% (Gicol * e) * (t(Mi) %*% X)
    A32i <- t((MX) * e * Gicol * Mi) %*% recD3i %*% (Gicol * e * dM) -
      t(e * Gicol * Mi) %*% (recD3i * M) %*% (Gicol * e * (MX))
    A33i <- t(dM * e * Gicol) %*% ((onesN %x% t(MX)) * recD3i) %*% (Gicol * Mi * e)
    A34i <- t((MX) * e * Gicol) %*% (M * recD3i) %*% (Gicol * Mi * e)
    A35i <- (t(Mi) %*% X) * (t(Gofficol^2 * recD2i) %*% (dM * e * e)) -
      (t(Gofficol^2 * Mi * recD2i) %*% (MX * e * e))

    A41i <- t(Gicol^2 * e * dM * recD2i * Me) %*% (D2D3i) %*% (Mi * X) -
      t(Gicol^2 * e * recD2i * Me) %*% (D2D3i * M * (onesN %x% t(X))) %*% (Mi)
    A42i <- t(Gicol^2 * e * dM * recD2i) %*% (recD3i * M) %*% (Mi^2 * X) -
      t(Gicol^2 * e * Mi * recD2i) %*% (recD3i * M * M) %*% (Mi * X) -
      t(Gicol^2 * e * dM * Mi * recD2i) %*% ((onesN %x% t(dM)) * recD3i) %*% (Mi * X) +
      t(Gicol^2 * e * Mi^2 * recD2i) %*% (recD3i * M) %*% (dM * X)
    A43i <- t(Gicol^2 * e * dM * Mi * recD2i) %*% (recD3i) %*% (Mi^2 * X * Me) -
      t(Gicol^2 * e * Mi^2 * recD2i) %*% (recD3i * M) %*% (Mi * X * Me) -
      M[i, i] * t(Gicol^2 * e * dM * recD2i) %*% (recD3i * M) %*% (Mi * X * Me) +
      M[i, i] * t(Gicol^2 * e * Mi * recD2i) %*% (recD3i * M * M) %*% (Me * X)
    A44i <- X[i] * M[i, i] * (t(Gofficol^2 * recD2i) %*% (Me * e)) -
      X[i] * (t(Mi) %*% e) * (t(Gofficol^2 * recD2i) %*% (Mi * e))

    A51i <- t(Girow * Gicol * e * dM * recD2i * MX) %*% (D2D3i) %*% (Mi * X) -
      t(Girow * Gicol * e * recD2i * MX) %*% (D2D3i * M * (onesN %x% t(X))) %*% (Mi)
    A52i <- t(Girow * Gicol * e * dM * recD2i) %*% (recD3i * M) %*% (Mi^2 * X) -
      t(Girow * Gicol * e * Mi * recD2i) %*% (recD3i * M * M) %*% (Mi * X) -
      t(Girow * Gicol * e * dM * Mi * recD2i) %*% ((onesN %x% t(dM)) * recD3i) %*% (Mi * X) +
      t(Girow * Gicol * e * Mi^2 * recD2i) %*% (recD3i * M) %*% (dM * X)
    A53i <- t(Girow * Gicol * e * dM * Mi * recD2i) %*% (recD3i) %*% (Mi^2 * X * MX) -
      t(Girow * Gicol * e * Mi^2 * recD2i) %*% (recD3i * M) %*% (Mi * X * MX) -
      M[i, i] * t(Girow * Gicol * e * dM * recD2i) %*% (recD3i * M) %*% (Mi * X * MX) +
      M[i, i] * t(Girow * Gicol * e * Mi * recD2i) %*% (recD3i * M * M) %*% (MX * X)
    A54i <- X[i] * M[i, i] * (t(Goffirow * Gofficol * recD2i) %*% (MX * e)) -
      X[i] * (t(Mi) %*% X) * (t(Goffirow * Gofficol * recD2i) %*% (Mi * e))

    A1vec[i] <- A11i - A12i - A13i + A14i + A15i
    A2vec[i] <- A21i - A22i - A23i + A24i + A25i
    A3vec[i] <- A31i - A32i - A33i + A34i + A35i
    A4vec[i] <- A41i + A42i * (M[i, ] %*% e) + A43i + A44i
    A5vec[i] <- A51i + A52i * (M[i, ] %*% X) + A53i + A54i

    if (noisy) {
      if (i %% 10 == 0) cat(i / n, " ")
    }
  }

  sum(A1vec * e + 2 * A2vec * e + A3vec * X - A4vec * X - A5vec * e)
}

#' Compute L3O Variance for Score Statistic (Grouped Design)
#'
#' @description
#' Computes the Leave-Three-Out (L3O) variance estimator (\eqn{\hat{V}_{LM}}) for the score statistic
#' in a grouped instrument design with covariates. This function simultaneously calculates
#' all five variance components by exploiting the block-diagonal symmetry of the projection matrices
#' inherent to discrete instruments.
#'
#' @param df Data frame. Contains the variables used in estimation.
#' @param group Column name (unquoted). The instrument grouping variable.
#' @param groupW Column name (unquoted). The covariate stratification variable.
#' @param X Column name (unquoted). The endogenous regressor.
#' @param e Column name (unquoted). Residuals under the null hypothesis (\eqn{Y - X\beta_0}).
#' @param MX Column name (unquoted). Leverage-adjusted regressor (\eqn{M_i X_i}).
#' @param Me Column name (unquoted). Leverage-adjusted residual (\eqn{M_i e_i}).
#'
#' @details
#' This function implements the variance estimator \eqn{\hat{V}_{LM}} for the Limited Information
#' Maximum Likelihood (LIML) or UJIVE estimator. It computes:
#' \deqn{\hat{V} = A_1 + 2A_2 + A_3 - A_4 - A_5}
#'
#' Because the design implies symmetric weighting matrices (\eqn{G_{ij} = G_{ji}}), the function
#' simplifies the asymmetric components.
#'
#' \strong{Components:}
#' \itemize{
#'   \item \strong{A1, A2, A3}: Variance and covariance of the score statistic. Calculated via 5 sub-components each.
#'   \item \strong{A4, A5}: Bias correction terms for "own-observation" variance contributions. Calculated via 4 sub-components each.
#' }
#'
#' @return Scalar. The estimated variance \eqn{\hat{V}_{LM}}.
#'
#' @export
L3Ovar_gloop_cov <- function(df, group, groupW, X, e, MX, Me) {
  df$group <- eval(substitute(group), df) # Find column named group in dataframe df, then create new column in the dataframe with the same name
  df$groupW <- eval(substitute(groupW), df)
  df$X <- eval(substitute(X), df)
  df$e <- eval(substitute(e), df)
  df$MX <- eval(substitute(MX), df)
  df$Me <- eval(substitute(Me), df)

  A11vecs <- A12vecs <- A13vecs <- A14vecs <- A15vecs <- rep(0, max(unique(df$group))) # Create empty/zero vectors
  A21vecs <- A22vecs <- A23vecs <- A24vecs <- A25vecs <- rep(0, max(unique(df$group)))
  A31vecs <- A32vecs <- A33vecs <- A34vecs <- A35vecs <- rep(0, max(unique(df$group)))
  A41vecs <- A42vecs <- A43vecs <- A44vecs <- rep(0, max(unique(df$group)))
  A51vecs <- A52vecs <- A53vecs <- A54vecs <- rep(0, max(unique(df$group)))
  for (s in unique(df$groupW)) { # Start to loop through every groupW
    ds <- df[df$groupW == s, ] # Temporary dataset for each s
    ZQ <- matrix(0, nrow = length(ds$group), ncol = length(unique(ds$group))) # Create empty matrix, col = number instrument groups, row = number of obs
    ds$groupidx <- Getgroupindex(ds, group) # Numbering group
    ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1 # Add 1 values in matrix ZQ where i=j
    ZW <- matrix(1, nrow = length(ds$groupW), ncol = length(unique(ds$groupW))) # Create matrix of ones, col = 1 (unique in ds$groupW), row=numberofgroupW

    PQ <- ZQ %*% solve(t(ZQ) %*% ZQ) %*% t(ZQ) # Projection matrix P = X (X'X)^-1 X'
    PW <- ZW %*% solve(t(ZW) %*% ZW) %*% t(ZW)

    # calculate values specific to this subset
    Gs <- solve(diag(nrow(ds)) - diag(diag(PQ))) %*% (PQ - diag(diag(PQ))) - # L1O adjustment: Remove diagonal, rescale weights, subtract groupW variation
      solve(diag(nrow(ds)) - diag(diag(PW))) %*% (PW - diag(diag(PW)))
    Ps <- PQ
    Ms <- diag(nrow(ds)) - Ps # Residuals M = I - P
    dMs <- matrix(diag(Ms), ncol = 1) # Extract diagonal of M
    D2s <- dMs %*% t(dMs) - Ms * Ms # D_ij
    recD2s <- 1 / D2s # Invert for weights, then remove diagonal
    diag(recD2s) <- 0

    for (g in unique(ds$group)) { # Start loop of single instrument group
      repidx <- min(which(ds$group == g)) # representative index
      Pis <- matrix(Ps[, repidx], ncol = 1) # Extract projection weights
      Pgs <- ifelse(Pis == 0, 0, 1) %*% matrix(1, ncol = length(Pis), nrow = 1) # identify which observations are linked to this group
      Gis <- matrix(Gs[, repidx], ncol = 1) # L1O weights
      Gis[repidx, 1] <- Gis[repidx + 1, 1] # Put the G back
      Mis <- matrix(-Ps[, repidx], ncol = 1) # Annihilator
      recD2is <- matrix(recD2s[, repidx], ncol = 1) # variance correction
      recD2is[repidx, 1] <- recD2is[repidx + 1, 1] # Put the G back
      D3is <- Ms[repidx, repidx] * D2s - (dMs %*% t(Mis)^2 + Mis^2 %*% t(dMs) - 2 * Ms * (Mis %*% t(Mis)))
      D2D3is <- D2s / D3is
      diag(D2D3is) <- 0
      recD3is <- 1 / D3is
      diag(recD3is) <- 0
      ones <- matrix(rep(1, nrow(ds)), ncol = 1)
      Mes <- matrix(ds$Me, ncol = 1)
      MXs <- matrix(ds$MX, ncol = 1)
      es <- matrix(ds$e, ncol = 1)

      A11vecs[g] <- t(ds$X * Gis) %*% (D2D3is * ivectomats(ds, ds$e * ds$Me, g)) %*% (ds$X * Gis)
      A12vecs[g] <- t(ds$X * ds$Me * Gis * Mis) %*% (recD3is * ivectomats(ds, ds$e, g)) %*% (Gis * ds$X * dMs) -
        t(ds$X * Gis * Mis) %*% (recD3is * Ms * ivectomats(ds, ds$e, g)) %*% (Gis * ds$X * ds$Me)
      A13vecs[g] <- t(ds$X * dMs * Gis) %*% (recD3is * (ones %*% t(Mes)) * ivectomats(ds, ds$e, g)) %*% (ds$X * Gis * Mis)
      A14vecs[g] <- t(ds$X * ds$Me * Gis) %*% (recD3is * Ms * ivectomats(ds, ds$e, g)) %*% (ds$X * Gis * Mis)
      A15vecs[g] <- t(ds$Me * ds$e) %*% (recD2s * Gs^2 * Pgs) %*% (ds$X^2 * dMs) -
        t(ds$e) %*% (Gs^2 * Ms * recD2s * Pgs) %*% (ds$Me * ds$X^2)

      A21vecs[g] <- t(ds$X * Gis) %*% (D2D3is * ivectomats(ds, ds$e * ds$MX, g)) %*% (ds$e * Gis)
      A22vecs[g] <- t(ds$X * ds$MX * Gis * Mis) %*% (recD3is * ivectomats(ds, ds$e, g)) %*% (Gis * ds$e * dMs) -
        t(ds$X * Gis * Mis) %*% (recD3is * Ms * ivectomats(ds, ds$e, g)) %*% (Gis * ds$e * ds$MX)
      A23vecs[g] <- t(ds$X * dMs * Gis) %*% (recD3is * (ones %*% t(MXs)) * ivectomats(ds, ds$e, g)) %*% (ds$e * Gis * Mis)
      A24vecs[g] <- t(ds$X * ds$MX * Gis) %*% (recD3is * Ms * ivectomats(ds, ds$e, g)) %*% (ds$e * Gis * Mis)
      A25vecs[g] <- t(ds$MX * ds$e) %*% (recD2s * Gs^2 * Pgs) %*% (ds$X * ds$e * dMs) -
        t(ds$e) %*% (Gs^2 * Ms * recD2s * Pgs) %*% (ds$MX * ds$X * ds$e)

      A31vecs[g] <- t(ds$e * Gis) %*% (D2D3is * ivectomats(ds, ds$X * ds$MX, g)) %*% (ds$e * Gis)
      A32vecs[g] <- t(ds$e * ds$MX * Gis * Mis) %*% (recD3is * ivectomats(ds, ds$X, g)) %*% (Gis * ds$e * dMs) -
        t(ds$e * Gis * Mis) %*% (recD3is * Ms * ivectomats(ds, ds$X, g)) %*% (Gis * ds$e * ds$MX)
      A33vecs[g] <- t(ds$e * dMs * Gis) %*% (recD3is * (ones %*% t(MXs)) * ivectomats(ds, ds$X, g)) %*% (ds$e * Gis * Mis)
      A34vecs[g] <- t(ds$e * ds$MX * Gis) %*% (recD3is * Ms * ivectomats(ds, ds$X, g)) %*% (ds$e * Gis * Mis)
      A35vecs[g] <- t(ds$MX * ds$X) %*% (recD2s * Gs^2 * Pgs) %*% (ds$e^2 * dMs) -
        t(ds$X) %*% (Gs^2 * Ms * recD2s * Pgs) %*% (ds$MX * ds$e^2)

      A41vecs[g] <- t(Gis^2 * ds$e * dMs * ds$Me * recD2is) %*% (D2D3is * ivectomats(ds, ds$X, g)) %*% (Mis * ds$X) -
        t(Gis^2 * ds$e * ds$Me * recD2is) %*% (D2D3is * Ms * (ones %x% t(ds$X)) * ivectomats(ds, ds$X, g)) %*% (Mis)
      A42vecs[g] <- t(Gis^2 * ds$e * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$X * ds$Me, g)) %*% (Mis^2 * ds$X) -
        t(Gis^2 * ds$e * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$X * ds$Me, g)) %*% (Mis * ds$X) -
        t(Gis^2 * ds$e * dMs * Mis * recD2is) %*% (recD3is * (ones %*% t(dMs)) * ivectomats(ds, ds$X * ds$Me, g)) %*% (Mis * ds$X) +
        t(Gis^2 * ds$e * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$X * ds$Me, g)) %*% (dMs * ds$X)
      A43vecs[g] <- t(Gis^2 * ds$e * dMs * Mis * recD2is) %*% (recD3is * ivectomats(ds, ds$X, g)) %*% (Mis^2 * ds$X * ds$Me) -
        t(Gis^2 * ds$e * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$X, g)) %*% (Mis * ds$X * ds$Me) -
        t(Gis^2 * ds$e * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$X * dMs, g)) %*% (Mis * ds$X * ds$Me) +
        t(Gis^2 * ds$e * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$X * dMs, g)) %*% (ds$X * ds$Me)
      A44vecs[g] <- t(ds$X^2 * dMs) %*% (recD2s * Gs^2 * Pgs) %*% (ds$Me * ds$e) -
        t(ds$X^2 * ds$Me) %*% (recD2s * Gs^2 * Pgs) %*% (Mis * ds$e)

      A51vecs[g] <- t(Gis^2 * ds$e * dMs * ds$MX * recD2is) %*% (D2D3is * ivectomats(ds, ds$e, g)) %*% (Mis * ds$X) -
        t(Gis^2 * ds$e * ds$MX * recD2is) %*% (D2D3is * Ms * (ones %x% t(ds$X)) * ivectomats(ds, ds$e, g)) %*% (Mis)
      A52vecs[g] <- t(Gis^2 * ds$e * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$e * ds$MX, g)) %*% (Mis^2 * ds$X) -
        t(Gis^2 * ds$e * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$e * ds$MX, g)) %*% (Mis * ds$X) -
        t(Gis^2 * ds$e * dMs * Mis * recD2is) %*% (recD3is * (ones %*% t(dMs)) * ivectomats(ds, ds$e * ds$MX, g)) %*% (Mis * ds$X) +
        t(Gis^2 * ds$e * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$e * ds$MX, g)) %*% (dMs * ds$X)
      A53vecs[g] <- t(Gis^2 * ds$e * dMs * Mis * recD2is) %*% (recD3is * ivectomats(ds, ds$e, g)) %*% (Mis^2 * ds$X * ds$MX) -
        t(Gis^2 * ds$e * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$e, g)) %*% (Mis * ds$X * ds$MX) -
        t(Gis^2 * ds$e * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$e * dMs, g)) %*% (Mis * ds$X * ds$MX) +
        t(Gis^2 * ds$e * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$e * dMs, g)) %*% (ds$X * ds$MX)
      A54vecs[g] <- t(ds$X * ds$e * dMs) %*% (recD2s * Gs^2 * Pgs) %*% (ds$MX * ds$e) -
        t(ds$X * ds$e * ds$MX) %*% (recD2s * Gs^2 * Pgs) %*% (Mis * ds$e)

    }
  }

  ret <- (A11vecs - A12vecs - A13vecs + A14vecs + A15vecs) +
    2 * (A21vecs - A22vecs - A23vecs + A24vecs + A25vecs) +
    (A31vecs - A32vecs - A33vecs + A34vecs + A35vecs) -
    (A41vecs + A42vecs + A43vecs + A44vecs) -
    (A51vecs + A52vecs + A53vecs + A54vecs)

  sum(ret)
}

#' Compute L3O Variance for Score Statistic (Grouped Design, No Covariates)
#'
#' @description
#' Computes the Leave-Three-Out (L3O) variance estimator (\eqn{\hat{V}_{LM}}) for the score statistic
#' in a grouped instrument design without covariates.
#' This function leverages the block-diagonal structure of the projection matrix to
#' compute all variance components locally within each group.
#'
#' @param df Data frame. Contains the variables used in estimation.
#' @param group Column name (unquoted). The instrument grouping variable.
#' @param X Column name (unquoted). The endogenous regressor.
#' @param e Column name (unquoted). Residuals under the null hypothesis (\eqn{Y - X\beta_0}).
#' @param MX Column name (unquoted). Leverage-adjusted regressor (\eqn{M_i X_i}).
#' @param Me Column name (unquoted). Leverage-adjusted residual (\eqn{M_i e_i}).
#'
#' @details
#' This function implements the variance estimator \eqn{\hat{V}_{LM}} for the specific case where
#' \eqn{G = P} (symmetric weights) and the design matrix is block-diagonal (no covariates).
#'
#' Instead of operating on \eqn{N \times N} matrices, the function iterates through groups.
#' Within each group subset, the projection matrix is dense but small. The function computes:
#' \deqn{\hat{V} = \sum_g (A_{1g} + 2A_{2g} + A_{3g} - A_{4g} - A_{5g})}
#'
#' \strong{Components:}
#' \itemize{
#'   \item \strong{A1, A2, A3}: Signal, covariance, and error variance terms calculated using
#'   the symmetric weighting matrix \eqn{G=P}.
#'   \item \strong{A4, A5}: Bias correction terms utilizing squared weights \eqn{G_{ij}^2}.
#' }
#'
#' @return Scalar. The estimated variance \eqn{\hat{V}_{LM}}.
#'
#' @export
L3Ovar_gloop_nocov <- function(df, group, X, e, MX, Me) {
  df$group <- eval(substitute(group), df)
  A11vecs <- A12vecs <- A13vecs <- A14vecs <- A15vecs <- rep(0, length(unique(df$group)))
  A21vecs <- A22vecs <- A23vecs <- A24vecs <- A25vecs <- rep(0, length(unique(df$group)))
  A31vecs <- A32vecs <- A33vecs <- A34vecs <- A35vecs <- rep(0, length(unique(df$group)))
  A41vecs <- A42vecs <- A43vecs <- A44vecs <- rep(0, length(unique(df$group)))
  A51vecs <- A52vecs <- A53vecs <- A54vecs <- rep(0, length(unique(df$group)))
  A11vecs <- A12vecs <- A13vecs <- A14vecs <- A15vecs <- rep(0, length(unique(df$group)))
  A21vecs <- A22vecs <- A23vecs <- A24vecs <- A25vecs <- rep(0, length(unique(df$group)))
  A31vecs <- A32vecs <- A33vecs <- A34vecs <- A35vecs <- rep(0, length(unique(df$group)))
  A41vecs <- A42vecs <- A43vecs <- A44vecs <- rep(0, length(unique(df$group)))
  A51vecs <- A52vecs <- A53vecs <- A54vecs <- rep(0, length(unique(df$group)))
  for (g in 1:length(unique(df$group))) {
    ds <- df[df$group == g, ]
    ZQ <- matrix(0, nrow = length(ds$group), ncol = length(unique(ds$group)))
    ds$groupidx <- Getgroupindex(ds, group)
    ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1

    PQ <- ZQ %*% solve(t(ZQ) %*% ZQ) %*% t(ZQ)

    # calculate values specific to this subset
    Gs <- PQ
    Ps <- PQ
    Psoff <- Ps - diag(nrow(ds)) * diag(Ps)
    Ms <- diag(nrow(ds)) - Ps
    dMs <- matrix(diag(Ms), ncol = 1)
    D2s <- dMs %*% t(dMs) - Ms * Ms
    recD2s <- 1 / D2s
    diag(recD2s) <- 0
    # diag0 <- matrix(1,nrow=nrow(ds),ncol=nrow(ds)); diag(diag0) <- 0

    repidx <- min(which(ds$group == g)) # representative index
    Gis <- matrix(Gs[, repidx], ncol = 1)
    Gis[repidx, 1] <- Gis[repidx + 1, 1] # Put the G back
    Mis <- matrix(-Ps[, repidx], ncol = 1)
    recD2is <- matrix(recD2s[, repidx], ncol = 1)
    recD2is[repidx, 1] <- recD2is[repidx + 1, 1] # Put the G back
    D3is <- Ms[repidx, repidx] * D2s - (dMs %*% t(Mis)^2 + Mis^2 %*% t(dMs) - 2 * Ms * (Mis %*% t(Mis)))
    D2D3is <- D2s / D3is
    diag(D2D3is) <- 0
    recD3is <- 1 / D3is
    diag(recD3is) <- 0
    ones <- matrix(rep(1, nrow(ds)), ncol = 1)
    Mes <- matrix(ds$Me, ncol = 1)
    MXs <- matrix(ds$MX, ncol = 1)
    es <- matrix(ds$e, ncol = 1)

    A11vecs[g] <- t(ds$X * Gis) %*% (D2D3is * ivectomats(ds, ds$e * ds$Me, g)) %*% (ds$X * Gis)
    A12vecs[g] <- t(ds$X * ds$Me * Gis) %*% (Ms * Psoff * recD3is * ivectomats(ds, ds$e, g)) %*% (ds$X * dMs) -
      t(ds$X * Gis) %*% (recD3is * Ms * Ms * Psoff * ivectomats(ds, ds$e, g)) %*% (ds$X * ds$Me)
    A13vecs[g] <- t(ds$X * dMs * Gis) %*% (recD3is * (ones %*% t(Mes)) * ivectomats(ds, ds$e, g)) %*% (ds$X * Gis * Mis)
    A14vecs[g] <- t(ds$X * ds$Me * Gis) %*% (recD3is * Ms * Ms * ivectomats(ds, ds$e, g)) %*% (ds$X * Gis)
    A15vecs[g] <- t(ds$Me * ds$e) %*% (recD2s * Psoff^2) %*% (ds$X^2 * dMs) -
      t(ds$e) %*% (Ms * recD2s * Psoff^2) %*% (ds$Me * ds$X^2)

    A21vecs[g] <- t(ds$X * Gis) %*% (D2D3is * ivectomats(ds, ds$e * ds$MX, g)) %*% (ds$e * Gis)
    A22vecs[g] <- t(ds$X * ds$MX * Gis) %*% (Ms * Psoff * recD3is * ivectomats(ds, ds$e, g)) %*% (ds$e * dMs) -
      t(ds$X * Gis) %*% (recD3is * Ms * Ms * Psoff * ivectomats(ds, ds$e, g)) %*% (ds$e * ds$MX)
    A23vecs[g] <- t(ds$X * dMs * Gis) %*% (recD3is * (ones %*% t(MXs)) * ivectomats(ds, ds$e, g)) %*% (ds$e * Gis * Mis)
    A24vecs[g] <- t(ds$X * ds$MX * Gis) %*% (recD3is * Ms * Ms * ivectomats(ds, ds$e, g)) %*% (ds$e * Gis)
    A25vecs[g] <- t(ds$MX * ds$e) %*% (recD2s * Psoff^2) %*% (ds$X * ds$e * dMs) -
      t(ds$e) %*% (Ms * recD2s * Psoff^2) %*% (ds$MX * ds$X * ds$e)

    A31vecs[g] <- t(ds$e * Gis) %*% (D2D3is * ivectomats(ds, ds$X * ds$MX, g)) %*% (ds$e * Gis)
    A32vecs[g] <- t(ds$e * ds$MX * Gis) %*% (Ms * Psoff * recD3is * ivectomats(ds, ds$X, g)) %*% (ds$e * dMs) -
      t(ds$e * Gis) %*% (recD3is * Ms * Ms * Psoff * ivectomats(ds, ds$X, g)) %*% (ds$e * ds$MX)
    A33vecs[g] <- t(ds$e * dMs * Gis) %*% (recD3is * (ones %*% t(MXs)) * ivectomats(ds, ds$X, g)) %*% (ds$e * Gis * Mis)
    A34vecs[g] <- t(ds$e * ds$MX * Gis) %*% (recD3is * Ms * Ms * ivectomats(ds, ds$X, g)) %*% (ds$e * Gis)
    A35vecs[g] <- t(ds$MX * ds$X) %*% (recD2s * Psoff^2) %*% (ds$e^2 * dMs) -
      t(ds$X) %*% (Ms * recD2s * Psoff^2) %*% (ds$MX * ds$e^2)

    A41vecs[g] <- t(Gis^2 * ds$e * dMs * ds$Me) %*% (recD3is * ivectomats(ds, ds$X, g)) %*% (Mis * ds$X) -
      t(Gis^2 * ds$e * ds$Me) %*% (recD3is * Ms * (ones %x% t(ds$X)) * ivectomats(ds, ds$X, g)) %*% (Mis)
    A42vecs[g] <- t(Gis^2 * ds$e * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$X * ds$Me, g)) %*% (Mis^2 * ds$X) -
      t(Gis^2 * ds$e * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$X * ds$Me, g)) %*% (Mis * ds$X) -
      t(Gis^2 * ds$e * dMs * Mis * recD2is) %*% (recD3is * (ones %*% t(dMs)) * ivectomats(ds, ds$X * ds$Me, g)) %*% (Mis * ds$X) +
      t(Gis^2 * ds$e * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$X * ds$Me, g)) %*% (dMs * ds$X)
    A43vecs[g] <- t(Gis^2 * ds$e * dMs * Mis * recD2is) %*% (recD3is * ivectomats(ds, ds$X, g)) %*% (Mis^2 * ds$X * ds$Me) -
      t(Gis^2 * ds$e * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$X, g)) %*% (Mis * ds$X * ds$Me) -
      t(Gis^2 * ds$e * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$X * dMs, g)) %*% (Mis * ds$X * ds$Me) +
      t(Gis^2 * ds$e * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$X * dMs, g)) %*% (ds$X * ds$Me)
    A44vecs[g] <- t(ds$X^2 * dMs) %*% (recD2s * Psoff^2) %*% (ds$Me * ds$e) -
      t(ds$X^2 * ds$Me) %*% (recD2s * Psoff^2) %*% (Mis * ds$e)

    A51vecs[g] <- t(Gis^2 * ds$e * dMs * ds$MX) %*% (recD3is * ivectomats(ds, ds$e, g)) %*% (Mis * ds$X) -
      t(Gis^2 * ds$e * ds$MX) %*% (recD3is * Ms * (ones %x% t(ds$X)) * ivectomats(ds, ds$e, g)) %*% (Mis)
    A52vecs[g] <- t(Gis^2 * ds$e * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$e * ds$MX, g)) %*% (Mis^2 * ds$X) -
      t(Gis^2 * ds$e * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$e * ds$MX, g)) %*% (Mis * ds$X) -
      t(Gis^2 * ds$e * dMs * Mis * recD2is) %*% (recD3is * (ones %*% t(dMs)) * ivectomats(ds, ds$e * ds$MX, g)) %*% (Mis * ds$X) +
      t(Gis^2 * ds$e * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$e * ds$MX, g)) %*% (dMs * ds$X)
    A53vecs[g] <- t(Gis^2 * ds$e * dMs * Mis * recD2is) %*% (recD3is * ivectomats(ds, ds$e, g)) %*% (Mis^2 * ds$X * ds$MX) -
      t(Gis^2 * ds$e * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$e, g)) %*% (Mis * ds$X * ds$MX) -
      t(Gis^2 * ds$e * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$e * dMs, g)) %*% (Mis * ds$X * ds$MX) +
      t(Gis^2 * ds$e * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$e * dMs, g)) %*% (ds$X * ds$MX)
    A54vecs[g] <- t(ds$X * ds$e * dMs) %*% (recD2s * Psoff^2) %*% (ds$MX * ds$e) -
      t(ds$X * ds$e * ds$MX) %*% (recD2s * Psoff^2) %*% (Mis * ds$e)
  }

  ret <- (A11vecs - A12vecs - A13vecs + A14vecs + A15vecs) +
    2 * (A21vecs - A22vecs - A23vecs + A24vecs + A25vecs) +
    (A31vecs - A32vecs - A33vecs + A34vecs + A35vecs) -
    (A41vecs + A42vecs + A43vecs + A44vecs) -
    (A51vecs + A52vecs + A53vecs + A54vecs)

  sum(ret)
}

#' Compute L3O Variance for Score Statistic (No Covariates, General Symmetric P)
#'
#' @description
#' Calculates the Leave-Three-Out (L3O) variance estimator (\eqn{\hat{V}_{LM}}) for the score statistic
#' in the "No Covariates" setting. Unlike the group-optimized functions, this implementation
#' loops over every observation \eqn{i}, making it suitable for any design where the
#' weighting matrix is symmetric (\eqn{G = P}), even if it lacks a strict block-diagonal structure.
#'
#' @param X Numeric vector of length n. The endogenous regressor.
#' @param e Numeric vector of length n. Residuals under the null hypothesis.
#' @param P Matrix of dimension n x n. The projection matrix (must be symmetric).
#'
#' @details
#' This function implements the variance estimator \eqn{\hat{V}_{LM}} under the assumption
#' that there are no covariates, implying \eqn{G = P} and \eqn{P} is symmetric.
#'
#' \strong{Specifically:}
#' The function iterates through each observation \eqn{i} (1 to \eqn{N}). Inside the loop, it:
#' \enumerate{
#'   \item Extracts the \eqn{i}-th column of \eqn{P} and \eqn{M}.
#'   \item Computes the L3O determinant adjustments \eqn{D_{ijk}} for all \eqn{k}.
#'   \item Calculates all five variance components (\eqn{A_1 \dots A_5}) simultaneously using
#'   pre-computed vector products to maximize efficiency.
#' }
#'
#' This implementation is slower than \code{L3Ovar_gloop_nocov} for simple grouped designs but
#' is more general.
#'
#' @return Scalar. The estimated variance \eqn{\hat{V}_{LM}}.
#'
#' @export
L3Ovar_iloop_nocov <- function(X, e, P) {
  n <- length(X)
  M <- diag(n) - P

  # Stuff that can be calculated once
  dM <- matrix(diag(M), ncol = 1) # force column vector
  D2 <- dM %*% t(dM) - M * M
  onesN <- matrix(rep(1, n), ncol = 1)
  recD2 <- 1 / D2
  diag(recD2) <- 0
  Me <- M %*% e
  MX <- M %*% X
  Poff <- P - diag(diag(P))

  A1vec <- A2vec <- A3vec <- A4vec <- A5vec <- rep(0, n)
  for (i in 1:n) {
    # Calculation conditioned on i
    Mi <- matrix(M[, i], ncol = 1) # force column vector
    Pi <- matrix(P[, i], ncol = 1) # force column vector
    D3i <- M[i, i] * D2 - (dM %*% t(Mi)^2 + Mi^2 %*% t(dM) - 2 * M * (Mi %*% t(Mi)))
    Di <- matrix(D2[, i], ncol = 1)
    D2D3i <- D2 / D3i
    D2D3i[i, ] <- 0
    D2D3i[, i] <- 0
    diag(D2D3i) <- 0
    recD3i <- 1 / D3i
    recD3i[i, ] <- 0
    recD3i[, i] <- 0
    diag(recD3i) <- 0
    recD2i <- matrix(1 / D2[, i], ncol = 1)
    recD2i[i] <- 0
    Poffi <- Poff[, i]
    tMie <- t(Mi) %*% e
    tMiX <- t(Mi) %*% X
    Pie <- Pi * e
    PiX <- Pi * X
    MXe <- MX * e
    Mie <- Mi * e
    MiX <- Mi * X
    MrecD3iPiMie <- (M * recD3i) %*% (Pie * Mi)
    recD3iMMiXMX <- (recD3i * M) %*% (MiX * MX)
    recD3iMMiXMe <- (recD3i * M) %*% (MiX * Me)
    recD3iPiedM <- recD3i %*% (Pie * dM)
    Pi2eMirecD2irecD3iMM <- t(Pi^2 * e * Mi * recD2i) %*% (recD3i * M * M)
    Pi2edMMirecD2irecD3i <- t(Pi^2 * e * dM * Mi * recD2i) %*% (recD3i)
    PiXD2D3i <- t(PiX) %*% D2D3i


    A11i <- PiXD2D3i %*% (PiX) * (tMie)
    A12i <- t((Me) * PiX * Mi) %*% recD3i %*% (Poffi * X * dM) -
      t(PiX * Mi) %*% (recD3i * M) %*% (Poffi * X * (Me))
    A13i <- t(dM * PiX) %*% ((onesN %x% t(Me)) * recD3i) %*% (Pi * MiX)
    A14i <- t((Me) * PiX) %*% (M * recD3i) %*% (Pi * MiX)
    A15i <- (tMie) * (t(Poffi^2 * recD2i) %*% (dM * X^2)) -
      (t(Poffi^2 * Mi * recD2i) %*% (Me * X^2))

    A21i <- PiXD2D3i %*% (Pie) * (tMiX)
    A22i <- t((MX) * PiX * Mi) %*% recD3iPiedM -
      t(PiX * Mi) %*% (recD3i * M) %*% (Pie * (MX))
    A23i <- t(dM * PiX) %*% ((onesN %x% t(MX)) * recD3i) %*% (Pie * Mi)
    A24i <- t((MX) * PiX) %*% MrecD3iPiMie
    A25i <- (tMiX) * (t(Poffi^2 * recD2i) %*% (dM * X * e)) -
      (t(Poffi^2 * Mi * recD2i) %*% (MXe * X))

    A31i <- t(Pie) %*% D2D3i %*% (Pie) * (tMiX)
    A32i <- t((MX) * Pie * Mi) %*% recD3iPiedM -
      t(Pie * Mi) %*% (recD3i * M) %*% (Pie * (MX))
    A33i <- t(dM * Pie) %*% ((onesN %x% t(MX)) * recD3i) %*% (Pie * Mi)
    A34i <- t((MX) * Pie) %*% MrecD3iPiMie
    A35i <- (tMiX) * (t(Poffi^2 * recD2i) %*% (dM * e * e)) -
      (t(Poffi^2 * Mi * recD2i) %*% (MXe * e))

    A41i <- t(Pi^2 * e * dM * recD2i * Me) %*% ((onesN %x% t(Di)) * recD3i) %*% (MiX) -
      t(Pi^2 * e * recD2i * Me) %*% ((onesN %x% t(Di)) * recD3i * M * (onesN %x% t(X))) %*% (Mi)
    A42i <- t(Pi^2 * e * dM * recD2i) %*% (recD3i * M) %*% (Mi^2 * X) -
      Pi2eMirecD2irecD3iMM %*% (MiX) -
      t(Pi^2 * e * dM * Mi * recD2i) %*% ((onesN %x% t(dM)) * recD3i) %*% (MiX) +
      t(Pi^2 * e * Mi^2 * recD2i) %*% (recD3i * M) %*% (dM * X)
    A43i <- Pi2edMMirecD2irecD3i %*% (Mi^2 * X * Me) -
      t(Pi^2 * e * Mi^2 * recD2i) %*% recD3iMMiXMe -
      M[i, i] * t(Pi^2 * e * dM * recD2i) %*% recD3iMMiXMe +
      M[i, i] * Pi2eMirecD2irecD3iMM %*% (Me * X)
    A44i <- X[i] * M[i, i] * (t(Poffi^2 * recD2i) %*% (Me * e)) -
      X[i] * (tMie) * (t(Poffi^2 * recD2i) %*% (Mie))

    A51i <- t(Pi^2 * e * dM * recD2i * MX) %*% ((onesN %x% t(Di)) * recD3i) %*% (Mi * X) -
      t(Pi^2 * e * recD2i * MX) %*% ((onesN %x% t(Di)) * recD3i * M * (onesN %x% t(X))) %*% (Mi)
    A53i <- Pi2edMMirecD2irecD3i %*% (Mi^2 * X * MX) -
      t(Pi^2 * e * Mi^2 * recD2i) %*% recD3iMMiXMX -
      M[i, i] * t(Pi^2 * e * dM * recD2i) %*% recD3iMMiXMX +
      M[i, i] * Pi2eMirecD2irecD3iMM %*% (MX * X)
    A54i <- X[i] * M[i, i] * (t(Poffi^2 * recD2i) %*% (MXe)) -
      X[i] * (tMiX) * (t(Poffi^2 * recD2i) %*% (Mie))

    A1vec[i] <- A11i - A12i - A13i + A14i + A15i
    A2vec[i] <- A21i - A22i - A23i + A24i + A25i
    A3vec[i] <- A31i - A32i - A33i + A34i + A35i
    A4vec[i] <- A41i + A42i * (M[i, ] %*% e) + A43i + A44i
    A5vec[i] <- A51i + A42i * (M[i, ] %*% X) + A53i + A54i

    if (i %% 10 == 0) cat(i / n, " ")
  }

  sum(A1vec * e + 2 * A2vec * e + A3vec * X - A4vec * X - A5vec * e)
}
