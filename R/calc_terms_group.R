#' Compute Group-Level A1-Type Variance Component
#'
#' @description
#' Calculates the \eqn{A_1} variance component by exploiting the discrete group structure of
#' the instruments. Instead of iterating over every observation \eqn{i}, this function iterates
#' over instrument groups (e.g., judges) and covariate strata, using representative observations
#' to efficiently compute the Leave-Three-Out (L3O) variance sums.
#'
#'
#' @param df Data frame. Contains the data vectors and grouping variables.
#' @param group Column name (unquoted). The instrument grouping variable (defines clusters).
#' @param groupW Column name (unquoted). The covariate stratification variable (defines strata).
#' @param ipos Column name (unquoted). Outer summation weight \eqn{v^{(I)}}.
#' @param jpos Column name (unquoted). Inner summation term \eqn{v^{(J)}}.
#' @param kpos Column name (unquoted). Inner summation term \eqn{v^{(K)}}.
#' @param lpos Column name (unquoted). Bias correction term \eqn{v^{(L)}}.
#' @param noisy Logical. If \code{TRUE}, prints progress of the stratification loop.
#'   Defaults to \code{FALSE}.
#'
#' @details
#' This function implements the \eqn{A_1} estimator (Signal Variance) specifically for
#' discrete instrument settings where projection matrices are block-diagonal with respect
#' to \code{groupW} and constant within \code{group}.
#'
#' \strong{The function works as follow:}
#' \enumerate{
#'   \item Splits data by \code{groupW}. Projection matrices \eqn{P} and \eqn{G}
#'   are constructed locally within each stratum to reduce computational cost.
#'   \item Removes groups with fewer than 3 observations (required for L3O validity).
#'   \item Iterates over groups \eqn{g}. Uses \code{\link{ivectomats}} to
#'   compute the L2O-adjusted sum of outer weights for all members of group \eqn{g} simultaneously.
#' }
#'
#' @return Numeric scalar. The total variance component.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @noRd
A1type_sum <- function(df, group, groupW, ipos, jpos, kpos, lpos, noisy = FALSE) {

  # 0. Setup with zero vectors
  A11vecs <- A12vecs <- A13vecs <- A14vecs <- A15vecs <- rep(0, max(df$group))

  # 1. Evaluate inputs from dataframe
  df$ipos <- eval(substitute(ipos), df)
  df$jpos <- eval(substitute(jpos), df)
  df$kpos <- eval(substitute(kpos), df)
  df$lpos <- eval(substitute(lpos), df)

  # 2. Looping through groupW
  iteration <- 1
  for (s in unique(df$groupW)) {
    ds <- df[df$groupW == s, ]

  # 2.1. Identify stratum groups within groupW and remove ones with < 3 obs
    ds$numingrp <- 0
    for (j in unique(ds$group)) {
      ds$numingrp[ds$group == j] <- sum(ds$group == j)
    }
    ds <- ds[ds$numingrp >= 3, ]

    # 3. Compute matrix P_Q for each stratum
    if (nrow(ds) == 0) {
      for (g in unique(ds$group)) {
        A11vecs[g] <- A12vecs[g] <- A13vecs[g] <- A14vecs[g] <- A15vecs[g] <- 0
      }
    } else {
      ZQ <- matrix(0, nrow = length(ds$group), ncol = length(unique(ds$group)))
      ds$groupidx <- Getgroupindex(ds, group)
      ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1
      # ZW <- matrix(1, nrow=length(ds$groupW), ncol=length(unique(ds$groupW)))

      PQ <- ZQ %*% solve(t(ZQ) %*% ZQ) %*% t(ZQ)
      # PW <- ZW %*% solve(t(ZW)%*% ZW) %*% t(ZW)

      ZWmat <- matrix(1, nrow = length(ds$groupW), ncol = length(ds$groupW))
      PW <- ZWmat / (length(ds$groupW))

      # calculate values specific to this subset
      Gs <- diag(1 / (diag(diag(nrow(ds))) - pmin(diag(PQ), .99))) %*% (PQ - diag(diag(PQ))) -
        diag(1 / (diag(diag(nrow(ds))) - pmin(diag(PW), .99))) %*% (PW - diag(diag(PW)))
      Ps <- PQ
      Ms <- diag(nrow(ds)) - Ps
      dMs <- matrix(diag(Ms), ncol = 1)
      D2s <- dMs %*% t(dMs) - Ms * Ms
      recD2s <- 1 / D2s
      diag(recD2s) <- 0

      # 4. Iterate through each group within stratum, find its contribution
      for (g in unique(ds$group)) {
        if (nrow(ds[ds$group == g, ]) <= 3) {
          A11vecs[g] <- A12vecs[g] <- A13vecs[g] <- A14vecs[g] <- A15vecs[g] <- 0
        } else {
          repidx <- min(which(ds$group == g)) # representative index
          Pis <- matrix(Ps[, repidx], ncol = 1)
          Pgs <- ifelse(Pis == 0, 0, 1) %*% matrix(1, ncol = length(Pis), nrow = 1)
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
          Mes <- matrix(ds$lpos, ncol = 1)

          A11vecs[g] <- t(ds$jpos * Gis) %*% (D2D3is * ivectomats(ds, ds$ipos * ds$lpos, g)) %*% (ds$kpos * Gis)
          A12vecs[g] <- t(ds$jpos * ds$lpos * Gis * Mis) %*% (recD3is * ivectomats(ds, ds$ipos, g)) %*% (Gis * ds$kpos * dMs) -
            t(ds$jpos * Gis * Mis) %*% (recD3is * Ms * ivectomats(ds, ds$ipos, g)) %*% (Gis * ds$kpos * ds$lpos)
          A13vecs[g] <- t(ds$jpos * dMs * Gis) %*% (recD3is * (ones %*% t(Mes)) * ivectomats(ds, ds$ipos, g)) %*% (ds$kpos * Gis * Mis)
          A14vecs[g] <- t(ds$jpos * ds$lpos * Gis) %*% (recD3is * Ms * ivectomats(ds, ds$ipos, g)) %*% (ds$kpos * Gis * Mis)
          A15vecs[g] <- t(ds$lpos * ds$ipos) %*% (recD2s * Gs^2 * Pgs) %*% (ds$jpos * ds$kpos * dMs) -
            t(ds$ipos) %*% (Gs^2 * Ms * recD2s * Pgs) %*% (ds$lpos * ds$jpos * ds$kpos)
        }
      }
    }
    if (noisy) {
      cat(iteration, "of", max(df$groupW), "done. ")
      iteration <- iteration + 1
    }
  }

  ret <- (A11vecs - A12vecs - A13vecs + A14vecs + A15vecs)

  sum(ret)
}

#' Compute Group-Level A4-Type Variance Component
#'
#' @description
#' Calculates the \eqn{A_4} variance component ("own-variance" bias correction) by exploiting
#' the discrete group structure of the instruments. This function iterates over instrument
#' groups (e.g., judges) and covariate strata to efficiently compute the bias arising from
#' squared diagonal weights \eqn{G_{ji}^2}.
#'
#' @param df Data frame. Contains the data vectors and grouping variables.
#' @param group Column name (unquoted). The instrument grouping variable.
#' @param groupW Column name (unquoted). The covariate stratification variable.
#' @param ipos Column name (unquoted). Outer summation weight \eqn{v^{(I)}}.
#' @param jpos Column name (unquoted). Inner summation term \eqn{v^{(J)}} weighted by \eqn{G_{ji}^2}.
#' @param kpos Column name (unquoted). Term interacting with leverage \eqn{v^{(K)}}.
#' @param lpos Column name (unquoted). Residual/bias interaction term \eqn{v^{(L)}}.
#' @param noisy Logical. If \code{TRUE}, prints progress of the stratification loop.
#'   Defaults to \code{FALSE}.
#'
#' @details
#' This function implements the \eqn{A_4} estimator for discrete instrument designs.
#' It estimates the positive bias component:
#' \deqn{S = \sum_{g} \sum_{i \in g} v_i^{(I)} \sum_{j \notin g} G_{ji}^2 \widehat{Var}(v_j) W_{ij}}
#'
#' \strong{The function works as follow:}
#' \enumerate{
#'   \item Computes projection matrices locally within \code{groupW} strata.
#'   \item For each group \eqn{g}, extracts the column
#'   vector of weights \eqn{G_{\cdot g}} (constant for all \eqn{i \in g}). Squaring this
#'   vector captures the "own-variance" influence \eqn{G_{jg}^2}.
#'   \item Uses \code{\link{ivectomats}} to aggregate outer weights
#'   over the group, applying L3O corrections efficiently.
#' }
#'
#' @return Numeric scalar.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @noRd
A4type_sum <- function(df, group, groupW, ipos, jpos, kpos, lpos, noisy = FALSE) {
  A41vecs <- A42vecs <- A43vecs <- A44vecs <- rep(0, max(df$group))

  df$ipos <- eval(substitute(ipos), df)
  df$jpos <- eval(substitute(jpos), df)
  df$kpos <- eval(substitute(kpos), df)
  df$lpos <- eval(substitute(lpos), df)

  iteration <- 1
  for (s in unique(df$groupW)) {
    ds <- df[df$groupW == s, ]
    ds$numingrp <- 0
    for (j in unique(ds$group)) {
      ds$numingrp[ds$group == j] <- sum(ds$group == j)
    }
    ds <- ds[ds$numingrp >= 3, ]
    if (nrow(ds) == 0) {
      for (g in unique(ds$group)) {
        A41vecs[g] <- A42vecs[g] <- A43vecs[g] <- A44vecs[g] <- 0
      }
    } else {
      ZQ <- matrix(0, nrow = length(ds$group), ncol = length(unique(ds$group)))
      ds$groupidx <- Getgroupindex(ds, group)
      ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1
      # ZW <- matrix(1, nrow=length(ds$groupW), ncol=length(unique(ds$groupW)))

      PQ <- ZQ %*% solve(t(ZQ) %*% ZQ) %*% t(ZQ)
      # PW <- ZW %*% solve(t(ZW)%*% ZW) %*% t(ZW)

      ZWmat <- matrix(1, nrow = length(ds$groupW), ncol = length(ds$groupW))
      PW <- ZWmat / (length(ds$groupW))

      # calculate values specific to this subset
      Gs <- diag(1 / (diag(diag(nrow(ds))) - pmin(diag(PQ), .99))) %*% (PQ - diag(diag(PQ))) -
        diag(1 / (diag(diag(nrow(ds))) - pmin(diag(PW), .99))) %*% (PW - diag(diag(PW)))
      Ps <- PQ
      Ms <- diag(nrow(ds)) - Ps
      dMs <- matrix(diag(Ms), ncol = 1)
      D2s <- dMs %*% t(dMs) - Ms * Ms
      recD2s <- 1 / D2s
      diag(recD2s) <- 0

      for (g in unique(ds$group)) {
        if (nrow(ds[ds$group == g, ]) <= 3) {
          A41vecs[g] <- A42vecs[g] <- A43vecs[g] <- A44vecs[g] <- 0
        } else {
          repidx <- min(which(ds$group == g)) # representative index
          Pis <- matrix(Ps[, repidx], ncol = 1)
          Pgs <- ifelse(Pis == 0, 0, 1) %*% matrix(1, ncol = length(Pis), nrow = 1)
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
          Mes <- matrix(ds$lpos, ncol = 1)

          A41vecs[g] <- t(Gis^2 * ds$jpos * dMs * ds$lpos * recD2is) %*% (D2D3is * ivectomats(ds, ds$ipos, g)) %*% (Mis * ds$kpos) -
            t(Gis^2 * ds$jpos * ds$lpos * recD2is) %*% (D2D3is * Ms * (ones %x% t(ds$kpos)) * ivectomats(ds, ds$ipos, g)) %*% (Mis)
          A42vecs[g] <- t(Gis^2 * ds$jpos * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$ipos * ds$lpos, g)) %*% (Mis^2 * ds$kpos) -
            t(Gis^2 * ds$jpos * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$ipos * ds$lpos, g)) %*% (Mis * ds$kpos) -
            t(Gis^2 * ds$jpos * dMs * Mis * recD2is) %*% (recD3is * (ones %*% t(dMs)) * ivectomats(ds, ds$ipos * ds$lpos, g)) %*% (Mis * ds$kpos) +
            t(Gis^2 * ds$jpos * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$ipos * ds$lpos, g)) %*% (dMs * ds$kpos)
          A43vecs[g] <- t(Gis^2 * ds$jpos * dMs * Mis * recD2is) %*% (recD3is * ivectomats(ds, ds$ipos, g)) %*% (Mis^2 * ds$kpos * ds$lpos) -
            t(Gis^2 * ds$jpos * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$ipos, g)) %*% (Mis * ds$kpos * ds$lpos) -
            t(Gis^2 * ds$jpos * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$ipos * dMs, g)) %*% (Mis * ds$kpos * ds$lpos) +
            t(Gis^2 * ds$jpos * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$ipos * dMs, g)) %*% (ds$kpos * ds$lpos)
          A44vecs[g] <- t(ds$ipos * ds$kpos * dMs) %*% (recD2s * Gs^2 * Pgs) %*% (ds$lpos * ds$jpos) -
            t(ds$ipos * ds$kpos * ds$lpos) %*% (recD2s * Gs^2 * Pgs) %*% (Mis * ds$jpos)
        }
      }
    }
    if (noisy) {
      cat(iteration, "of", max(df$groupW), "done. ")
      iteration <- iteration + 1
    }
  }

  ret <- (A41vecs + A42vecs + A43vecs + A44vecs)

  sum(ret)
}

#' Compute A1 Variance Component (No Covariates)
#'
#' @description
#' Calculates the \eqn{A_1} variance component (Signal Variance) for the specific case
#' where there are no covariates. This function exploits the resulting block-diagonal
#' structure of the projection matrix to stratify calculations by instrument group,
#' significantly improving computational efficiency.
#'
#' @param df Data frame. Contains the data vectors and the grouping variable.
#' @param groupZ Column name (unquoted). The instrument grouping variable.
#' @param ipos Column name (unquoted). Outer summation weight \eqn{v^{(I)}}.
#' @param jpos Column name (unquoted). Inner summation term \eqn{v^{(J)}}.
#' @param kpos Column name (unquoted). Inner summation term \eqn{v^{(K)}}.
#' @param lpos Column name (unquoted). Bias correction term \eqn{v^{(L)}}.
#' @param noisy Logical. If \code{TRUE}, prints progress of the group iteration.
#'   Defaults to \code{FALSE}.
#'
#' @details
#' This function implements the \eqn{A_1} estimator for the "Many Means" model or
#' simple One-Way Layout designs without additional controls.
#'
#' \strong{The function works as follow:}
#' \itemize{
#'   \item Because there are no covariates linking observations
#'   across groups, the projection matrix \eqn{P} is strictly block-diagonal. The code
#'   iterates through each group \eqn{g}, treating it as an independent dataset.
#'   \item The UJIVE weighting matrix simplifies to \eqn{G = U(P_Z)},
#'   as there is no covariate projection \eqn{P_W} to partial out.
#' }
#'
#' The term computed corresponds to the variance of the first-stage fitted values (signal),
#' corrected for finite-sample bias using the Leave-Three-Out (L3O) approach.
#'
#' @return Numeric scalar.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @noRd
A1type_sum_nocov <- function(df, groupZ, ipos, jpos, kpos, lpos, noisy = FALSE) {
  df$groupZ <- eval(substitute(groupZ), df)
  A11vecs <- A12vecs <- A13vecs <- A14vecs <- A15vecs <- rep(0, max(df$groupZ))

  df$ipos <- eval(substitute(ipos), df)
  df$jpos <- eval(substitute(jpos), df)
  df$kpos <- eval(substitute(kpos), df)
  df$lpos <- eval(substitute(lpos), df)

  iteration <- 1
  for (s in unique(df$groupZ)) {
    ds <- df[df$groupZ == s, ]
    ZQ <- matrix(0, nrow = length(ds$groupZ), ncol = length(unique(ds$groupZ)))
    ds$groupidx <- Getgroupindex(ds, groupZ)
    ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1

    PQ <- ZQ %*% solve(t(ZQ) %*% ZQ) %*% t(ZQ)

    # calculate values specific to this subset
    Gs <- diag(1 / (diag(diag(nrow(ds))) - pmin(diag(PQ), .99))) %*% (PQ - diag(diag(PQ)))
    Ps <- PQ
    Ms <- diag(nrow(ds)) - Ps
    dMs <- matrix(diag(Ms), ncol = 1)
    D2s <- dMs %*% t(dMs) - Ms * Ms
    recD2s <- 1 / D2s
    diag(recD2s) <- 0

    g <- s

    repidx <- min(which(ds$groupZ == g)) # representative index
    Pis <- matrix(Ps[, repidx], ncol = 1)
    Pgs <- ifelse(Pis == 0, 0, 1) %*% matrix(1, ncol = length(Pis), nrow = 1)
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
    Mes <- matrix(ds$lpos, ncol = 1)

    A11vecs[g] <- t(ds$jpos * Gis) %*% (D2D3is * ivectomats(ds, ds$ipos * ds$lpos, g)) %*% (ds$kpos * Gis)
    A12vecs[g] <- t(ds$jpos * ds$lpos * Gis * Mis) %*% (recD3is * ivectomats(ds, ds$ipos, g)) %*% (Gis * ds$kpos * dMs) -
      t(ds$jpos * Gis * Mis) %*% (recD3is * Ms * ivectomats(ds, ds$ipos, g)) %*% (Gis * ds$kpos * ds$lpos)
    A13vecs[g] <- t(ds$jpos * dMs * Gis) %*% (recD3is * (ones %*% t(Mes)) * ivectomats(ds, ds$ipos, g)) %*% (ds$kpos * Gis * Mis)
    A14vecs[g] <- t(ds$jpos * ds$lpos * Gis) %*% (recD3is * Ms * ivectomats(ds, ds$ipos, g)) %*% (ds$kpos * Gis * Mis)
    A15vecs[g] <- t(ds$lpos * ds$ipos) %*% (recD2s * Gs^2 * Pgs) %*% (ds$jpos * ds$kpos * dMs) -
      t(ds$ipos) %*% (Gs^2 * Ms * recD2s * Pgs) %*% (ds$lpos * ds$jpos * ds$kpos)
    if (noisy) {
      cat(iteration, "of", max(df$groupZ), "done. ")
      iteration <- iteration + 1
    }
  }

  ret <- (A11vecs - A12vecs - A13vecs + A14vecs + A15vecs)

  sum(ret)
}

#' Compute A4 Variance Component (No Covariates)
#'
#' @description
#' Calculates the \eqn{A_4} variance component ("own-variance" bias correction) for the
#' specific case where there are no covariates. This function exploits the block-diagonal
#' structure of the projection matrix in the absence of global controls to stratify
#' calculations by instrument group, optimizing performance.
#'
#' @param df Data frame. Contains the data vectors and the grouping variable.
#' @param groupZ Column name (unquoted). The instrument grouping variable.
#' @param ipos Column name (unquoted). Outer summation weight \eqn{v^{(I)}}.
#' @param jpos Column name (unquoted). Inner summation term \eqn{v^{(J)}} weighted by \eqn{G_{ji}^2}.
#' @param kpos Column name (unquoted). Term interacting with leverage \eqn{v^{(K)}}.
#' @param lpos Column name (unquoted). Residual/bias interaction term \eqn{v^{(L)}}.
#' @param noisy Logical. If \code{TRUE}, prints progress of the group iteration.
#'   Defaults to \code{FALSE}.
#'
#' @details
#' This function implements the \eqn{A_4} estimator for the "Many Means" model. It estimates
#' the positive bias arising from the squared diagonal weights \eqn{G_{ji}^2}, which must be
#' subtracted from the final variance estimate.
#'
#' \strong{The function works as follow:}
#' \itemize{
#'   \item Since no covariates link observations across groups,
#'   the dataset is split by \code{groupZ}, and calculations are performed locally on subsets.
#'   \item The UJIVE weight \eqn{G} is computed as \eqn{U(P_Z)},
#'   omitting the covariate projection term \eqn{U(P_W)}.
#' }
#'
#' @return Numeric scalar.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @noRd
A4type_sum_nocov <- function(df, groupZ, ipos, jpos, kpos, lpos, noisy = FALSE) {
  df$groupZ <- eval(substitute(groupZ), df)
  A41vecs <- A42vecs <- A43vecs <- A44vecs <- rep(0, max(df$groupZ))

  df$ipos <- eval(substitute(ipos), df)
  df$jpos <- eval(substitute(jpos), df)
  df$kpos <- eval(substitute(kpos), df)
  df$lpos <- eval(substitute(lpos), df)

  iteration <- 1
  for (s in unique(df$groupZ)) {
    ds <- df[df$groupZ == s, ]
    ZQ <- matrix(0, nrow = length(ds$groupZ), ncol = length(unique(ds$groupZ)))
    ds$groupidx <- Getgroupindex(ds, groupZ)
    ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1

    PQ <- ZQ %*% solve(t(ZQ) %*% ZQ) %*% t(ZQ)

    # calculate values specific to this subset
    Gs <- diag(1 / (diag(diag(nrow(ds))) - pmin(diag(PQ), .99))) %*% (PQ - diag(diag(PQ)))
    Ps <- PQ
    Ms <- diag(nrow(ds)) - Ps
    dMs <- matrix(diag(Ms), ncol = 1)
    D2s <- dMs %*% t(dMs) - Ms * Ms
    recD2s <- 1 / D2s
    diag(recD2s) <- 0

    g <- s

    repidx <- min(which(ds$groupZ == g)) # representative index
    Pis <- matrix(Ps[, repidx], ncol = 1)
    Pgs <- ifelse(Pis == 0, 0, 1) %*% matrix(1, ncol = length(Pis), nrow = 1)
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
    Mes <- matrix(ds$lpos, ncol = 1)

    A41vecs[g] <- t(Gis^2 * ds$jpos * dMs * ds$lpos * recD2is) %*% (D2D3is * ivectomats(ds, ds$ipos, g)) %*% (Mis * ds$kpos) -
      t(Gis^2 * ds$jpos * ds$lpos * recD2is) %*% (D2D3is * Ms * (ones %x% t(ds$kpos)) * ivectomats(ds, ds$ipos, g)) %*% (Mis)
    A42vecs[g] <- t(Gis^2 * ds$jpos * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$ipos * ds$lpos, g)) %*% (Mis^2 * ds$kpos) -
      t(Gis^2 * ds$jpos * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$ipos * ds$lpos, g)) %*% (Mis * ds$kpos) -
      t(Gis^2 * ds$jpos * dMs * Mis * recD2is) %*% (recD3is * (ones %*% t(dMs)) * ivectomats(ds, ds$ipos * ds$lpos, g)) %*% (Mis * ds$kpos) +
      t(Gis^2 * ds$jpos * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$ipos * ds$lpos, g)) %*% (dMs * ds$kpos)
    A43vecs[g] <- t(Gis^2 * ds$jpos * dMs * Mis * recD2is) %*% (recD3is * ivectomats(ds, ds$ipos, g)) %*% (Mis^2 * ds$kpos * ds$lpos) -
      t(Gis^2 * ds$jpos * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$ipos, g)) %*% (Mis * ds$kpos * ds$lpos) -
      t(Gis^2 * ds$jpos * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$ipos * dMs, g)) %*% (Mis * ds$kpos * ds$lpos) +
      t(Gis^2 * ds$jpos * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$ipos * dMs, g)) %*% (ds$kpos * ds$lpos)
    A44vecs[g] <- t(ds$ipos * ds$kpos * dMs) %*% (recD2s * Gs^2 * Pgs) %*% (ds$lpos * ds$jpos) -
      t(ds$ipos * ds$kpos * ds$lpos) %*% (recD2s * Gs^2 * Pgs) %*% (Mis * ds$jpos)
    if (noisy) {
      cat(iteration, "of", max(df$groupZ), "done. ")
      iteration <- iteration + 1
    }
  }

  ret <- (A41vecs + A42vecs + A43vecs + A44vecs)

  sum(ret)
}
