#' Compute Quadratic Coefficients for CI (No Covariates)
#'
#' @description
#' Calculates the coefficients \eqn{(a, b, c)} of the quadratic inequality \eqn{a\beta^2 + b\beta + c \leq 0}
#' used to construct confidence intervals for \eqn{\beta} in the grouped design setting (no covariates).
#' This function inverts the UJIVE/LIML score test using optimized variance estimators for
#' block-diagonal designs.
#'
#' @param df Data frame. Contains the observable variables \eqn{X, Y} and their projections.
#' @param groupZ Column name (unquoted). The instrument grouping variable.
#' @param X Column name (unquoted). The endogenous regressor.
#' @param Y Column name (unquoted). The outcome variable.
#' @param MX Column name (unquoted). Leverage-adjusted regressor (\eqn{M X}).
#' @param MY Column name (unquoted). Leverage-adjusted outcome (\eqn{M Y}).
#' @param q Numeric scalar. Critical value for the test inversion (e.g., \eqn{1.96^2}).
#'   Defaults to \code{qnorm(.975)^2} (approx. 3.84) for a 95 percent confidence interval.
#' @param noisy Logical. If \code{TRUE}, prints progress dots during calculation.
#'   Defaults to \code{FALSE}.
#'
#' @details
#' This function performs the same logic as \code{\link{GetCIcoef}} but is specifically
#' appropriate for designs where instruments form mutually exclusive groups
#' and no global covariates link them.
#'
#' The returned coefficients define the confidence set:
#' \deqn{\{ \beta : (P_{XX}^2 - q C_2)\beta^2 + (-2 P_{XY} P_{XX} - q C_1)\beta + (P_{XY}^2 - q C_0) \leq 0 \}}
#'
#'
#' @return Numeric vector of length 3 containing \code{c(a, b, c)}.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @export
GetCIcoef_nocov <- function(df, groupZ, X, Y, MX, MY, q = qnorm(.975)^2, noisy = FALSE) {
  df$groupZ <- eval(substitute(groupZ), df)
  df$X <- eval(substitute(X), df)
  df$Y <- eval(substitute(Y), df)
  df$MX <- eval(substitute(MX), df)
  df$MY <- eval(substitute(MY), df)
  ## Calculate Components
  C0 <- A1type_sum_nocov(df, groupZ, ipos = Y, jpos = X, kpos = X, lpos = MY, noisy = noisy) +
    2 * A1type_sum_nocov(df, groupZ, ipos = Y, jpos = X, kpos = Y, lpos = MX, noisy = noisy) +
    A1type_sum_nocov(df, groupZ, ipos = X, jpos = Y, kpos = Y, lpos = MX, noisy = noisy) -
    A4type_sum_nocov(df, groupZ, ipos = X, jpos = Y, kpos = X, lpos = MY, noisy = noisy) -
    A4type_sum_nocov(df, groupZ, ipos = Y, jpos = Y, kpos = X, lpos = MX, noisy = noisy)
  C1 <- -(A1type_sum_nocov(df, groupZ, ipos = X, jpos = X, kpos = X, lpos = MY, noisy = noisy) +
            2 * A1type_sum_nocov(df, groupZ, ipos = X, jpos = X, kpos = Y, lpos = MX, noisy = noisy) +
            A1type_sum_nocov(df, groupZ, ipos = X, jpos = X, kpos = Y, lpos = MX, noisy = noisy) -
            A4type_sum_nocov(df, groupZ, ipos = X, jpos = X, kpos = X, lpos = MY, noisy = noisy) -
            A4type_sum_nocov(df, groupZ, ipos = X, jpos = Y, kpos = X, lpos = MX, noisy = noisy) +
            A1type_sum_nocov(df, groupZ, ipos = Y, jpos = X, kpos = X, lpos = MX, noisy = noisy) +
            2 * A1type_sum_nocov(df, groupZ, ipos = Y, jpos = X, kpos = X, lpos = MX, noisy = noisy) +
            A1type_sum_nocov(df, groupZ, ipos = X, jpos = Y, kpos = X, lpos = MX, noisy = noisy) -
            A4type_sum_nocov(df, groupZ, ipos = X, jpos = Y, kpos = X, lpos = MX, noisy = noisy) -
            A4type_sum_nocov(df, groupZ, ipos = Y, jpos = X, kpos = X, lpos = MX, noisy = noisy))
  C2 <- A1type_sum_nocov(df, groupZ, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy) +
    2 * A1type_sum_nocov(df, groupZ, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy) +
    A1type_sum_nocov(df, groupZ, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy) -
    A4type_sum_nocov(df, groupZ, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy) -
    A4type_sum_nocov(df, groupZ, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy)
  PXY <- GetLM_nocov(df, X, Y, groupZ)
  PXX <- GetLM_nocov(df, X, X, groupZ)


  acon <- PXX^2 - q * C2
  bcon <- -2 * PXY * PXX - q * C1
  ccon <- PXY^2 - q * C0

  c(acon, bcon, ccon)
}

#' Compute Quadratic Coefficients for CI (General Design)
#'
#' @description
#' Computes the coefficients \eqn{(a, b, c)} of the quadratic inequality \eqn{a\beta^2 + b\beta + c \leq 0}
#' used to construct confidence intervals for \eqn{\beta} in general instrumental variable designs.
#' This function supports asymmetric weighting matrices (e.g., UJIVE with continuous covariates)
#' by using the fully generalized \code{_iloop} variance estimators.
#'
#' @param df Data frame. Contains the observable variables \eqn{X, Y} and their projections.
#' @param P Matrix of dimension n x n. The full projection matrix \eqn{P}.
#' @param G Matrix of dimension n x n. The UJIVE weighting matrix \eqn{G}.
#' @param X Column name (unquoted). The endogenous regressor.
#' @param Y Column name (unquoted). The outcome variable.
#' @param MX Column name (unquoted). Leverage-adjusted regressor (\eqn{M X}).
#' @param MY Column name (unquoted). Leverage-adjusted outcome (\eqn{M Y}).
#' @param Z Matrix of instruments.
#' @param W Matrix of covariates.
#' @param q Numeric scalar. Critical value for the test inversion (typically \eqn{1.96^2}).
#'   Defaults to \code{qnorm(.975)^2} (approx. 3.84) for a 95 percent confidence interval.
#' @param noisy Logical. If \code{TRUE}, prints progress dots during calculation.
#'   Defaults to \code{FALSE}.
#'
#' @details
#' This is the most general version of the confidence interval coefficient calculator.
#' It does not assume symmetry of the weighting matrix \eqn{G}. Consequently, it computes
#' all five distinct variance components (\eqn{A_1} through \eqn{A_5}) for each term in the
#' polynomial expansion of \eqn{\hat{V}(\beta)}.
#'
#' The returned coefficients define the curvature and position of the confidence set parabola:
#'
#'
#' The function explicitly constructs the diagonal adjustments for the UJIVE signal calculation
#' using the matrices \eqn{Z} and \eqn{W}.
#'
#' @return Numeric vector of length 3 containing \code{c(a, b, c)}.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @export
GetCIcoef_iloop <- function(df, P, G, X, Y, MX, MY, Z, W, q = qnorm(.975)^2, noisy = FALSE) {
  df$X <- eval(substitute(X), df)
  df$Y <- eval(substitute(Y), df)
  df$MX <- eval(substitute(MX), df)
  df$MY <- eval(substitute(MY), df)
  ## Calculate Components
  C0 <- A1type_iloop_sum(df, P, G, ipos = Y, jpos = X, kpos = X, lpos = MY, noisy = noisy) +
    2 * A2type_iloop_sum(df, P, G, ipos = Y, jpos = X, kpos = Y, lpos = MX, noisy = noisy) +
    A3type_iloop_sum(df, P, G, ipos = X, jpos = Y, kpos = Y, lpos = MX, noisy = noisy) -
    A4type_iloop_sum(df, P, G, ipos = X, jpos = Y, kpos = X, lpos = MY, noisy = noisy) -
    A5type_iloop_sum(df, P, G, ipos = Y, jpos = Y, kpos = X, lpos = MX, noisy = noisy)
  C1 <- -(A1type_iloop_sum(df, P, G, ipos = X, jpos = X, kpos = X, lpos = MY, noisy = noisy) +
            2 * A2type_iloop_sum(df, P, G, ipos = X, jpos = X, kpos = Y, lpos = MX, noisy = noisy) +
            A3type_iloop_sum(df, P, G, ipos = X, jpos = X, kpos = Y, lpos = MX, noisy = noisy) -
            A4type_iloop_sum(df, P, G, ipos = X, jpos = X, kpos = X, lpos = MY, noisy = noisy) -
            A5type_iloop_sum(df, P, G, ipos = X, jpos = Y, kpos = X, lpos = MX, noisy = noisy) +
            A1type_iloop_sum(df, P, G, ipos = Y, jpos = X, kpos = X, lpos = MX, noisy = noisy) +
            2 * A2type_iloop_sum(df, P, G, ipos = Y, jpos = X, kpos = X, lpos = MX, noisy = noisy) +
            A3type_iloop_sum(df, P, G, ipos = X, jpos = Y, kpos = X, lpos = MX, noisy = noisy) -
            A4type_iloop_sum(df, P, G, ipos = X, jpos = Y, kpos = X, lpos = MX, noisy = noisy) -
            A5type_iloop_sum(df, P, G, ipos = Y, jpos = X, kpos = X, lpos = MX, noisy = noisy))
  C2 <- A1type_iloop_sum(df, P, G, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy) +
    2 * A2type_iloop_sum(df, P, G, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy) +
    A3type_iloop_sum(df, P, G, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy) -
    A4type_iloop_sum(df, P, G, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy) -
    A5type_iloop_sum(df, P, G, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy)

  ## Generate dPQ
  n <- nrow(df)
  Q <- cbind(Z, W)
  dPQ <- dPW <- rep(0, n)
  WW <- t(W) %*% W
  WWinv <- solve(WW)
  QQ <- t(Q) %*% Q
  QQinv <- solve(QQ)
  for (i in 1:n) {
    dPW[i] <- W[i, ] %*% WWinv %*% W[i, ]
    dPQ[i] <- Q[i, ] %*% QQinv %*% Q[i, ]
  }
  IdPW <- pmax(1 - dPW, .01)
  IdPQ <- pmax(1 - dPQ, .01)

  PXY <- GetLM_WQ(df, IdPW, IdPQ, dPW, dPQ, W, Q, X, Y)
  PXX <- GetLM_WQ(df, IdPW, IdPQ, dPW, dPQ, W, Q, X, X)


  acon <- PXX^2 - q * C2
  bcon <- -2 * PXY * PXX - q * C1
  ccon <- PXY^2 - q * C0

  c(acon, bcon, ccon)
}

#' Compute Quadratic Coefficients for CI (No Covariates, General P)
#'
#' @description
#' Calculates the coefficients \eqn{(a, b, c)} for the confidence interval quadratic inequality
#' in the "No Covariates" setting (\eqn{G = P}), but for \strong{general symmetric projection matrices}.
#' This function performs a highly optimized single-pass loop to compute all polynomial coefficients
#' of the variance estimator \eqn{\hat{V}(\beta)} simultaneously.
#'
#' @param X Numeric vector of length n. The endogenous regressor.
#' @param Y Numeric vector of length n. The outcome variable.
#' @param P Matrix of dimension n x n. The symmetric projection matrix.
#' @param q Numeric scalar. The critical value for the test inversion (typically \eqn{1.96^2}).
#'   Defaults to \code{qnorm(0.975)^2}.
#' @param noisy Logical. If \code{TRUE}, prints progress through the N loops.
#'   Defaults to \code{FALSE}.
#'
#' @details
#' This function is the solver for generic symmetric designs.
#' It inverts the test statistic:
#' \deqn{\frac{(P_{XY} - \beta P_{XX})^2}{C_2 \beta^2 + C_1 \beta + C_0} \leq q}
#'
#' \strong{Optimization:}
#' Rather than calling variance estimators multiple times, it decomposes the variance formula
#' into geometric components (depending only on \eqn{P}) and data components (\eqn{X, Y}).
#' It iterates through observations \eqn{i} once, accumulating the weighted contributions
#' for \eqn{C_0}, \eqn{C_1}, and \eqn{C_2} in parallel.
#'
#' @return Numeric vector of length 3 containing \code{c(a, b, c)}.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @export
GetCIcoef_iloop_nocov <- function(X, Y, P, q = qnorm(0.975)^2, noisy = FALSE) {
  n <- length(X)
  M <- diag(n) - P
  MM <- M * M

  # Stuff that can be calculated once
  dM <- matrix(diag(M), ncol = 1) # force column vector
  D2 <- dM %*% t(dM) - MM
  onesN <- matrix(rep(1, n), ncol = 1)
  recD2 <- 1 / D2
  diag(recD2) <- 0
  MY <- M %*% Y
  MX <- M %*% X
  Poff <- P - diag(diag(P))
  XdM <- dM * X
  YdM <- dM * Y
  XX <- X * X
  XY <- X * Y
  YY <- Y * Y
  onesNMX <- tcrossprod(onesN, MX)
  onesNMY <- tcrossprod(onesN, MY)
  onesNX <- tcrossprod(onesN, X)
  onesNdM <- tcrossprod(onesN, dM)
  MXYY <- MX * YY
  MXXY <- MX * XY
  MYXX <- MY * XX
  MXX <- MX * X
  MYX <- MY * X
  MYY <- MY * Y
  MXY <- MX * Y
  XYdM <- dM * XY
  YYdM <- dM * YY
  XXdM <- dM * XX
  MXXX <- MX * XX
  MonesNX <- M * onesNX

  C0A1vec <- C0A2vec <- C0A3vec <- C0A4vec <- C0A5vec <- rep(0, n)
  C11A1vec <- C11A2vec <- C11A3vec <- C11A4vec <- C11A5vec <- rep(0, n)
  C12A1vec <- C12A2vec <- C12A3vec <- C12A4vec <- C12A5vec <- rep(0, n)

  for (i in 1:n) {
    # Calculation conditioned on i
    Mi <- matrix(M[, i], ncol = 1) # force column vector
    Pi <- matrix(P[, i], ncol = 1) # force column vector
    Pi2 <- Pi^2
    Mi2 <- Mi^2
    D3i <- M[i, i] * D2 - (dM %*% t(Mi)^2 + Mi2 %*% t(dM) - 2 * M * (Mi %*% t(Mi)))
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

    # Calculate things once
    MXi <- MX[i]
    MYi <- MY[i]
    Xi <- X[i]
    XiMii <- Xi * M[i, i]
    XPi <- X * Pi
    XMi <- X * Mi
    XMi2 <- XMi * Mi
    XPiMi <- XPi * Mi
    YPi <- Y * Pi
    YMi <- Y * Mi
    YMi2 <- YMi * Mi
    YPiMi <- YPi * Mi
    Pi2recD2i <- Pi2 * recD2i
    MrecD3i <- M * recD3i
    MMrecD3i <- recD3i * MM
    tXPiD2D3i <- crossprod(XPi, D2D3i)
    D2D3iXMi <- (D2D3i) %*% (XMi)
    Pi2XMi2 <- Pi2 * XMi2
    Pi2XMi <- Pi2 * XMi
    XMiMX <- XMi * MX
    dMXPi <- dM * XPi
    XMi2MY <- XMi2 * MY
    XMi2MX <- XMi2 * MX
    XMiMY <- XMi * MY
    XPiMX <- MX * XPi
    YPiMX <- YPi * MX
    YPidM <- YPi * dM
    XPiMY <- MY * XPi
    D2D3iMonesNX <- D2D3i * MonesNX
    Pi2XdMrecD2i <- Pi2recD2i * XdM
    Pi2YMirecD2i <- Pi2recD2i * YMi
    MiPi2recD2i <- Mi * Pi2recD2i
    Pi2YMi2recD2i <- Pi2recD2i * YMi2
    Pi2YdMrecD2i <- Pi2recD2i * YdM
    onesNdMrecD3i <- onesNdM * recD3i
    onesNMXrecD3i <- onesNMX * recD3i
    Pi2XMirecD2i <- Pi2XMi * recD2i
    Pi2YrecD2i <- Pi2recD2i * Y
    Pi2YdMMirecD2i <- Pi2YdMrecD2i * Mi
    tPi2XdMrecD2iMrecD3i <- crossprod(Pi2XdMrecD2i, MrecD3i)
    tPi2YdMrecD2iMrecD3i <- crossprod(Pi2YdMrecD2i, MrecD3i)
    tXPiMiMrecD3i <- crossprod(XPiMi, MrecD3i)
    tPi2XMi2recD2iMrecD3i <- crossprod(Pi2XMi2 * recD2i, MrecD3i)
    tMXXPiMirecD3i <- crossprod(MX * XPiMi, recD3i)
    tPi2YMirecD2iMMrecD3i <- crossprod(Pi2YMirecD2i, MMrecD3i)
    tPi2YMi2recD2iMrecD3i <- crossprod(Pi2YMi2recD2i, MrecD3i)
    tPi2YdMMirecD2ionesNdMrecD3i <- crossprod(Pi2YdMMirecD2i, onesNdMrecD3i)
    tPi2XMirecD2iMMrecD3i <- crossprod(Pi2XMirecD2i, MMrecD3i)
    tPi2XdMMirecD2irecD3i <- crossprod(Pi2XdMrecD2i * Mi, recD3i)
    tPi2recD2iXMi <- sum(Pi2recD2i * XMi)
    MrecD3iYPiMi <- MrecD3i %*% (YPiMi)
    tXPiD2D3iXPi <- tXPiD2D3i %*% (XPi)
    tPi2recD2iYMi <- sum(Pi2YMirecD2i)
    tPi2recD2iXXdM <- crossprod(Pi2recD2i, XXdM)
    tPi2YdMMirecD2irecD3i <- crossprod(Pi2YdMMirecD2i, recD3i)
    onesNMXrecD3iYPiMi <- (onesNMXrecD3i) %*% (YPiMi)
    D2D3iMonesNXMi <- (D2D3iMonesNX) %*% (Mi)
    MrecD3iXPiMi <- t(tXPiMiMrecD3i)
    Pi2recD2iMYX <- Pi2recD2i * MYX

    C0A11i <- tXPiD2D3iXPi * MYi
    C0A12i <- crossprod(MY * XPiMi, recD3i) %*% (dMXPi) -
      tXPiMiMrecD3i %*% (XPiMY)
    C0A13i <- crossprod(dMXPi, onesNMY * recD3i) %*% (XPiMi)
    C0A14i <- sum(XPiMY * MrecD3iXPiMi)
    C0A15i <- MYi * (tPi2recD2iXXdM) -
      sum(MiPi2recD2i * MYXX)

    C0A21i <- tXPiD2D3i %*% (YPi) * (MXi)
    C0A22i <- tMXXPiMirecD3i %*% (YPidM) -
      tXPiMiMrecD3i %*% (YPiMX)
    C0A23i <- crossprod(dMXPi, onesNMXrecD3iYPiMi)
    C0A24i <- crossprod(MX * XPi, MrecD3iYPiMi)
    C0A25i <- (MXi) * sum(Pi2recD2i * XYdM) -
      sum(MiPi2recD2i * MXXY)

    C0A31i <- crossprod(YPi, D2D3i) %*% (YPi) * (MXi)
    C0A32i <- crossprod(MX * YPiMi, recD3i) %*% (YPidM) -
      crossprod(YPiMi, MrecD3i) %*% (YPiMX)
    C0A33i <- crossprod(YPidM, onesNMXrecD3iYPiMi)
    C0A34i <- crossprod(YPiMX, MrecD3iYPiMi)
    C0A35i <- (MXi) * sum(Pi2recD2i * YYdM) -
      sum(MiPi2recD2i * MXYY)

    C0A41i <- crossprod(Pi2YdMrecD2i * MY, D2D3iXMi) -
      crossprod(Pi2YrecD2i * MY, D2D3iMonesNXMi)
    C0A42i <- tPi2YdMrecD2iMrecD3i %*% (XMi2) -
      tPi2YMirecD2iMMrecD3i %*% (XMi) -
      tPi2YdMMirecD2ionesNdMrecD3i %*% (XMi) +
      tPi2YMi2recD2iMrecD3i %*% (XdM)
    C0A43i <- tPi2YdMMirecD2irecD3i %*% (XMi2MY) -
      tPi2YMi2recD2iMrecD3i %*% (XMiMY) -
      M[i, i] * tPi2YdMrecD2iMrecD3i %*% (XMiMY) +
      M[i, i] * tPi2YMirecD2iMMrecD3i %*% (MYX)
    C0A44i <- XiMii * sum(Pi2recD2i * MYY) -
      Xi * (MYi) * (tPi2recD2iYMi)

    C0A51i <- crossprod(Pi2YdMrecD2i * MX, D2D3iXMi) -
      crossprod(Pi2YrecD2i * MX, D2D3iMonesNXMi)
    C0A52i <- C0A42i
    C0A53i <- tPi2YdMMirecD2irecD3i %*% (XMi2MX) -
      tPi2YMi2recD2iMrecD3i %*% (XMiMX) -
      M[i, i] * tPi2YdMrecD2iMrecD3i %*% (XMiMX) +
      M[i, i] * tPi2YMirecD2iMMrecD3i %*% (MXX)
    C0A54i <- XiMii * sum(Pi2recD2i * MXY) -
      Xi * (MXi) * (tPi2recD2iYMi)

    C11A41i <- crossprod(Pi2XdMrecD2i * MY, D2D3iXMi) -
      crossprod(Pi2recD2iMYX, D2D3iMonesNXMi)
    C11A42i <- tPi2XdMrecD2iMrecD3i %*% (XMi2) -
      tPi2XMirecD2iMMrecD3i %*% (XMi) -
      crossprod(Pi2 * XdM * Mi * recD2i, onesNdMrecD3i) %*% (XMi) +
      tPi2XMi2recD2iMrecD3i %*% (XdM)
    C11A43i <- tPi2XdMMirecD2irecD3i %*% (XMi2MY) -
      tPi2XMi2recD2iMrecD3i %*% (XMiMY) -
      M[i, i] * tPi2XdMrecD2iMrecD3i %*% (XMiMY) +
      M[i, i] * tPi2XMirecD2iMMrecD3i %*% (MYX)
    C11A44i <- XiMii * sum(Pi2recD2iMYX) -
      Xi * (MYi) * (tPi2recD2iXMi)

    C12A11i <- tXPiD2D3iXPi * (MXi)
    C12A12i <- tMXXPiMirecD3i %*% (dMXPi) -
      tXPiMiMrecD3i %*% (XPiMX)
    C12A13i <- crossprod(dMXPi, onesNMX * recD3i) %*% (XPiMi)
    C12A14i <- crossprod(XPiMX, MrecD3iXPiMi)
    C12A15i <- (MXi) * (tPi2recD2iXXdM) -
      sum(MiPi2recD2i * MXXX)

    C12A51i <- crossprod(Pi2XdMrecD2i * MX, D2D3iXMi) -
      crossprod(Pi2 * X * recD2i * MX, D2D3iMonesNXMi)
    C12A52i <- C11A42i
    C12A53i <- tPi2XdMMirecD2irecD3i %*% (XMi2MX) -
      tPi2XMi2recD2iMrecD3i %*% (XMiMX) -
      M[i, i] * tPi2XdMrecD2iMrecD3i %*% (XMiMX) +
      M[i, i] * tPi2XMirecD2iMMrecD3i %*% (MXX)
    C12A54i <- XiMii * sum(Pi2recD2i * MXX) -
      Xi * (MXi) * (tPi2recD2iXMi)

    C0A1vec[i] <- C0A11i - C0A12i - C0A13i + C0A14i + C0A15i
    C0A2vec[i] <- C0A21i - C0A22i - C0A23i + C0A24i + C0A25i
    C0A3vec[i] <- C0A31i - C0A32i - C0A33i + C0A34i + C0A35i
    C0A4vec[i] <- C0A41i + C0A42i * MYi + C0A43i + C0A44i
    C0A5vec[i] <- C0A51i + C0A52i * MXi + C0A53i + C0A54i

    C11A4vec[i] <- C11A41i + C11A42i * MYi + C11A43i + C11A44i
    C12A1vec[i] <- C12A11i - C12A12i - C12A13i + C12A14i + C12A15i
    C12A5vec[i] <- C12A51i + C12A52i * MXi + C12A53i + C12A54i

    if (noisy) {
      cat(i, "of", n, "done. ")
    }
  }

  C11A1vec <- C0A1vec
  C11A3vec <- C11A2vec <- C0A2vec
  C11A5vec <- C0A5vec
  C12A2vec <- C12A1vec
  C12A3vec <- C0A2vec
  C12A4vec <- C0A5vec
  C2A1vec <- C2A2vec <- C2A3vec <- C12A1vec
  C2A4vec <- C2A5vec <- C12A5vec

  C0 <- sum(C0A1vec * Y + 2 * C0A2vec * Y + C0A3vec * X - C0A4vec * X - C0A5vec * Y)
  C11 <- sum(C11A1vec * X + 2 * C11A2vec * X + C11A3vec * X - C11A4vec * X - C11A5vec * X)
  C12 <- sum(C12A1vec * Y + 2 * C12A2vec * Y + C12A3vec * X - C12A4vec * X - C12A5vec * Y)
  C1 <- -C11 - C12
  C2 <- sum(C2A1vec * X + 2 * C2A2vec * X + C2A3vec * X - C2A4vec * X - C12A5vec * X)

  PXY <- t(Y) %*% Poff %*% X
  PXX <- t(X) %*% Poff %*% X


  acon <- PXX^2 - q * C2
  bcon <- -2 * PXY * PXX - q * C1
  ccon <- PXY^2 - q * C0

  c(acon, bcon, ccon)
}

#' Compute Quadratic Roots for Confidence Sets
#'
#' @description
#' Computes the roots of the quadratic equation \eqn{a\beta^2 + b\beta + c = 0}.
#' These roots serve as the boundaries for the confidence sets constructed by inverting
#' the UJIVE/LIML score test.
#'
#' @param CIcoef Numeric vector of length 3. The coefficients \eqn{(a, b, c)} of the quadratic
#'   inequality, typically obtained from \code{\link{GetCIcoef}}.
#'
#' @details
#' This function applies the standard quadratic formula:
#' \deqn{\beta_{1,2} = \frac{-b \pm \sqrt{b^2 - 4ac}}{2a}}
#'
#'
#' \strong{Note:} This function does not check the sign of the discriminant. It is intended
#' to be called by a wrapper function (like \code{\link{GetCItypebd}}) that first verifies
#' the existence of real roots (i.e., \eqn{b^2 - 4ac \ge 0}).
#'
#' Depending on the sign of \eqn{a} (convexity), these values represent either the endpoints
#' of a bounded confidence interval or the inner boundaries of a disjoint ("donut") confidence set.
#'
#' @return Numeric vector of length 2. Contains the two roots.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @export
GetCIvals <- function(CIcoef) {
  acon <- CIcoef[1]
  bcon <- CIcoef[2]
  ccon <- CIcoef[3]
  det <- CIcoef[2]^2 - 4 * CIcoef[1] * CIcoef[3]
  CILB <- (-bcon - sqrt(det)) / (2 * acon)
  CIUB <- (-bcon + sqrt(det)) / (2 * acon)
  c(CILB, CIUB)
}

#' Classify and Compute Confidence Interval Bounds
#'
#' @description
#' Solves the quadratic inequality \eqn{a\beta^2 + b\beta + c \leq 0} derived from the
#' score test inversion to determine the topology and boundaries of the confidence set
#' for \eqn{\beta}.
#'
#' @param CIcoef Numeric vector of length 3. The coefficients \eqn{(a, b, c)} obtained
#'   from \code{\link{GetCIcoef}}.
#'
#' @details
#' The confidence set is defined as \eqn{\{ \beta : a\beta^2 + b\beta + c \leq 0 \}}.
#' Depending on the sign of \eqn{a} and the discriminant \eqn{\Delta = b^2 - 4ac},
#' this set can take one of four forms:
#'
#'
#' \itemize{
#'   \item \strong{Type 1: Bounded Interval} (\eqn{a \ge 0, \Delta \ge 0}).
#'   The parabola opens upward with real roots. The CI is the closed interval \eqn{[\beta_{min}, \beta_{max}]}.
#'
#'   \item \strong{Type 2: Disjoint Union (Donut)} (\eqn{a < 0, \Delta \ge 0}).
#'   The parabola opens downward with real roots. The CI is the union of two infinite rays:
#'   \eqn{(-\infty, \beta_{min}] \cup [\beta_{max}, \infty)}.
#'
#'   \item \strong{Type 3: Real Line} (\eqn{a < 0, \Delta < 0}).
#'   The parabola is always negative. The CI includes the entire real line.
#'   Bounds are returned as placeholders \code{c(-100, 100)}.
#'
#'   \item \strong{Type 4: Empty Set} (\eqn{a \ge 0, \Delta < 0}).
#'   The parabola is always positive. The confidence set is empty, implying the model
#'   is rejected at the specified significance level for all \eqn{\beta}.
#' }
#'
#' @return Numeric vector of length 3. Format: \code{c(CItype, LowerBound, UpperBound)}.
#'
#' @export
GetCItypebd <- function(CIcoef) {
  det <- CIcoef[2]^2 - 4 * CIcoef[1] * CIcoef[3]
  if (CIcoef[1] >= 0 & det >= 0) {
    CItype <- 1 # convex interval
    CIbounds <- GetCIvals(CIcoef)
  } else if (CIcoef[1] < 0 & det >= 0) {
    CItype <- 2 # donut
    CIbounds <- GetCIvals(CIcoef)
  } else if (CIcoef[1] < 0 & det < 0) {
    CItype <- 3 # accept everything
    CIbounds <- c(-100, 100)
  } else {
    CItype <- 4 # reject everything
    CIbounds <- c(NA, NA)
  }
  c(CItype, CIbounds)
}

#' Compute Quadratic Coefficients for CI (Stratified Asymmetric)
#'
#' @description
#' Calculates the coefficients \eqn{(a, b, c)} for the confidence interval quadratic inequality
#' in a general design that includes both grouping (stratification) and covariates.
#' This function is optimized for designs where projection matrices are block-diagonal
#' with respect to \code{groupW} but potentially asymmetric within blocks due to covariate adjustments.
#'
#' @param df Data frame. Contains the observable variables and grouping indicators.
#' @param groupW Column name (unquoted). The covariate stratification variable.
#' @param groupQ Column name (unquoted). The instrument grouping variable.
#' @param X Column name (unquoted). The endogenous regressor.
#' @param Y Column name (unquoted). The outcome variable.
#' @param MX Column name (unquoted). Leverage-adjusted regressor (\eqn{M X}).
#' @param MY Column name (unquoted). Leverage-adjusted outcome (\eqn{M Y}).
#' @param q Numeric scalar. The critical value for test inversion (typically \eqn{1.96^2}).
#'   Defaults to \code{qnorm(.975)^2}.
#' @param noisy Logical. If \code{TRUE}, prints progress.
#'   Defaults to \code{FALSE}.
#'
#' @details
#' This function performs a loop:
#' \enumerate{
#'   \item Splits data by \code{groupW}. Computes local projection matrices
#'   \eqn{P_Q} and \eqn{P_W} for each stratum.
#'   \item Calculates the full UJIVE weighting matrix \eqn{G = U(P_Q) - U(P_W)} locally.
#'   Since \eqn{G} is asymmetric, it computes all 5 variance components.
#'   \item Accumulates all polynomial coefficients for \eqn{\hat{V}(\beta)}
#'   simultaneously using extensive pre-calculation of vector products.
#' }
#'
#' @return Numeric vector of length 3 containing \code{c(a, b, c)}.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @export
GetL3OCIcoef_fast <- function(df, groupW, groupQ, X, Y, MX, MY, q = qnorm(.975)^2, noisy = FALSE) {
  df$X <- eval(substitute(X), df)
  df$Y <- eval(substitute(Y), df)
  df$MX <- eval(substitute(MX), df)
  df$MY <- eval(substitute(MY), df)
  df$groupW <- eval(substitute(groupW), df)
  df$groupQ <- eval(substitute(groupQ), df)

  C0A11vecs <- C0A12vecs <- C0A13vecs <- C0A14vecs <- C0A15vecs <- rep(0, max(df$groupQ))
  C0A21vecs <- C0A22vecs <- C0A23vecs <- C0A24vecs <- C0A25vecs <- rep(0, max(df$groupQ))
  C0A31vecs <- C0A32vecs <- C0A33vecs <- C0A34vecs <- C0A35vecs <- rep(0, max(df$groupQ))
  C0A41vecs <- C0A42vecs <- C0A43vecs <- C0A44vecs <- rep(0, max(df$groupQ))
  C0A51vecs <- C0A52vecs <- C0A53vecs <- C0A54vecs <- rep(0, max(df$groupQ))
  C11A11vecs <- C11A12vecs <- C11A13vecs <- C11A14vecs <- C11A15vecs <- rep(0, max(df$groupQ))
  C11A21vecs <- C11A22vecs <- C11A23vecs <- C11A24vecs <- C11A25vecs <- rep(0, max(df$groupQ))
  C11A41vecs <- C11A42vecs <- C11A43vecs <- C11A44vecs <- rep(0, max(df$groupQ))
  C11A51vecs <- C11A52vecs <- C11A53vecs <- C11A54vecs <- rep(0, max(df$groupQ))
  C12A11vecs <- C12A12vecs <- C12A13vecs <- C12A14vecs <- C12A15vecs <- rep(0, max(df$groupQ))
  C12A31vecs <- C12A32vecs <- C12A33vecs <- C12A34vecs <- C12A35vecs <- rep(0, max(df$groupQ))
  C12A51vecs <- C12A52vecs <- C12A53vecs <- C12A54vecs <- rep(0, max(df$groupQ))
  C2A11vecs <- C2A12vecs <- C2A13vecs <- C2A14vecs <- C2A15vecs <- rep(0, max(df$groupQ))
  C2A41vecs <- C2A42vecs <- C2A43vecs <- C2A44vecs <- rep(0, max(df$groupQ))

  iteration <- 1

  # outer loop
  for (s in unique(df$groupW)) {
    ds <- df[df$groupW == s, ]

    # ds <- ds %>%  group_by(groupQ) %>% mutate(numingrp = length(ds$groupQ))
    # ds <- ds[ds$numingrp>=3,]

    ZWmat <- matrix(1, nrow = length(ds$groupW), ncol = length(ds$groupW))
    PW <- ZWmat / (length(ds$groupW))
    if (nrow(ds) == 0) {
      for (g in unique(ds$groupQ)) {
        C0A11vecs[g] <- C0A12vecs[g] <- C0A13vecs[g] <- C0A14vecs[g] <- C0A15vecs[g] <- 0
        C0A21vecs[g] <- C0A22vecs[g] <- C0A23vecs[g] <- C0A24vecs[g] <- C0A25vecs[g] <- 0
        C0A31vecs[g] <- C0A32vecs[g] <- C0A33vecs[g] <- C0A34vecs[g] <- C0A35vecs[g] <- 0
        C0A41vecs[g] <- C0A42vecs[g] <- C0A43vecs[g] <- C0A44vecs[g] <- 0
        C0A51vecs[g] <- C0A52vecs[g] <- C0A53vecs[g] <- C0A54vecs[g] <- 0
        C11A11vecs[g] <- C11A12vecs[g] <- C11A13vecs[g] <- C11A14vecs[g] <- C11A15vecs[g] <- 0
        C11A21vecs[g] <- C11A22vecs[g] <- C11A23vecs[g] <- C11A24vecs[g] <- C11A25vecs[g] <- 0
        C11A41vecs[g] <- C11A42vecs[g] <- C11A43vecs[g] <- C11A44vecs[g] <- 0
        C11A51vecs[g] <- C11A52vecs[g] <- C11A53vecs[g] <- C11A54vecs[g] <- 0
        C12A11vecs[g] <- C12A12vecs[g] <- C12A13vecs[g] <- C12A14vecs[g] <- C12A15vecs[g] <- 0
        C12A31vecs[g] <- C12A32vecs[g] <- C12A33vecs[g] <- C12A34vecs[g] <- C12A35vecs[g] <- 0
        C12A51vecs[g] <- C12A52vecs[g] <- C12A53vecs[g] <- C12A54vecs[g] <- 0
        C2A11vecs[g] <- C2A12vecs[g] <- C2A13vecs[g] <- C2A14vecs[g] <- C2A15vecs[g] <- 0
        C2A41vecs[g] <- C2A42vecs[g] <- C2A43vecs[g] <- C2A44vecs[g] <- 0
      }
    } else {
      ZQ <- matrix(0, nrow = length(ds$groupQ), ncol = length(unique(ds$groupQ)))
      ds$groupidx <- Getgroupindex(ds, groupQ)
      ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1

      PQ <- ZQ %*% solve(t(ZQ) %*% ZQ) %*% t(ZQ)

      # calculate values specific to this subset
      Gs <- diag(1 / (diag(diag(nrow(ds))) - pmin(diag(PQ), .99))) %*% (PQ - diag(diag(PQ))) -
        diag(1 / (diag(diag(nrow(ds))) - pmin(diag(PW), .99))) %*% (PW - diag(diag(PW)))
      Ps <- PQ
      Ms <- diag(nrow(ds)) - Ps
      dMs <- matrix(diag(Ms), ncol = 1)
      D2s <- dMs %*% t(dMs) - Ms * Ms
      recD2s <- 1 / D2s
      diag(recD2s) <- 0

      # Stuff to calculate once
      XMY <- ds$X * ds$MY
      XMX <- ds$X * ds$MX
      YMX <- ds$Y * ds$MX
      YMY <- ds$Y * ds$MY
      XdMs <- ds$X * dMs
      YdMs <- ds$Y * dMs
      XMYX <- XMY * ds$X
      ones <- matrix(rep(1, nrow(ds)), ncol = 1)
      onestdMs <- ones %*% t(dMs)
      onestMX <- (ones %*% t(ds$MX))
      onestMY <- (ones %*% t(ds$MY))
      onesX <- ones %x% t(ds$X)

      # inner loop
      for (g in unique(ds$groupQ)) {
        repidx <- min(which(ds$groupQ == g)) # representative index
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
        # Mes <- matrix(ds$lpos,ncol=1)

        ## Things to calculate only once
        XGis <- ds$X * Gis
        XMis <- ds$X * Mis
        YGis <- ds$Y * Gis
        YMis <- ds$Y * Mis
        YMis2 <- YMis * Mis
        XMis2 <- XMis * Mis
        Gs2 <- Gs^2
        Gis2 <- Gis^2
        Mis2 <- Mis^2
        Ms2 <- Ms^2
        GisMis <- Gis * Mis
        YGisMis <- YGis * Mis
        XGisMis <- XGis * Mis
        Mis2X <- Mis2 * ds$X
        XGisdMs <- XGis * dMs
        YGisdMs <- YGis * dMs
        XMYGis <- XMY * Gis
        XMYGisMis <- XMYGis * Mis
        Gis2YdMs <- Gis2 * YdMs
        Gis2XdMs <- Gis2 * XdMs
        recD2isMX <- ds$MX * recD2is
        YGisMX <- YGis * ds$MX
        Mis2XMX <- Mis2X * ds$MX
        XGisMX <- XGis * ds$MX
        MisXMX <- Mis * XMX
        Gis2XdMsMisrecD2is <- Gis2XdMs * Mis * recD2is
        Gis2YdMsMisrecD2is <- Gis2YdMs * Mis * recD2is
        ivectomatsY <- ivectomats(ds, ds$Y, g)
        ivectomatsX <- ivectomats(ds, ds$X, g)
        ivectomatsXdMs <- ivectomats(ds, XdMs, g)
        ivectomatsYMX <- ivectomats(ds, YMX, g)
        ivectomatsXMY <- ivectomats(ds, XMY, g)
        ivectomatsXMX <- ivectomats(ds, XMX, g)
        D2D3isivectomatsX <- D2D3is * ivectomatsX
        D2D3isivectomatsXMY <- D2D3is * ivectomatsXMY
        D2D3isivectomatsXMX <- D2D3is * ivectomatsXMX
        recD3isivectomatsY <- recD3is * ivectomatsY
        recD3isMsivectomatsY <- recD3isivectomatsY * Ms
        recD3isivectomatsX <- recD3is * ivectomatsX
        recD3isMsivectomatsX <- recD3isivectomatsX * Ms
        recD3isonestMXivectomatsX <- recD3is * onestMX * ivectomatsX
        recD2sGs2 <- recD2s * Gs2
        recD3isMs <- recD3is * Ms
        recD3isMs2 <- recD3isMs * Ms
        recD3isMsivectomatsXdMs <- recD3isMs * ivectomatsXdMs
        recD3isMsivectomatsXMY <- recD3isMs * ivectomatsXMY
        recD3isMs2ivectomatsXMY <- recD3isMs2 * ivectomatsXMY
        recD3isMsivectomatsYMX <- recD3isMs * ivectomatsYMX
        Gs2MsrecD2sPgs <- Gs2 * Ms * recD2s * Pgs
        recD2sGs2Pgs <- recD2sGs2 * Pgs
        D2D3isMsonesX <- D2D3is * Ms * (onesX)
        D2D3isMsonesXivectomatsX <- D2D3isMsonesX * ivectomatsX
        recD3isonestdMs <- recD3is * (onestdMs)
        recD3isonestdMsivectomatsXMY <- recD3isonestdMs * ivectomatsXMY
        recD3isonestdMsivectomatsXMX <- recD3isonestdMs * ivectomatsXMX

        # Calculations for every groupQ
        C0A11vecs[g] <- t(XGis) %*% (D2D3is * ivectomats(ds, YMY, g)) %*% (XGis)
        C0A12vecs[g] <- t(XMYGisMis) %*% (recD3isivectomatsY) %*% (XGisdMs) -
          t(XGisMis) %*% (recD3isMsivectomatsY) %*% (Gis * XMY)
        C0A13vecs[g] <- t(XdMs * Gis) %*% (recD3isivectomatsY * onestMY) %*% (XGisMis)
        C0A14vecs[g] <- t(XMYGis) %*% (recD3isMsivectomatsY) %*% (XGisMis)
        C0A15vecs[g] <- t(YMY) %*% (recD2sGs2Pgs) %*% (ds$X * XdMs) -
          t(ds$Y) %*% (Gs2MsrecD2sPgs) %*% (XMYX)

        C0A21vecs[g] <- t(XGis) %*% (D2D3is * ivectomatsYMX) %*% (YGis)
        C0A22vecs[g] <- t(XMX * GisMis) %*% (recD3isivectomatsY) %*% (YGisdMs) -
          t(XGisMis) %*% (recD3isMsivectomatsY) %*% (YGisMX)
        C0A23vecs[g] <- t(XdMs * Gis) %*% (recD3isivectomatsY * onestMX) %*% (YGisMis)
        C0A24vecs[g] <- t(XMX * Gis) %*% (recD3isMsivectomatsY) %*% (YGisMis)
        C0A25vecs[g] <- t(YMX) %*% (recD2sGs2Pgs) %*% (ds$X * YdMs) -
          t(ds$Y) %*% (Gs2MsrecD2sPgs) %*% (XMX * ds$Y)

        C0A31vecs[g] <- t(YGis) %*% (D2D3isivectomatsXMX) %*% (YGis)
        C0A32vecs[g] <- t(YMX * GisMis) %*% (recD3isivectomatsX) %*% (YGisdMs) -
          t(YGisMis) %*% (recD3isMsivectomatsX) %*% (YGisMX)
        C0A33vecs[g] <- t(YdMs * Gis) %*% (recD3isonestMXivectomatsX) %*% (YGisMis)
        C0A34vecs[g] <- t(YMX * Gis) %*% (recD3isMsivectomatsX) %*% (YGisMis)
        C0A35vecs[g] <- t(XMX) %*% (recD2sGs2Pgs) %*% (ds$Y * YdMs) -
          t(ds$X) %*% (Gs2MsrecD2sPgs) %*% (YMX * ds$Y)

        C0A41vecs[g] <- t(Gis2YdMs * ds$MY * recD2is) %*% (D2D3isivectomatsX) %*% (XMis) -
          t(Gis2 * YMY * recD2is) %*% (D2D3isMsonesXivectomatsX) %*% (Mis)
        C0A42vecs[g] <- t(Gis2YdMs * recD2is) %*% (recD3isMsivectomatsXMY) %*% (Mis2X) -
          t(Gis2 * YMis * recD2is) %*% (recD3isMs2ivectomatsXMY) %*% (XMis) -
          t(Gis2YdMsMisrecD2is) %*% (recD3isonestdMsivectomatsXMY) %*% (XMis) +
          t(Gis2 * YMis2 * recD2is) %*% (recD3isMsivectomatsXMY) %*% (dMs * ds$X)
        C0A43vecs[g] <- t(Gis2YdMsMisrecD2is) %*% (recD3isivectomatsX) %*% (Mis2 * XMY) -
          t(Gis2 * YMis2 * recD2is) %*% (recD3isMsivectomatsX) %*% (Mis * XMY) -
          t(Gis2YdMs * recD2is) %*% (recD3isMsivectomatsXdMs) %*% (Mis * XMY) +
          t(Gis2 * YMis * recD2is) %*% (recD3isMs2 * ivectomatsXdMs) %*% (XMY)
        C0A44vecs[g] <- t(ds$X * XdMs) %*% (recD2sGs2Pgs) %*% (YMY) -
          t(ds$X * XMY) %*% (recD2sGs2Pgs) %*% (YMis)

        C0A51vecs[g] <- t(Gis2YdMs * recD2isMX) %*% (D2D3is * ivectomatsY) %*% (XMis) -
          t(Gis2 * YMX * recD2is) %*% (D2D3isMsonesX * ivectomatsY) %*% (Mis)
        C0A52vecs[g] <- t(Gis2YdMs * recD2is) %*% (recD3isMsivectomatsYMX) %*% (Mis2X) -
          t(Gis2 * YMis * recD2is) %*% (recD3isMs2 * ivectomatsYMX) %*% (XMis) -
          t(Gis2YdMsMisrecD2is) %*% (recD3isonestdMs * ivectomatsYMX) %*% (XMis) +
          t(Gis2 * YMis2 * recD2is) %*% (recD3isMsivectomatsYMX) %*% (dMs * ds$X)
        C0A53vecs[g] <- t(Gis2YdMsMisrecD2is) %*% (recD3isivectomatsY) %*% (Mis2XMX) -
          t(Gis2 * YMis2 * recD2is) %*% (recD3isMsivectomatsY) %*% (MisXMX) -
          t(Gis2YdMs * recD2is) %*% (recD3isMs * ivectomats(ds, YdMs, g)) %*% (MisXMX) +
          t(Gis2 * YMis * recD2is) %*% (recD3isMs2 * ivectomats(ds, YdMs, g)) %*% (XMX)
        C0A54vecs[g] <- t(ds$Y * XdMs) %*% (recD2sGs2Pgs) %*% (YMX) -
          t(ds$Y * XMX) %*% (recD2sGs2Pgs) %*% (YMis)

        C11A11vecs[g] <- t(XGis) %*% (D2D3isivectomatsXMY) %*% (XGis)
        C11A12vecs[g] <- t(XMYGisMis) %*% (recD3isivectomatsX) %*% (XGisdMs) -
          t(XGisMis) %*% (recD3isMsivectomatsX) %*% (Gis * XMY)
        C11A13vecs[g] <- t(XdMs * Gis) %*% (recD3is * onestMY * ivectomatsX) %*% (XGisMis)
        C11A14vecs[g] <- t(XMYGis) %*% (recD3isMsivectomatsX) %*% (XGisMis)
        C11A15vecs[g] <- t(XMY) %*% (recD2sGs2Pgs) %*% (ds$X * XdMs) -
          t(ds$X) %*% (Gs2MsrecD2sPgs) %*% (XMYX)

        C11A21vecs[g] <- t(XGis) %*% (D2D3isivectomatsXMX) %*% (YGis)
        C11A22vecs[g] <- t(XMX * GisMis) %*% (recD3isivectomatsX) %*% (YGisdMs) -
          t(XGisMis) %*% (recD3isMsivectomatsX) %*% (YGisMX)
        C11A23vecs[g] <- t(XdMs * Gis) %*% (recD3isonestMXivectomatsX) %*% (YGisMis)
        C11A24vecs[g] <- t(XMX * Gis) %*% (recD3isMsivectomatsX) %*% (YGisMis)
        C11A25vecs[g] <- t(XMX) %*% (recD2sGs2Pgs) %*% (ds$X * YdMs) -
          t(ds$X) %*% (Gs2MsrecD2sPgs) %*% (XMX * ds$Y)

        C11A41vecs[g] <- t(Gis2XdMs * ds$MY * recD2is) %*% (D2D3isivectomatsX) %*% (XMis) -
          t(Gis2 * XMY * recD2is) %*% (D2D3isMsonesXivectomatsX) %*% (Mis)
        C11A42vecs[g] <- t(Gis2XdMs * recD2is) %*% (recD3isMsivectomatsXMY) %*% (Mis2X) -
          t(Gis2 * XMis * recD2is) %*% (recD3isMs2ivectomatsXMY) %*% (XMis) -
          t(Gis2XdMsMisrecD2is) %*% (recD3isonestdMsivectomatsXMY) %*% (XMis) +
          t(Gis2 * XMis2 * recD2is) %*% (recD3isMsivectomatsXMY) %*% (dMs * ds$X)
        C11A43vecs[g] <- t(Gis2XdMsMisrecD2is) %*% (recD3isivectomatsX) %*% (Mis2 * XMY) -
          t(Gis2 * XMis2 * recD2is) %*% (recD3isMsivectomatsX) %*% (Mis * XMY) -
          t(Gis2XdMs * recD2is) %*% (recD3isMsivectomatsXdMs) %*% (Mis * XMY) +
          t(Gis2 * XMis * recD2is) %*% (recD3isMs2 * ivectomatsXdMs) %*% (XMY)
        C11A44vecs[g] <- t(ds$X * XdMs) %*% (recD2sGs2Pgs) %*% (XMY) -
          t(ds$X * XMY) %*% (recD2sGs2Pgs) %*% (XMis)

        C11A51vecs[g] <- t(Gis2YdMs * recD2isMX) %*% (D2D3isivectomatsX) %*% (XMis) -
          t(Gis2 * YMX * recD2is) %*% (D2D3isMsonesXivectomatsX) %*% (Mis)
        C11A52vecs[g] <- t(Gis2YdMs * recD2is) %*% (recD3isMs * ivectomatsXMX) %*% (Mis2X) -
          t(Gis2 * YMis * recD2is) %*% (recD3isMs2 * ivectomatsXMX) %*% (XMis) -
          t(Gis2YdMsMisrecD2is) %*% (recD3isonestdMsivectomatsXMX) %*% (XMis) +
          t(Gis2 * YMis2 * recD2is) %*% (recD3isMs * ivectomatsXMX) %*% (dMs * ds$X)
        C11A53vecs[g] <- t(Gis2YdMsMisrecD2is) %*% (recD3isivectomatsX) %*% (Mis2XMX) -
          t(Gis2 * YMis2 * recD2is) %*% (recD3isMsivectomatsX) %*% (MisXMX) -
          t(Gis2YdMs * recD2is) %*% (recD3isMsivectomatsXdMs) %*% (MisXMX) +
          t(Gis2 * YMis * recD2is) %*% (recD3isMs2 * ivectomatsXdMs) %*% (XMX)
        C11A54vecs[g] <- t(ds$X * XdMs) %*% (recD2sGs2Pgs) %*% (YMX) -
          t(ds$X * XMX) %*% (recD2sGs2Pgs) %*% (YMis)

        C12A11vecs[g] <- t(XGis) %*% (D2D3is * ivectomatsYMX) %*% (XGis)
        C12A12vecs[g] <- t(XMX * GisMis) %*% (recD3isivectomatsY) %*% (XGisdMs) -
          t(XGisMis) %*% (recD3isMsivectomatsY) %*% (XGisMX)
        C12A13vecs[g] <- t(XdMs * Gis) %*% (recD3is * onestMX * ivectomatsY) %*% (XGisMis)
        C12A14vecs[g] <- t(XMX * Gis) %*% (recD3isMsivectomatsY) %*% (XGisMis)
        C12A15vecs[g] <- t(YMX) %*% (recD2sGs2Pgs) %*% (ds$X * XdMs) -
          t(ds$Y) %*% (Gs2MsrecD2sPgs) %*% (XMX * ds$X)

        C12A31vecs[g] <- t(YGis) %*% (D2D3isivectomatsXMX) %*% (XGis)
        C12A32vecs[g] <- t(YMX * GisMis) %*% (recD3isivectomatsX) %*% (XGisdMs) -
          t(YGisMis) %*% (recD3isMsivectomatsX) %*% (XGisMX)
        C12A33vecs[g] <- t(YdMs * Gis) %*% (recD3isonestMXivectomatsX) %*% (XGisMis)
        C12A34vecs[g] <- t(YMX * Gis) %*% (recD3isMsivectomatsX) %*% (XGisMis)
        C12A35vecs[g] <- t(XMX) %*% (recD2sGs2Pgs) %*% (ds$Y * XdMs) -
          t(ds$X) %*% (Gs2MsrecD2sPgs) %*% (YMX * ds$X)

        C12A51vecs[g] <- t(Gis2XdMs * recD2isMX) %*% (D2D3is * ivectomatsY) %*% (XMis) -
          t(Gis2 * XMX * recD2is) %*% (D2D3isMsonesX * ivectomatsY) %*% (Mis)
        C12A52vecs[g] <- t(Gis2XdMs * recD2is) %*% (recD3isMsivectomatsYMX) %*% (Mis2X) -
          t(Gis2 * XMis * recD2is) %*% (recD3isMs2 * ivectomatsYMX) %*% (XMis) -
          t(Gis2XdMsMisrecD2is) %*% (recD3isonestdMs * ivectomatsYMX) %*% (XMis) +
          t(Gis2 * XMis2 * recD2is) %*% (recD3isMsivectomatsYMX) %*% (dMs * ds$X)
        C12A53vecs[g] <- t(Gis2XdMsMisrecD2is) %*% (recD3isivectomatsY) %*% (Mis2XMX) -
          t(Gis2 * XMis2 * recD2is) %*% (recD3isMsivectomatsY) %*% (MisXMX) -
          t(Gis2XdMs * recD2is) %*% (recD3isMs * ivectomats(ds, YdMs, g)) %*% (MisXMX) +
          t(Gis2 * XMis * recD2is) %*% (recD3isMs2 * ivectomats(ds, YdMs, g)) %*% (XMX)
        C12A54vecs[g] <- t(ds$Y * XdMs) %*% (recD2sGs2Pgs) %*% (XMX) -
          t(ds$Y * XMX) %*% (recD2sGs2Pgs) %*% (XMis)

        C2A11vecs[g] <- t(XGis) %*% (D2D3isivectomatsXMX) %*% (XGis)
        C2A12vecs[g] <- t(XMX * GisMis) %*% (recD3isivectomatsX) %*% (XGisdMs) -
          t(XGisMis) %*% (recD3isMsivectomatsX) %*% (XGisMX)
        C2A13vecs[g] <- t(XdMs * Gis) %*% (recD3isonestMXivectomatsX) %*% (XGisMis)
        C2A14vecs[g] <- t(XMX * Gis) %*% (recD3isMsivectomatsX) %*% (XGisMis)
        C2A15vecs[g] <- t(XMX) %*% (recD2sGs2Pgs) %*% (ds$X * XdMs) -
          t(ds$X) %*% (Gs2MsrecD2sPgs) %*% (XMX * ds$X)

        C2A41vecs[g] <- t(Gis2XdMs * recD2isMX) %*% (D2D3isivectomatsX) %*% (XMis) -
          t(Gis2 * XMX * recD2is) %*% (D2D3isMsonesXivectomatsX) %*% (Mis)
        C2A42vecs[g] <- t(Gis2XdMs * recD2is) %*% (recD3isMs * ivectomatsXMX) %*% (Mis2X) -
          t(Gis2 * XMis * recD2is) %*% (recD3isMs2 * ivectomatsXMX) %*% (XMis) -
          t(Gis2XdMsMisrecD2is) %*% (recD3isonestdMsivectomatsXMX) %*% (XMis) +
          t(Gis2 * XMis2 * recD2is) %*% (recD3isMs * ivectomatsXMX) %*% (dMs * ds$X)
        C2A43vecs[g] <- t(Gis2XdMsMisrecD2is) %*% (recD3isivectomatsX) %*% (Mis2XMX) -
          t(Gis2 * XMis2 * recD2is) %*% (recD3isMsivectomatsX) %*% (MisXMX) -
          t(Gis2XdMs * recD2is) %*% (recD3isMsivectomatsXdMs) %*% (MisXMX) +
          t(Gis2 * XMis * recD2is) %*% (recD3isMs2 * ivectomatsXdMs) %*% (XMX)
        C2A44vecs[g] <- t(ds$X * XdMs) %*% (recD2sGs2Pgs) %*% (XMX) -
          t(ds$X * XMX) %*% (recD2sGs2Pgs) %*% (XMis)
      }
    }
    if (noisy) {
      cat(iteration, "of", max(df$groupW), "done. ")
      iteration <- iteration + 1
    }
  }

  C0A1 <- C0A11vecs - C0A12vecs - C0A13vecs + C0A14vecs + C0A15vecs
  C0A2 <- C0A21vecs - C0A22vecs - C0A23vecs + C0A24vecs + C0A25vecs
  C0A3 <- C0A31vecs - C0A32vecs - C0A33vecs + C0A34vecs + C0A35vecs
  C0A4 <- C0A41vecs + C0A42vecs + C0A43vecs + C0A44vecs
  C0A5 <- C0A51vecs + C0A52vecs + C0A53vecs + C0A54vecs


  C11A1 <- C11A11vecs - C11A12vecs - C11A13vecs + C11A14vecs + C11A15vecs
  C11A2 <- C11A21vecs - C11A22vecs - C11A23vecs + C11A24vecs + C11A25vecs
  C11A4 <- C11A41vecs + C11A42vecs + C11A43vecs + C11A44vecs
  C11A5 <- C11A51vecs + C11A52vecs + C11A53vecs + C11A54vecs

  C12A1 <- C12A11vecs - C12A12vecs - C12A13vecs + C12A14vecs + C12A15vecs
  C12A3 <- C12A31vecs - C12A32vecs - C12A33vecs + C12A34vecs + C12A35vecs
  C12A4 <- C11A5
  C12A5 <- C12A51vecs + C12A52vecs + C12A53vecs + C12A54vecs

  C2A1 <- C2A11vecs - C2A12vecs - C2A13vecs + C2A14vecs + C2A15vecs
  C2A4 <- C2A41vecs + C2A42vecs + C2A43vecs + C2A44vecs

  C0 <- sum(C0A1 + 2 * C0A2 + C0A3 - C0A4 - C0A5)
  C11 <- C11A1 + 3 * C11A2 - C11A4 - C11A5
  C12 <- 3 * C12A1 + C12A3 - C12A4 - C12A5
  C1 <- -sum(C11 + C12)
  C2 <- sum(4 * C2A1 - 2 * C2A4)

  PXY <- GetLM(df, X, Y, groupW, group)
  PXX <- GetLM(df, X, X, groupW, group)


  acon <- PXX^2 - q * C2
  bcon <- -2 * PXY * PXX - q * C1
  ccon <- PXY^2 - q * C0

  c(acon, bcon, ccon)
}

#' Compute UJIVE Signal Component (General Design)
#'
#' @description
#' Calculates the UJIVE "signal" or cross-product term \eqn{X' G Y} for a general design with
#' covariates. The weighting matrix \eqn{G} is defined as \eqn{U(P_Q) - U(P_W)}, where \eqn{U(P)}
#' represents the projection matrix with diagonal elements set to zero and rows rescaled by \eqn{1/(1-P_{ii})}.
#'
#' @param df Data frame. Contains the observable variables.
#' @param IdPW Numeric vector. The annihilator diagonals \eqn{1 - P_{W,ii}}.
#' @param IdPQ Numeric vector. The annihilator diagonals \eqn{1 - P_{Q,ii}}.
#' @param dPW Numeric vector. The projection diagonals \eqn{P_{W,ii}}.
#' @param dPQ Numeric vector. The projection diagonals \eqn{P_{Q,ii}}.
#' @param W Matrix. The covariate matrix.
#' @param Q Matrix. The combined instrument and covariate matrix \eqn{[Z, W]}.
#' @param X Column name (unquoted). The first variable (e.g., regressor).
#' @param Y Column name (unquoted). The second variable (e.g., outcome).
#'
#' @details
#' This function computes the scalar:
#' \deqn{S = X' [ (I-D_{P_Q})^{-1}(P_Q - D_{P_Q}) - (I-D_{P_W})^{-1}(P_W - D_{P_W}) ] Y}
#'
#'
#' It uses an efficient matrix algebra expansion that avoids constructing the \eqn{N \times N}
#' projection matrices explicitly. The computation complexity is linear in \eqn{N} (given pre-computed
#' basis matrices), making it suitable for large datasets where \eqn{G} is dense.
#'
#' The result corresponds to the cross-term \eqn{P_{XY}} (if \eqn{X \neq Y}) or the self-term \eqn{P_{XX}} (if \eqn{X = Y})
#' used in the confidence interval inequality.
#'
#' @return Numeric scalar. The computed cross-product.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @export
GetLM_WQ <- function(df, IdPW, IdPQ, dPW, dPQ, W, Q, X, Y) {
  df$Xpos <- eval(substitute(X), df)
  df$Ypos <- eval(substitute(Y), df)

  WW <- t(W) %*% W
  WWinv <- solve(WW)
  QQ <- t(Q) %*% Q
  QQinv <- solve(QQ)
  QX <- t(Q) %*% df$Xpos
  QY <- t(Q) %*% df$Ypos
  WX <- t(W) %*% df$Xpos
  WY <- t(W) %*% df$Ypos

  QXdQ <- t(Q) %*% (df$Xpos / IdPQ)
  QYdQ <- t(Q) %*% (df$Ypos / IdPQ)
  WXdW <- t(W) %*% (df$Xpos / IdPW)
  WYdW <- t(W) %*% (df$Ypos / IdPW)


  (t(QXdQ) %*% QQinv %*% QY - sum(df$Xpos * dPQ / IdPQ * df$Ypos)) -
    (t(WXdW) %*% WWinv %*% WY - sum(df$Xpos * dPW / IdPW * df$Ypos))
}

#' Compute UJIVE Signal Component (No Covariates)
#'
#' @description
#' Calculates the UJIVE "signal" or cross-product term \eqn{X' G e} for the grouped instrument design.
#' This function exploits the block-diagonal structure of the projection matrix implied by
#' mutually exclusive instrument groups to compute the quadratic form via efficient group-wise summation.
#'
#'
#' @param df Data frame. Contains the observable variables and grouping indicator.
#' @param X Column name (unquoted). The first variable (e.g., endogenous regressor).
#' @param e Column name (unquoted). The second variable (e.g., outcome or residual).
#' @param groupZ Column name (unquoted). The instrument grouping variable.
#' @param noisy Logical. If \code{TRUE}, prints progress of the group iteration.
#'   Defaults to \code{FALSE}.
#'
#' @details
#' This function computes the scalar:
#' \deqn{S = \sum_{g=1}^J e_g' [ (I-D_{P_g})^{-1}(P_g - D_{P_g}) ] X_g}
#' where \eqn{P_g} is the projection matrix onto the intercept for group \eqn{g} (i.e., the group mean).
#'
#' This calculation corresponds to the numerator components (\eqn{P_{XY}, P_{XX}}) of the
#' score statistic variance estimator in designs without covariates.
#'
#' @return Numeric scalar. The total sum of group-specific quadratic forms.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @export
GetLM_nocov <- function(df, X, e, groupZ, noisy = FALSE) {
  df$groupZ <- eval(substitute(groupZ), df)
  LMvecs <- rep(0, max(df$groupZ))
  df$Xpos <- eval(substitute(X), df)
  df$epos <- eval(substitute(e), df)

  iteration <- 1
  for (s in unique(df$groupZ)) {
    ds <- df[df$groupZ == s, ]
    ZQ <- matrix(0, nrow = length(ds$groupZ), ncol = length(unique(ds$groupZ)))
    ds$groupidx <- Getgroupindex(ds, groupZ)
    ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1

    PQ <- ZQ %*% solve(t(ZQ) %*% ZQ) %*% t(ZQ)

    # calculate values specific to this subset
    Gs <- diag(1 / (diag(diag(nrow(ds))) - pmin(diag(PQ), .99))) %*% (PQ - diag(diag(PQ)))
    LMvecs[s] <- t(ds$epos) %*% Gs %*% ds$Xpos

    if (noisy) {
      cat(iteration, "of", max(df$groupZ), "done. ")
      iteration <- iteration + 1
    }
  }
  sum(LMvecs)
}

#' Compute Covariance Matrix of UJIVE Quadratic Forms
#'
#' @description
#' Estimates the joint variance-covariance matrix of the three core quadratic forms used in
#' UJIVE/LIML estimation: \eqn{Y'GY}, \eqn{X'GY}, and \eqn{X'GX}. This function is highly
#' optimized for large-scale datasets by leveraging block-diagonal geometries and scalar algebra.
#'
#' @param df Data frame. Contains the variables used in estimation.
#' @param groupW Column name (unquoted). The covariate stratification variable.
#' @param group Column name (unquoted). The instrument grouping variable.
#' @param X Column name (unquoted). The endogenous regressor.
#' @param Y Column name (unquoted). The outcome variable.
#' @param MX Column name (unquoted). Leverage-adjusted regressor (\eqn{M X}).
#' @param MY Column name (unquoted). Leverage-adjusted outcome (\eqn{M Y}).
#' @param noisy Logical. If \code{TRUE}, prints progress during variance component calculation.
#'   Defaults to \code{FALSE}.
#'
#' @details
#' The function estimates the covariance components for the vector of quadratic forms:
#' \deqn{\Psi = [Y'GY, \quad X'GY, \quad X'GX]^T}
#'
#' \strong{Algorithmic Implementation:}
#' To achieve high performance, the function processes the data using a two-level nested loop
#' over covariate strata (\code{groupW}) and instrument groups (\code{group}).
#' The heavy \eqn{N \times N} geometry matrices (such as the projection matrix \eqn{P}, the
#' leverage-adjusted weight matrix \eqn{G}, and the residual variance matrix \eqn{D_2}) are
#' pre-computed exactly once per stratum.
#'
#' Within each instrument group, the function relies on localized extraction and exact scalar
#' algebra sub-routines (via \eqn{A_1} and \eqn{A_4} component helpers) to evaluate the
#' leave-three-out (L3O) variance interactions. This approach minimizes computational complexity
#' to \eqn{O(N)} at the group level and completely eliminates redundant matrix allocations and inversions.
#'
#' The returned vector contains the unique elements of the symmetric covariance matrix \eqn{\Sigma_\Psi}:
#' \itemize{
#'   \item \code{sig11}: \eqn{Var(Y'GY)}
#'   \item \code{sig22}: \eqn{Var(X'GY)}
#'   \item \code{sig33}: \eqn{Var(X'GX)}
#'   \item \code{sig12}: \eqn{Cov(Y'GY, X'GY)}
#'   \item \code{sig23}: \eqn{Cov(X'GY, X'GX)}
#'   \item \code{sig13}: \eqn{Cov(Y'GY, X'GX)}
#' }
#'
#' @return Numeric vector of length 6. Contains \code{c(sig11, sig22, sig33, sig12, sig23, sig13)}.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @export
GetSigMx <- function(df, groupW, group, X, Y, MX, MY, noisy = FALSE) {

  # 1. Evaluation
  df$X <- eval(substitute(X), df)
  df$Y <- eval(substitute(Y), df)
  df$MX <- eval(substitute(MX), df)
  df$MY <- eval(substitute(MY), df)
  df$groupW <- eval(substitute(groupW), df)
  df$group <- eval(substitute(group), df)

  s11_val <- s22_val <- s33_val <- 0
  s12_val <- s23_val <- s13_val <- 0
  iteration <- 1

  # 2. Outer Loop (Covariate Blocks)
  for (s in unique(df$groupW)) {
    ds <- df[df$groupW == s, ]
    ds$numingrp <- 0
    for (j in unique(ds$group)) ds$numingrp[ds$group == j] <- sum(ds$group == j)
    ds <- ds[ds$numingrp >= 3, ]

    if (nrow(ds) > 0) {
      N <- nrow(ds)

      # Block Pre-computation
      ZQ <- matrix(0, nrow = N, ncol = length(unique(ds$group)))
      ds$groupidx <- as.integer(factor(ds$group))
      ZQ[cbind(seq_len(N), ds$groupidx)] <- 1
      PQ <- ZQ %*% solve(t(ZQ) %*% ZQ) %*% t(ZQ)

      ZWmat <- matrix(1, nrow = length(ds$groupW), ncol = length(ds$groupW))
      PW <- ZWmat / (length(ds$groupW))

      Gs <- diag(1 / (diag(diag(N)) - pmin(diag(PQ), .999))) %*% (PQ - diag(diag(PQ))) -
        diag(1 / (diag(diag(N)) - pmin(diag(PW), .999))) %*% (PW - diag(diag(PW)))

      Ps <- PQ
      Ms <- diag(N) - Ps
      dMs <- matrix(diag(Ms), ncol = 1)
      D2s <- dMs %*% t(dMs) - Ms * Ms
      recD2s <- 1 / D2s
      diag(recD2s) <- 0

      recD2sGs2Pgs_base <- recD2s * Gs^2
      Gs2MsrecD2sPgs_base <- Gs^2 * Ms * recD2s

      vX <- ds$X; vY <- ds$Y; vMX <- ds$MX; vMY <- ds$MY

      # 3. Inner Loop (Instrument Groups)
      for (g in unique(ds$group)) {
        idx_g <- which(ds$group == g)

        if (length(idx_g) > 3) {
          repidx <- idx_g[1]

          # Group Vectors
          Pis <- matrix(Ps[, repidx], ncol = 1)
          Pgs <- ifelse(Pis == 0, 0, 1) %*% matrix(1, ncol = length(Pis), nrow = 1)

          # Gis[j] stands for G[i', j] for any i' in the instrument group g.
          # Claim 3 (Appendix H.1 rewrite): within a Q-group, G[i', j] is the
          # same constant g_g for every i' in g and any fixed j in g, j != i'.
          # G[repidx, repidx] = 0 by UJIVE construction, so we must overwrite
          # that entry with the in-group off-diagonal value g_g. Any other
          # element of the group works; idx_g[2] is guaranteed to be a valid
          # in-group neighbour regardless of how rows of `ds` are ordered.
          # (The original form `Gis[repidx + 1]` silently assumes that `ds`
          #  stores rows of the same instrument group contiguously.)
          Gis <- matrix(Gs[, repidx], ncol = 1)
          Gis[repidx, 1] <- Gis[idx_g[2], 1]

          # Mis is defined as -Ps[, repidx], i.e., the repidx column of -P.
          # Away from the repidx row, this coincides with M[, repidx]. At the
          # repidx row itself, Mis[repidx] = -P[repidx, repidx] differs from
          # the true M[repidx, repidx] = 1 - P[repidx, repidx] by exactly -1.
          # This offset is DELIBERATE: it makes the downstream D3is matrix
          # satisfy the identity
          #   D3is[repidx, k] = D_{i', repidx, k}  for every i' in idx_g \ {repidx, k}
          # (Lemma in Appendix H.4.1 of the rewrite), so a single D3is built
          # from repidx correctly encodes the per-i weight for every outer i
          # in the group, including at the repidx row/column.
          Mis <- matrix(-Ps[, repidx], ncol = 1)

          # Determinant Matrices
          D3is <- Ms[repidx, repidx] * D2s - (dMs %*% t(Mis)^2 + Mis^2 %*% t(dMs) - 2 * Ms * (Mis %*% t(Mis)))
          recD3is <- 1 / D3is
          diag(recD3is) <- 0
          D2D3is <- D2s / D3is
          diag(D2D3is) <- 0

          # Pre-multiply Ms into weights
          D2D3is_Ms <- D2D3is * Ms
          recD3is_Ms <- recD3is * Ms
          recD3is_Ms2 <- recD3is_Ms * Ms

          recD2is <- matrix(recD2s[, repidx], ncol = 1)
          # Same rationale as the Gis fix-up: use idx_g[2] rather than
          # repidx + 1 so that correctness does not depend on row ordering.
          recD2is[repidx, 1] <- recD2is[idx_g[2], 1]
          Gis2 <- Gis^2

          recD2sGs2Pgs <- recD2sGs2Pgs_base * Pgs
          Gs2MsrecD2sPgs <- Gs2MsrecD2sPgs_base * Pgs

          # --- Call External Helpers ---
          s11_val <- s11_val + (4 * compute_A1_scalar_ext(vY, vY, vY, vMY, idx_g, Gis, Mis, dMs, D2D3is, recD3is, recD3is_Ms, recD2sGs2Pgs, Gs2MsrecD2sPgs)
                                - 2 * compute_A4_scalar_ext(vY, vY, vY, vMY, idx_g, Gis2, Mis, dMs, recD2is, D2D3is, D2D3is_Ms, recD3is, recD3is_Ms, recD3is_Ms2, recD2sGs2Pgs))

          s33_val <- s33_val + (4 * compute_A1_scalar_ext(vX, vX, vX, vMX, idx_g, Gis, Mis, dMs, D2D3is, recD3is, recD3is_Ms, recD2sGs2Pgs, Gs2MsrecD2sPgs)
                                - 2 * compute_A4_scalar_ext(vX, vX, vX, vMX, idx_g, Gis2, Mis, dMs, recD2is, D2D3is, D2D3is_Ms, recD3is, recD3is_Ms, recD3is_Ms2, recD2sGs2Pgs))

          s22_val <- s22_val + (compute_A1_scalar_ext(vY, vX, vX, vMY, idx_g, Gis, Mis, dMs, D2D3is, recD3is, recD3is_Ms, recD2sGs2Pgs, Gs2MsrecD2sPgs)
                                + 2*compute_A1_scalar_ext(vY, vX, vY, vMX, idx_g, Gis, Mis, dMs, D2D3is, recD3is, recD3is_Ms, recD2sGs2Pgs, Gs2MsrecD2sPgs)
                                + compute_A1_scalar_ext(vX, vY, vY, vMX, idx_g, Gis, Mis, dMs, D2D3is, recD3is, recD3is_Ms, recD2sGs2Pgs, Gs2MsrecD2sPgs)
                                - compute_A4_scalar_ext(vY, vX, vY, vMX, idx_g, Gis2, Mis, dMs, recD2is, D2D3is, D2D3is_Ms, recD3is, recD3is_Ms, recD3is_Ms2, recD2sGs2Pgs)
                                - compute_A4_scalar_ext(vY, vY, vX, vMX, idx_g, Gis2, Mis, dMs, recD2is, D2D3is, D2D3is_Ms, recD3is, recD3is_Ms, recD3is_Ms2, recD2sGs2Pgs))

          # 12, 23, 13
          s12_val <- s12_val + (compute_A1_scalar_ext(vX, vY, vY, vMY, idx_g, Gis, Mis, dMs, D2D3is, recD3is, recD3is_Ms, recD2sGs2Pgs, Gs2MsrecD2sPgs)
                                + 3*compute_A1_scalar_ext(vY, vX, vY, vMY, idx_g, Gis, Mis, dMs, D2D3is, recD3is, recD3is_Ms, recD2sGs2Pgs, Gs2MsrecD2sPgs)
                                + compute_A4_scalar_ext(vX, vY, vY, vMY, idx_g, Gis2, Mis, dMs, recD2is, D2D3is, D2D3is_Ms, recD3is, recD3is_Ms, recD3is_Ms2, recD2sGs2Pgs)
                                - 3*compute_A4_scalar_ext(vY, vX, vY, vMY, idx_g, Gis2, Mis, dMs, recD2is, D2D3is, D2D3is_Ms, recD3is, recD3is_Ms, recD3is_Ms2, recD2sGs2Pgs))

          s23_val <- s23_val + (compute_A1_scalar_ext(vX, vX, vX, vMY, idx_g, Gis, Mis, dMs, D2D3is, recD3is, recD3is_Ms, recD2sGs2Pgs, Gs2MsrecD2sPgs)
                                + 3*compute_A1_scalar_ext(vX, vX, vY, vMX, idx_g, Gis, Mis, dMs, D2D3is, recD3is, recD3is_Ms, recD2sGs2Pgs, Gs2MsrecD2sPgs)
                                + compute_A4_scalar_ext(vX, vX, vX, vMY, idx_g, Gis2, Mis, dMs, recD2is, D2D3is, D2D3is_Ms, recD3is, recD3is_Ms, recD3is_Ms2, recD2sGs2Pgs)
                                - 3*compute_A4_scalar_ext(vX, vX, vY, vMX, idx_g, Gis2, Mis, dMs, recD2is, D2D3is, D2D3is_Ms, recD3is, recD3is_Ms, recD3is_Ms2, recD2sGs2Pgs))

          s13_val <- s13_val + (4*compute_A1_scalar_ext(vX, vX, vY, vMY, idx_g, Gis, Mis, dMs, D2D3is, recD3is, recD3is_Ms, recD2sGs2Pgs, Gs2MsrecD2sPgs)
                                - 2*compute_A4_scalar_ext(vX, vX, vY, vMY, idx_g, Gis2, Mis, dMs, recD2is, D2D3is, D2D3is_Ms, recD3is, recD3is_Ms, recD3is_Ms2, recD2sGs2Pgs))
        }
      }
    }
    if (noisy) { cat(iteration, "of", length(unique(df$groupW)), "done.\n"); iteration <- iteration + 1 }
  }
  c(s11_val, s22_val, s33_val, s12_val, s23_val, s13_val)
}

#' Compute Quadratic Coefficients for Confidence Sets (Grouped Data)
#'
#' @description
#' Estimates the coefficients \eqn{a}, \eqn{b}, and \eqn{c} for the quadratic inequality
#' \eqn{a\beta^2 + b\beta + c \le 0}, which defines the \eqn{1-\alpha} confidence set for the
#' structural parameter \eqn{\beta}. This function is highly optimized for large-scale datasets,
#' relying on block-diagonal geometries and scalar algebra.
#'
#' @param df Data frame. Contains the observable variables and their projections.
#' @param groupW Column name (unquoted). The covariate stratification variable.
#' @param group Column name (unquoted). The instrument grouping variable.
#' @param X Column name (unquoted). The endogenous regressor.
#' @param Y Column name (unquoted). The outcome variable.
#' @param MX Column name (unquoted). Leverage-adjusted regressor (\eqn{M X}).
#' @param MY Column name (unquoted). Leverage-adjusted outcome (\eqn{M Y}).
#' @param q Numeric scalar. Critical value for the test statistic inversion (e.g., \eqn{\chi^2_{1, 1-\alpha}}).
#'   Defaults to \code{qnorm(.975)^2} (approx. 3.84) for a 95 percent confidence interval.
#' @param noisy Logical. If \code{TRUE}, prints progress dots during calculation.
#'   Defaults to \code{FALSE}.
#'
#' @details
#' The confidence set is constructed by inverting a test statistic based on the quadratic form
#' \eqn{Q(\beta) = (\mathbf{Y} - \beta \mathbf{X})' G (\mathbf{Y} - \beta \mathbf{X})}.
#' The coefficients are derived from the variance estimator \eqn{\hat{V}(\beta)} of this quadratic form,
#' decomposed into interactions between the outcome and the regressor.
#'
#' \strong{Algorithmic Implementation:}
#' To achieve high performance, the function processes the data using a two-level nested loop
#' over covariate strata (\code{groupW}) and instrument groups (\code{group}).
#' The heavy \eqn{N \times N} geometry matrices (such as the projection matrices \eqn{P}, \eqn{M},
#' and the leverage-adjusted weight matrix \eqn{G}) are pre-computed exactly once per stratum.
#'
#' Furthermore, the target inner products (\eqn{P_{XY}} and \eqn{P_{XX}}) are aggregated inline
#' during the stratum loop to avoid redundant passes over the data. Within each instrument group,
#' unique Leave-Three-Out (L3O) variance interactions (\eqn{A_1} and \eqn{A_4} terms) are evaluated
#' using exact scalar algebra helpers. This approach minimizes computational complexity
#' to \eqn{O(N)} at the group level and completely eliminates redundant matrix allocations.
#'
#' The returned coefficients correspond to:
#' \deqn{a = P_{XX}^2 - q \cdot C_2}
#' \deqn{b = -2 P_{XY} P_{XX} - q \cdot C_1}
#' \deqn{c = P_{XY}^2 - q \cdot C_0}
#'
#' Where \eqn{P_{XY}} and \eqn{P_{XX}} are the UJIVE estimators for the cross-products, and
#' \eqn{C_0, C_1, C_2} are the variance components compiled via the L3O adjustment framework.
#'
#' @return Numeric vector of length 3: \code{c(a, b, c)}.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @export
GetCIcoef <- function(df, groupW, group, X, Y, MX, MY,
                           q = qnorm(.975)^2, noisy = FALSE) {

  # 1. Evaluation  --------------------------------------------------------
  df$X      <- eval(substitute(X),      df)
  df$Y      <- eval(substitute(Y),      df)
  df$MX     <- eval(substitute(MX),     df)
  df$MY     <- eval(substitute(MY),     df)
  df$groupW <- eval(substitute(groupW), df)
  df$group  <- eval(substitute(group),  df)

  # Accumulators for the ten distinct A1/A4 scalars that appear across
  # C0, C1, C2. Naming: a1_ijkl / a4_ijkl with i,j,k in {X,Y} and l in
  # {MX,MY}. We only instantiate the ones actually needed.
  a1_YXX_MY <- 0; a1_YXY_MX <- 0; a1_XYY_MX <- 0     # C0 A1 terms
  a4_XYX_MY <- 0; a4_YYX_MX <- 0                     # C0 A4 terms

  a1_XXX_MY <- 0; a1_XXY_MX <- 0; a1_YXX_MX <- 0; a1_XYX_MX <- 0  # C1 A1 terms
  a4_XXX_MY <- 0; a4_XYX_MX <- 0; a4_YXX_MX <- 0                  # C1 A4 terms

  a1_XXX_MX <- 0; a4_XXX_MX <- 0                                  # C2 terms

  # GetLM accumulators: PXY = sum_s y' Gs x,   PXX = sum_s x' Gs x
  PXY <- 0
  PXX <- 0

  iteration <- 1

  # 2. Outer Loop (Covariate Blocks)  -------------------------------------
  for (s in unique(df$groupW)) {
    ds <- df[df$groupW == s, ]
    ds$numingrp <- 0
    for (j in unique(ds$group)) ds$numingrp[ds$group == j] <- sum(ds$group == j)
    ds <- ds[ds$numingrp >= 4, ]

    if (nrow(ds) > 0) {
      N <- nrow(ds)

      # Block Pre-computation  (shared with GetSigMx)
      ZQ <- matrix(0, nrow = N, ncol = length(unique(ds$group)))
      ds$groupidx <- as.integer(factor(ds$group))
      ZQ[cbind(seq_len(N), ds$groupidx)] <- 1
      PQ <- ZQ %*% solve(t(ZQ) %*% ZQ) %*% t(ZQ)

      ZWmat <- matrix(1, nrow = length(ds$groupW), ncol = length(ds$groupW))
      PW <- ZWmat / (length(ds$groupW))

      Gs <- diag(1 / (diag(diag(N)) - pmin(diag(PQ), .99))) %*% (PQ - diag(diag(PQ))) -
        diag(1 / (diag(diag(N)) - pmin(diag(PW), .99))) %*% (PW - diag(diag(PW)))

      Ps  <- PQ
      Ms  <- diag(N) - Ps
      dMs <- matrix(diag(Ms), ncol = 1)
      D2s <- dMs %*% t(dMs) - Ms * Ms
      recD2s <- 1 / D2s
      diag(recD2s) <- 0

      recD2sGs2Pgs_base  <- recD2s * Gs^2
      Gs2MsrecD2sPgs_base <- Gs^2 * Ms * recD2s

      vX <- ds$X; vY <- ds$Y; vMX <- ds$MX; vMY <- ds$MY

      # GetLM contribution from this block: y' Gs x  and  x' Gs x
      # (drop(...) guarantees a scalar even though R would otherwise
      # return a 1x1 matrix.)
      Gs_vX <- Gs %*% vX
      PXY <- PXY + drop(crossprod(vY, Gs_vX))
      PXX <- PXX + drop(crossprod(vX, Gs_vX))

      # 3. Inner Loop (Instrument Groups)  --------------------------------
      for (g in unique(ds$group)) {
        idx_g <- which(ds$group == g)

        if (length(idx_g) > 3) {
          repidx <- idx_g[1]

          # Group Vectors
          Pis <- matrix(Ps[, repidx], ncol = 1)
          Pgs <- ifelse(Pis == 0, 0, 1) %*% matrix(1, ncol = length(Pis), nrow = 1)

          # Gis fix-up: overwrite the structural zero at the repidx row
          # with the in-group off-diagonal constant. Using idx_g[2]
          # (rather than repidx+1) makes this independent of row order
          # within the dataframe.
          Gis <- matrix(Gs[, repidx], ncol = 1)
          Gis[repidx, 1] <- Gis[idx_g[2], 1]

          # Mis = -Ps[,repidx]; the -1 offset at the repidx row is
          # deliberate so that D3is encodes per-i weights for every
          # outer i in the group (see Appendix H.4.1 of the rewrite).
          Mis <- matrix(-Ps[, repidx], ncol = 1)

          # Determinant Matrices
          D3is <- Ms[repidx, repidx] * D2s -
            (dMs %*% t(Mis)^2 + Mis^2 %*% t(dMs) - 2 * Ms * (Mis %*% t(Mis)))
          recD3is <- 1 / D3is
          diag(recD3is) <- 0
          D2D3is <- D2s / D3is
          diag(D2D3is) <- 0

          # Pre-multiply Ms into weights
          D2D3is_Ms    <- D2D3is * Ms
          recD3is_Ms   <- recD3is * Ms
          recD3is_Ms2  <- recD3is_Ms * Ms

          recD2is <- matrix(recD2s[, repidx], ncol = 1)
          recD2is[repidx, 1] <- recD2is[idx_g[2], 1]
          Gis2 <- Gis^2

          recD2sGs2Pgs    <- recD2sGs2Pgs_base    * Pgs
          Gs2MsrecD2sPgs  <- Gs2MsrecD2sPgs_base  * Pgs

          # Bundle args so each A1 / A4 call stays readable
          a1 <- function(i, j, k, l) {
            compute_A1_scalar_ext(i, j, k, l, idx_g,
                                  Gis, Mis, dMs,
                                  D2D3is, recD3is, recD3is_Ms,
                                  recD2sGs2Pgs, Gs2MsrecD2sPgs)
          }
          a4 <- function(i, j, k, l) {
            compute_A4_scalar_ext(i, j, k, l, idx_g,
                                  Gis2, Mis, dMs, recD2is,
                                  D2D3is, D2D3is_Ms,
                                  recD3is, recD3is_Ms, recD3is_Ms2,
                                  recD2sGs2Pgs)
          }

          # --- C0 contributions -----------------------------------------
          a1_YXX_MY <- a1_YXX_MY + a1(vY, vX, vX, vMY)
          a1_YXY_MX <- a1_YXY_MX + a1(vY, vX, vY, vMX)
          a1_XYY_MX <- a1_XYY_MX + a1(vX, vY, vY, vMX)
          a4_XYX_MY <- a4_XYX_MY + a4(vX, vY, vX, vMY)
          a4_YYX_MX <- a4_YYX_MX + a4(vY, vY, vX, vMX)

          # --- C1 contributions -----------------------------------------
          # Note the reference computes two of these terms twice (with
          # coefficients 2 and 1); we collect each term once and apply
          # the combined coefficient (3) at the end.
          a1_XXX_MY <- a1_XXX_MY + a1(vX, vX, vX, vMY)
          a1_XXY_MX <- a1_XXY_MX + a1(vX, vX, vY, vMX)
          a1_YXX_MX <- a1_YXX_MX + a1(vY, vX, vX, vMX)
          a1_XYX_MX <- a1_XYX_MX + a1(vX, vY, vX, vMX)
          a4_XXX_MY <- a4_XXX_MY + a4(vX, vX, vX, vMY)
          a4_XYX_MX <- a4_XYX_MX + a4(vX, vY, vX, vMX)
          a4_YXX_MX <- a4_YXX_MX + a4(vY, vX, vX, vMX)

          # --- C2 contributions -----------------------------------------
          a1_XXX_MX <- a1_XXX_MX + a1(vX, vX, vX, vMX)
          a4_XXX_MX <- a4_XXX_MX + a4(vX, vX, vX, vMX)
        }
      }
    }
    if (noisy) {
      cat(iteration, "of", length(unique(df$groupW)), "done.\n")
      iteration <- iteration + 1
    }
  }

  # 4. Assemble C0, C1, C2  ------------------------------------------------
  C0 <-  a1_YXX_MY + 2 * a1_YXY_MX + a1_XYY_MX -
    a4_XYX_MY - a4_YYX_MX

  C1 <- -( a1_XXX_MY +
             3 * a1_XXY_MX -                 # coeffs 2 + 1 in the reference
             a4_XXX_MY - a4_XYX_MX +
             3 * a1_YXX_MX +                 # coeffs 1 + 2 in the reference
             a1_XYX_MX -
             a4_XYX_MX - a4_YXX_MX )

  C2 <- 4 * a1_XXX_MX - 2 * a4_XXX_MX

  # 5. Final coefficients  -------------------------------------------------
  acon <-  PXX^2        - q * C2
  bcon <- -2 * PXY * PXX - q * C1
  ccon <-  PXY^2        - q * C0

  c(acon, bcon, ccon)
}



#' Compute UJIVE Signal Component (Stratified Design)
#'
#' @description
#' Calculates the UJIVE "signal" or cross-product term \eqn{X' G e} for a design where
#' instruments are nested within discrete covariate strata (e.g., Judges within Years).
#' This function iterates through covariate groups to compute the quadratic form locally,
#' handling the centering of instruments within each block.
#'
#' @param df Data frame. Contains the observable variables and grouping indicators.
#' @param X Column name (unquoted). The first variable (e.g., endogenous regressor).
#' @param e Column name (unquoted). The second variable (e.g., outcome or residual).
#' @param groupW Column name (unquoted). The covariate stratification variable (defines blocks).
#' @param group Column name (unquoted). The instrument grouping variable (defines treatments).
#' @param noisy Logical. If \code{TRUE}, prints progress of the stratum iteration.
#'   Defaults to \code{FALSE}.
#'
#' @details
#' This function implements the estimator for the signal component \eqn{S} in a stratified design:
#' \deqn{S = \sum_{s} e_s' [ U(P_{Q,s}) - U(P_{W,s}) ] X_s}
#'
#' Within each stratum \eqn{s}:
#' \itemize{
#'   \item \eqn{P_{Q,s}} is the projection onto the instrument groups.
#'   \item \eqn{P_{W,s}} is the projection onto the stratum intercept (local mean).
#'   \item \eqn{U(P)} denotes the projection matrix with its diagonal elements removed,
#'   rescaled by the inverse of the annihilator diagonal \eqn{(1 - P_{ii})^{-1}}.
#' }
#'
#' This corresponds to the numerator terms (\eqn{P_{XY}, P_{XX}}) for test inversion in
#' designs with discrete controls.
#'
#' \strong{Algorithmic Implementation:}
#' To ensure high performance and low memory overhead on large datasets, this function computes
#' the projection transformation \eqn{G_s X_s} using strictly \eqn{O(N)} vector operations.
#' It computes group means via optimized aggregation, derives the diagonal leverage adjustments
#' inline, and directly applies the leave-one-out transformation without ever constructing
#' the dense \eqn{N \times N} projection matrices.
#'
#' @return Numeric scalar. The sum of stratum-specific quadratic forms.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @export
GetLM <- function(df, X, e, groupW, group, noisy = FALSE) {

  df$Xpos   <- eval(substitute(X),      df)
  df$epos   <- eval(substitute(e),      df)
  df$groupW <- eval(substitute(groupW), df)
  df$group  <- eval(substitute(group),  df)

  total     <- 0
  uniqW     <- unique(df$groupW)
  iteration <- 1L

  for (s in uniqW) {
    idx_s <- which(df$groupW == s)
    N     <- length(idx_s)
    if (N == 0L) next

    x  <- df$Xpos[idx_s]
    ev <- df$epos[idx_s]

    # --- Instrument-group structure within this W block --------------------
    gvec   <- df$group[idx_s]
    gfac   <- factor(gvec)                 # stable mapping to 1..J
    gidx   <- as.integer(gfac)             # length N, values in 1..J
    nj     <- tabulate(gidx)               # group sizes, length J
    n_i    <- nj[gidx]                     # group size for each row, length N

    # --- Group means of x -------------------------------------------------
    # sum_j = sum of x over rows with group index j; length J
    sum_j  <- as.numeric(tapply(x, gidx, sum))
    # If a group is empty in this block it won't appear in gidx; tapply
    # returns length = length(unique(gidx)) in the order of sort(unique(gidx))
    # which equals 1..J here because gidx was built from factor(). So
    # sum_j is aligned with 1..J by construction.
    PQx_i  <- (sum_j / nj)[gidx]           # (PQ x)_i,  length N

    # --- diag(PQ) and diag(PW) with the .99 clamp -------------------------
    dPQ_i  <- 1 / n_i                      # diag(PQ) for each row
    dPW    <- 1 / N                        # scalar, same for every row
    dPQ_c  <- pmin(dPQ_i, .99)
    dPW_c  <- min(dPW,   .99)

    DQ_i   <- 1 / (1 - dPQ_c)              # length N
    DW     <- 1 / (1 - dPW_c)              # scalar

    # --- Grand mean of x (= PW x) -----------------------------------------
    mean_x <- sum(x) / N                   # scalar

    # --- (Gs x) in O(N), then contract with e -----------------------------
    # Gs x = DQ * (PQ x - diag(PQ) * x) - DW * (PW x - diag(PW) * x)
    Gs_x   <- DQ_i * (PQx_i - dPQ_i * x) - DW * (mean_x - dPW * x)

    total  <- total + sum(ev * Gs_x)

    if (noisy) {
      cat(iteration, "of", length(uniqW), "done. ")
      iteration <- iteration + 1L
    }
  }

  total
}
