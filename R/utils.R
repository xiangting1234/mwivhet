#' Create a Matrix by Repeating a Vector Row-wise
#'
#' @description
#' A helper function that replicates a vector \code{X} into a matrix with \code{n} identical rows.
#' In linear algebra terms, this creates the rank-1 matrix \eqn{\mathbf{1}_n \otimes X'} (Kronecker product).
#'
#' @details
#' This is primarily used for "broadcasting" a vector across a matrix operation without using
#' slower `apply` loops. For example, subtracting a mean vector from every row of a data matrix.
#'
#' @param X Numeric vector. The vector to be repeated.
#' @param n Integer. The number of rows in the resulting matrix.
#'
#' @returns A numeric matrix of dimension \code{c(n, length(X))}, where every row is identical to \code{X}.
#' @noRd
bindrowvecs <- function(X, n) {
  matrix(rep(X, n), nrow = n, byrow = TRUE)
}

#' Create a Matrix by Repeating a Vector Column-wise
#'
#' @description
#' A helper function that replicates a vector \code{X} into a matrix of \code{n} rows.
#' Unlike \code{bindrowvecs}, this fills the matrix column-wise, ensuring that every column
#' is a recycled copy of \code{X}.
#'
#' @details
#' In the specific context of the MWIV package (where \code{length(X) == n}), this creates
#' a square \eqn{n \times n} matrix where every column is identical to the vector \code{X}.
#' Algebraically, this corresponds to the outer product \eqn{X \mathbf{1}_n'}.
#'
#' If \code{length(X) != n}, the function creates a matrix of dimension \eqn{n \times \text{length}(X)},
#' where each column consists of \code{X} repeated/recycled to fill \code{n} rows.
#'
#' @param X Numeric vector. The vector to be repeated (typically of length \eqn{n}).
#' @param n Integer. The number of rows in the resulting matrix.
#'
#' @returns A numeric matrix where every column contains the elements of \code{X}.
#' @noRd
bindcolvecs <- function(X, n) {
  matrix(rep(X, n), nrow = n, byrow = FALSE)
}

#' Compute Grouped Leave-Two-Out Sums
#'
#' @description
#' Calculates a matrix of leave-two-out sums for a specific group. This is a helper
#' function used to construct variance estimators that require removing the contributions
#' of observation pairs \eqn{(i,j)} from group-level totals.
#'
#' @param ds Data frame containing the grouping variable. The function assumes
#'   non-standard evaluation of a column named `group`.
#' @param X Numeric vector of length \eqn{n}. The values to be summed (e.g., instrument values).
#' @param g Scalar. The specific group identifier (level) to compute the sums for.
#'
#' @details
#' For a given group \eqn{g}, this function computes an \eqn{n \times n} matrix where the
#' \eqn{(i,j)}-th element is:
#' \deqn{S_{ij} = \sum_{k \in g} X_k - X_i \mathbb{I}(i \in g) - X_j \mathbb{I}(j \in g)}
#' This efficiently calculates the sum of \code{X} for group \code{g} while excluding
#' observations \eqn{i} and \eqn{j}.
#'
#' @return A matrix of dimension \eqn{n \times n}.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity".
#' Working Paper.
#'
#' @noRd
ivectomats <- function(ds, X, g) {
  ds$group <- eval(substitute(group), ds)
  gidx <- which(ds$group == g)
  sumq <- sum(X[gidx])
  subtractq <- matrix(rep(0, nrow(ds)^2), nrow = nrow(ds))
  subtractq[gidx, ] <- X[gidx]
  subtractq[, gidx] <- subtractq[, gidx] + t(subtractq[gidx, ])
  out <- (matrix(rep(1, nrow(ds)^2), nrow = nrow(ds)) - diag(nrow(ds))) * sumq - subtractq
  out
}

#' Generate Integer Group Indices
#'
#' @description
#' Converts a grouping variable into a vector of contiguous integer indices ranging
#' from 1 to \eqn{G}, where \eqn{G} is the number of unique groups. This is a helper
#' function used to facilitate matrix indexing for group-level variance calculations.
#'
#' @param ds Data frame. The dataset containing the grouping variable.
#' @param group Name of the grouping variable (passed as an unquoted symbol) within \code{ds}.
#'   This variable defines the clusters (e.g., judges, examiners, or time periods).
#'
#' @details
#' The function extracts the specified \code{group} column from the data frame using
#' non-standard evaluation. It then maps each unique value of the group to an integer
#' index corresponding to its order of appearance in \code{unique(group)}.
#'
#' This transformation is functionally equivalent to \code{as.numeric(as.factor(group))},
#' but ensures explicit handling within the provided data frame context.
#'
#' @return An integer vector of length \eqn{n}, containing values in \eqn{\{1, \dots, G\}}.
#'
#' @export
#'
#' @examples
#' data <- data.frame(judge_id = c("JudgeA", "JudgeB", "JudgeA", "JudgeC"))
#' Getgroupindex(data, judge_id)
Getgroupindex <- function(ds, group) {
  ds$group <- eval(substitute(group), ds)
  ds$groupidx <- 0
  gvals <- unique(ds$group)
  for (g in unique(ds$group)) {
    ds$groupidx[ds$group == g] <- which(unique(ds$group) == g)
  }
  ds$groupidx
}

#' Construct UJIVE Weighting and Projection Matrices
#'
#' @description
#' Generates the weighting matrix \eqn{G} and projection matrix \eqn{P} required for the
#' Unbiased Jackknife Instrumental Variables Estimator (UJIVE). It constructs these
#' matrices based on group/cluster indicators for instruments (\eqn{Z}) and
#' stratification covariates (\eqn{W}).
#'
#' @param group Integer vector. Primary grouping variable used to define the instrument structure (e.g., judge or examiner IDs).
#' @param groupW Integer vector. Grouping variable for stratification or covariates (e.g., time periods or court locations).
#' @param n Integer. The total sample size (number of observations).
#'
#' @details
#' The function constructs the design matrices for instruments (\eqn{Z}) and covariates (\eqn{W})
#' based on the provided grouping vectors. It creates dummy variable matrices where
#' \eqn{Z_{ij} = 1} if observation \eqn{i} belongs to group \eqn{j}.
#'
#' It calculates the standard projection matrices:
#' \deqn{P_Z = Z(Z'Z)^{-1}Z'}
#' \deqn{P_W = W(W'W)^{-1}W'}
#' \deqn{P = P_{[Z,W]} \quad (\text{Projection onto both } Z \text{ and } W)}
#'
#' The UJIVE weighting matrix \eqn{G} is then computed as the difference between the
#' "leave-one-out" adjusted projection of the full set and the covariates:
#' \deqn{G = D_P^{-1}(P - \text{diag}(P)) - D_W^{-1}(P_W - \text{diag}(P_W))}
#' where \eqn{D_P} and \eqn{D_W} are diagonal matrices containing the annihilator diagonals
#' (\eqn{1 - P_{ii}}) for the respective projections.
#'
#' Note: This function assumes that the resulting design matrices are full rank.
#'
#' @return A list containing four matrices:
#' \itemize{
#'   \item \code{G}: The \eqn{N \times N} UJIVE weighting matrix.
#'   \item \code{P}: The \eqn{N \times N} full projection matrix on \eqn{[Z, W]}.
#'   \item \code{Z}: The \eqn{N \times K} matrix of instrument indicators.
#'   \item \code{W}: The \eqn{N \times L} matrix of covariate indicators.
#' }
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @export
GetGP <- function(group, groupW, n) {
  groupZ <- groupW * (group %% 2)
  # Create Z matrix of indicators
  Z <- matrix(0, nrow = length(groupZ), ncol = length(unique(groupZ)) - 1)
  Z[cbind(seq_along(groupZ), groupZ)] <- 1

  ## Calculations related to Z
  ZZ <- t(Z) %*% Z
  ZZ_inv <- solve(ZZ)
  ZZ_inv2 <- chol(ZZ_inv)
  ZZZ_inv <- Z %*% ZZ_inv
  ZZZ_inv2 <- Z %*% ZZ_inv2 # n x k mx
  HZ <- Z %*% ZZ_inv %*% t(Z)

  ## Calculations related to W
  Wd <- matrix(0, nrow = length(groupW), ncol = length(unique(groupW)))
  Wd[cbind(seq_along(groupW), groupW)] <- 1
  # W <- cbind(Wd,Wcfix)
  W <- Wd
  WW <- t(W) %*% W
  WW_inv <- solve(WW)
  WW_inv2 <- chol(WW_inv)
  WWW_inv <- W %*% WW_inv
  WWW_inv2 <- W %*% WW_inv2
  HW <- W %*% WW_inv %*% t(W)
  MW <- diag(n) - HW

  ## Combine Z and W
  Q <- cbind(Z, W)
  QQ <- t(Q) %*% Q
  QQ_inv <- solve(QQ)
  QQ_inv2 <- chol(QQ_inv)
  QQQ_inv <- Q %*% QQ_inv
  QQQ_inv2 <- Q %*% QQ_inv2

  # P for full projection on Z and W
  P <- Q %*% QQ_inv %*% t(Q)
  M <- diag(n) - P

  ## G for UJIVE
  G <- solve(diag(n) - diag(diag(P))) %*% (P - diag(diag(P))) -
    solve(diag(n) - diag(diag(HW))) %*% (HW - diag(diag(HW)))

  list(G = G, P = P, Z = Z, W = W)
}



#' Fast Algebraic Computation of Group-Level Quadratic Forms
#'
#' @description
#' A high-performance helper function that analytically computes the quadratic form
#' \code{t(u) \%*\% (W_mat * ivectomats(ds, x_vec, g)) \%*\% v}. It circumvents the need
#' to construct the memory-intensive \eqn{n \times n} matrix produced by \code{ivectomats}.
#'
#' @details
#' In the estimation of variance components (e.g., $A_1$ through $A_5$), the original code
#' constructs a dense \eqn{n \times n} matrix via \code{ivectomats} for every group, leading
#' to severe memory bottlenecks and slow matrix multiplications.
#'
#' Because the \code{ivectomats} matrix has a highly specific rank-2 structure defined by
#' group sums and recycled vectors, the full quadratic form can be algebraically simplified
#' into scalar and vector dot-products.
#'
#' Assuming the diagonal of \code{W_mat} is zero (which is strictly enforced in the parent
#' variance functions for \code{D2D3is} and \code{recD3is}), the trace computation simplifies to:
#' \deqn{ \left( \sum_{k \in g} X_k \right) (u' W v) - \sum_{k \in g} X_k \left[ u_k (W v)_k + v_k (W u)_k \right] }
#'
#' This reduces the spatial complexity to \eqn{O(1)} (no matrix allocation) and the computational
#' complexity to the cost of the matrix-vector products \eqn{W u} and \eqn{W v}.
#'
#' @param u Numeric vector of length \eqn{n}. The left-hand vector in the quadratic form.
#' @param v Numeric vector of length \eqn{n}. The right-hand vector in the quadratic form.
#' @param x_vec Numeric vector of length \eqn{n}. The variable vector used to form the
#'   leave-two-out group sums (the \code{X} argument that would normally be passed to \code{ivectomats}).
#' @param W_mat Numeric matrix of dimension \eqn{n \times n}. The weighting matrix
#'   (typically \code{D2D3is} or \code{recD3is} derived from the stratum geometry). It is
#'   assumed to have zeros on its diagonal.
#' @param idx_g Integer vector. The row indices of the observations belonging to the current group \eqn{g}.
#'
#' @returns A numeric scalar representing the exact result of the matrix multiplication.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @noRd
calc_quad_fast <- function(u, v, x_vec, W_mat, idx_g) {
  # 1. Pre-calculate weighted vectors
  Wu <- c(W_mat %*% u)
  Wv <- c(W_mat %*% v)

  # 2. Term 1: The Summation Part
  sum_x <- sum(x_vec[idx_g])
  term1 <- sum_x * sum(u * Wv)

  # 3. Term 2: The Cross Part
  # Sum over k in group: x_k * [ u_k * (W v)_k + v_k * (W u)_k ]
  term2 <- sum(x_vec[idx_g] * (u[idx_g] * Wv[idx_g] + v[idx_g] * Wu[idx_g]))

  return(term1 - term2)
}


#' Compute the A1 Variance Component Scalar
#'
#' @description
#' A highly optimized helper function that calculates the \eqn{A_1} variance component
#' for a single group. It utilizes the fast algebraic quadratic form evaluator
#' \code{calc_quad_fast} to bypass the construction of large \eqn{n \times n} dense matrices.
#'
#' @details
#' In the estimation of the asymptotic variance matrix (the "meat" of the sandwich estimator,
#' \eqn{\Sigma}) under heterogeneous treatment effects, the term \eqn{A_1} captures the third-moment
#' or linear-interaction properties of the error terms.
#'
#' The original implementation calculated five complex quadratic forms (\code{v1} through \code{v5})
#' using \code{ivectomats}, which resulted in extreme memory overhead and \eqn{O(N^3)} time complexity
#' due to repeated dense matrix multiplications.
#'
#' This extended helper re-expresses these five terms using scalar math and pre-computed,
#' stratum-level geometric matrices. Specifically, it maps permutations of the data vectors
#' (\code{i, j, k, l})—which typically represent \eqn{X, Y, M_X, M_Y}—to the corresponding
#' weighting matrices (\code{D2D3is}, \code{recD3is}) generated by the leave-three-out (L3O)
#' geometry.
#'
#' Note on \code{v3}: The original code utilized an outer product \code{ones \%*\% t(l)} inside the trace.
#' By properties of the trace operator, this is algebraically equivalent to absorbing the vector
#' \code{l} into the right-hand vector \code{vR} (i.e., multiplying the right vector by \code{l}
#' element-wise), which avoids creating the \eqn{n \times n} outer product matrix.
#'
#' @param i Numeric vector. Data vector corresponding to \code{ipos} (e.g., \eqn{Y} or \eqn{X}).
#' @param j Numeric vector. Data vector corresponding to \code{jpos} (e.g., \eqn{X} or \eqn{Y}).
#' @param k Numeric vector. Data vector corresponding to \code{kpos}.
#' @param l Numeric vector. Data vector corresponding to \code{lpos} (typically a residual vector like \eqn{M_Y}).
#' @param idx_g Integer vector. The row indices belonging to the current group.
#' @param Gis Numeric matrix (single column). The group-specific leverage-adjusted weight vector.
#' @param Mis Numeric matrix (single column). The group-specific residual maker vector (\eqn{-P_s} for the group representative).
#' @param dMs Numeric matrix (single column). The diagonal elements of the residual maker matrix \eqn{M_s}.
#' @param D2D3is Numeric matrix. The group-specific inverse variance weight matrix for the L3O estimator.
#' @param recD3is Numeric matrix. The group-specific inverse weight matrix derived from \eqn{D_3}.
#' @param recD3is_Ms Numeric matrix. The element-wise product of \code{recD3is} and \eqn{M_s}.
#' @param recD2sGs2Pgs Numeric matrix. Pre-computed matrix: \code{recD2s * Gs^2 * Pgs}.
#' @param Gs2MsrecD2sPgs Numeric matrix. Pre-computed matrix: \code{Gs^2 * Ms * recD2s * Pgs}.
#'
#' @returns A numeric scalar representing the sum of the five sub-components
#'   (\eqn{v_1 - v_2 - v_3 + v_4 + v_5}) for the specified group.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @seealso \code{\link{calc_quad_fast}}
#'
#' @noRd
compute_A1_scalar_ext <- function(i, j, k, l, idx_g,
                                  Gis, Mis, dMs,
                                  D2D3is, recD3is, recD3is_Ms,
                                  recD2sGs2Pgs, Gs2MsrecD2sPgs) {

  # v1: Standard D2D3
  v1 <- calc_quad_fast(j * Gis, k * Gis, i * l, D2D3is, idx_g)

  # v2: Using recD3is
  v2 <- calc_quad_fast(j * l * Gis * Mis, Gis * k * dMs, i, recD3is, idx_g) -
    calc_quad_fast(j * Gis * Mis, Gis * k * l, i, recD3is_Ms, idx_g)

  # v3: "ones %*% t(l)" logic vectorized
  v3 <- calc_quad_fast(j * dMs * Gis, k * Gis * Mis * l, i, recD3is, idx_g)

  # v4: recD3is_Ms term
  v4 <- calc_quad_fast(j * l * Gis, k * Gis * Mis, i, recD3is_Ms, idx_g)

  # v5: Vector math
  v5 <- t(l * i) %*% recD2sGs2Pgs %*% (j * k * dMs) -
    t(i) %*% Gs2MsrecD2sPgs %*% (l * j * k)

  return(v1 - v2 - v3 + v4 + v5)
}


#' Compute the A4 Variance Component Scalar
#'
#' @description
#' A highly optimized helper function that calculates the \eqn{A_4} variance component
#' for a single group. It utilizes the fast algebraic quadratic form evaluator
#' \code{calc_quad_fast} to bypass the memory-intensive \eqn{n \times n} dense matrices
#' typically generated by \code{ivectomats}.
#'
#' @details
#' In the theoretical framework of Inference with Many Weak Instruments and Heterogeneity
#' (Yap, 2025), the asymptotic variance of the estimator requires corrections for non-constant
#' error variances (kurtosis). While the \eqn{A_1} term handles linear interaction properties,
#' the \eqn{A_4} (and \eqn{A_5}) terms represent the higher-order heterogeneity corrections.
#'
#' Mathematically, the \eqn{A_4} term is distinguished by its use of the squared leverage weights
#' (\eqn{G_{ij}^2}), passed into this function as \code{Gis2}.
#'
#' This extended helper evaluates the four sub-components (\code{v1} through \code{v4}) of \eqn{A_4}.
#' By substituting matrix-based \code{ivectomats} calls with \code{calc_quad_fast} and utilizing
#' pre-computed stratum-level matrices (\code{D2D3is_Ms}, \code{recD3is_Ms2}, etc.), it reduces
#' time complexity from \eqn{O(N^3)} to \eqn{O(N)} and spatial complexity to \eqn{O(1)} for each group loop.
#'
#' @param i Numeric vector. Data vector corresponding to \code{ipos} (e.g., \eqn{Y} or \eqn{X}).
#' @param j Numeric vector. Data vector corresponding to \code{jpos} (e.g., \eqn{X} or \eqn{Y}).
#' @param k Numeric vector. Data vector corresponding to \code{kpos}.
#' @param l Numeric vector. Data vector corresponding to \code{lpos} (typically a residual vector).
#' @param idx_g Integer vector. The row indices belonging to the current group.
#' @param Gis2 Numeric matrix (single column). The squared group-specific leverage-adjusted weight vector (\eqn{G_{is}^2}).
#' @param Mis Numeric matrix (single column). The group-specific residual maker vector.
#' @param dMs Numeric matrix (single column). The diagonal elements of the residual maker matrix \eqn{M_s}.
#' @param recD2is Numeric matrix (single column). Group-specific inverse variance weight vector from \eqn{D_2}.
#' @param D2D3is Numeric matrix. The group-specific inverse variance weight matrix for the L3O estimator.
#' @param D2D3is_Ms Numeric matrix. The element-wise product of \code{D2D3is} and \eqn{M_s}.
#' @param recD3is Numeric matrix. The group-specific inverse weight matrix derived from \eqn{D_3}.
#' @param recD3is_Ms Numeric matrix. The element-wise product of \code{recD3is} and \eqn{M_s}.
#' @param recD3is_Ms2 Numeric matrix. The element-wise product of \code{recD3is} and \eqn{M_s^2}.
#' @param recD2sGs2Pgs Numeric matrix. Pre-computed matrix: \code{recD2s * Gs^2 * Pgs}.
#'
#' @returns A numeric scalar representing the sum of the four sub-components
#'   (\eqn{v_1 + v_2 + v_3 + v_4}) for the specified group.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @seealso \code{\link{calc_quad_fast}}, \code{\link{compute_A1_scalar_ext}}
#'
#' @noRd
compute_A4_scalar_ext <- function(i, j, k, l, idx_g,
                                  Gis2, Mis, dMs, recD2is,
                                  D2D3is, D2D3is_Ms,
                                  recD3is, recD3is_Ms, recD3is_Ms2,
                                  recD2sGs2Pgs) {

  # v1
  v1 <- calc_quad_fast(Gis2 * j * dMs * l * recD2is, Mis * k, i, D2D3is, idx_g) -
    calc_quad_fast(Gis2 * j * l * recD2is, Mis * k, i, D2D3is_Ms, idx_g)

  # v2
  v2 <- calc_quad_fast(Gis2 * j * dMs * recD2is, Mis^2 * k, i * l, recD3is_Ms, idx_g) -
    calc_quad_fast(Gis2 * j * Mis * recD2is, Mis * k, i * l, recD3is_Ms2, idx_g) -
    calc_quad_fast(Gis2 * j * dMs * Mis * recD2is, Mis * k * dMs, i * l, recD3is, idx_g) +
    calc_quad_fast(Gis2 * j * Mis^2 * recD2is, dMs * k, i * l, recD3is_Ms, idx_g)

  # v3
  v3 <- calc_quad_fast(Gis2 * j * dMs * Mis * recD2is, Mis^2 * k * l, i, recD3is, idx_g) -
    calc_quad_fast(Gis2 * j * Mis^2 * recD2is, Mis * k * l, i, recD3is_Ms, idx_g) -
    calc_quad_fast(Gis2 * j * dMs * recD2is, Mis * k * l, i * dMs, recD3is_Ms, idx_g) +
    calc_quad_fast(Gis2 * j * Mis * recD2is, k * l, i * dMs, recD3is_Ms2, idx_g)

  # v4
  v4 <- t(i * k * dMs) %*% recD2sGs2Pgs %*% (l * j) -
    t(i * k * l) %*% recD2sGs2Pgs %*% (Mis * j)

  return(v1 + v2 + v3 + v4)
}
