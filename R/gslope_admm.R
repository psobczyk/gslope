#' Graphical SLOPE
#'
#' Method for estimating sparse inverse covariance matrix
#'
#' Implementation using alternating direction method of multipliers
#'
#' @param C sample covariance matrix
#' @param Z initial estimate of inverse covariance
#' @param U initial estimate of covariance
#' @param lambda_seq sequence of lambda values of ordered L1 norm
#' @param rho constant
#' @param max_iter maximum number of ADMM iterations
#' @param tol_infeas tolerance to feasibility rules
#' @param verbose should additional info be printed out
#' @export
gslope_admm <- function(C, Z, U, lambda_seq, rho,
                        max_iter = 100, tol_infeas = 1e-3,
                        verbose = FALSE){

  for(iter in 1:max_iter){
    C_new = rho*(Z - U) - C
    eig <- eigen(C_new, symmetric = TRUE)

    X_new <- eig$vectors %*% diag((eig$values + sqrt(eig$values^2 + 4*rho))/(2*rho)) %*% t(eig$vectors)
    Z_new <- SLOPE::prox_sorted_L1(x = as.vector(X_new + U), lambda = lambda_seq/rho)
    Z_new <- matrix(Z_new, nrow = nrow(U))
    U_new <-  U + (X_new - Z_new)

    dual_feasibility <- norm(as.vector(rho*(Z_new-Z)), type = "2")
    primal_feasibility <- norm(as.vector(Z_new - X_new), type = "2")

    X <- (X_new)
    Z <- (Z_new)
    U <- (U_new)

    if(verbose){
      message(sprintf("Iter %i\nprimal: %f\ndual: %f\n",
                      iter, primal_feasibility, dual_feasibility))
    }


    if(dual_feasibility < tol_infeas & primal_feasibility < tol_infeas){
      break;
    }
  }
  return(list(X = X, Z = Z, U = U,
              primal_feasibility = primal_feasibility,
              dual_feasibility = dual_feasibility))
}
