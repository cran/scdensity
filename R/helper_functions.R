# FUNCTIONS IN THIS FILE:
#    SequentialLineMin
#    getQuantileFunction
#    getQuantile
#    NormalGridFcn
#    ConvolutionMatrix
#    SpectralShift
#    QPsolve
#    CleanWeights
#    MakeCenters


#' Minimize a function of r variables by sequential univariate searches.
#'
#' The function seeks to minimize \code{fcn}, a scalar function of \eqn{r} variables.  \code{v0}
#' is a starting solution and \code{bounds} is a 2-vector giving upper and lower limits for
#' elements of the solution.
#'
#' This algorithm is designed to search for solutions of the form \eqn{v = [v_1 v_2 \ldots v_r]},
#' where \code{bounds(1)} \eqn{< v_1 < v_2 < ... < v_r <} \code{ bounds(2)}. It loops through the solution vector
#' one variable at a time, and does a 1-D line search using \code{optimize()} for an improving
#' value of that variable.  So when optimizing \eqn{v_i}, it searches the interval \eqn{(v_{i-1},
#' v_{i+1})} to maintain the increasing nature of \eqn{v}. The overall search terminates once a
#' pass through all \eqn{r} elements of \eqn{v} fails to produce any changes to \eqn{v}.
#'
#' @param fcn A function with taking an r-vector as its first argument: call as \code{fcn(v,...)}.
#' @param bounds A 2-vector giving the upper and lower limits for elements of a solution.
#' @param v0 A starting solution, with increasing elements. An r-vector. Not used if r == 1.
#' @param tol Tolerance passed to \code{optimize}.
#'
#' @return a list with elements:
#' \item{\code{minimizer}}{An r-vector containing the solution.}
#' \item{\code{minimum}}{The objective function value at the solution.}
#'
#' @importFrom stats optimize
#' @export
#'
#' @examples
#' fcn <- function(v) (v[1]+1)^2 + (v[2]-1)^2
#' SequentialLineMin(fcn, c(-5,5), c(-3,3))
#'
SequentialLineMin <- function(fcn, bounds, v0, tol=.Machine$double.eps^0.25) {

  v <- v0
  r <- length(v)

  if (r == 1) {
    soln <- optimize(fcn, bounds, tol=tol)
    v <- soln$minimum
    fv <- soln$objective
  } else {
    stop <- FALSE
    while (!stop) {
      vold <- v
      fv <- fcn(vold)
      for (i in 1:r) {
        fcni <- function(u) {
          vv <- v
          vv[i] <- u
          fcn(vv)
        }
        bnd <- c(bounds[1], v, bounds[2])
        soln <- optimize(fcni, c(bnd[i], bnd[i+2]), tol=tol)
        if (soln$objective < fv) {
          v[i] <- soln$minimum
          fv <- soln$objective
        }
      }
      stop <- all(v==vold)
    }
  }

  out <- list()
  out$minimizer <- v
  out$minimum <- fv
  return(out)

}




#' Build a quantile function for a given constrained density estimator.
#'
#' This function implements a crude but numerically reliable method to approximate the
#' quantile function of an \code{scdensity} estimate.
#'
#' @param x An object of S3 class \code{scdensity}.
#' @return A function that takes a fraction and returns the quantile.
#'
#' @keywords internal
getQuantileFunction <- function(x) {

  ng <- length(x$x) - 1
  delta <- diff(range(x$x))/ng
  g <- x$x[1] + delta/2 + delta*(0:(ng-1))
  cdf <- cumsum( delta * (x$y[2:(ng+1)] + x$y[1:ng]) / 2 )
  return(approxfun(cdf, g, ties=mean))

}


#' Estimate a specific quantile of a pdf given abscissa and ordinate values.
#'
#' @param x The abscissa values of the pdf
#' @param y The ordinate values of the pdf
#' @param p The probability at which to evaluate the quantile.
#'
#' @return The estimated quantile.
#'
#' @keywords internal
getQuantile <- function(x, y, p) {
  ng <- length(x) - 1
  delta <- diff(range(x))/ng
  g <- x[1] + delta/2 + delta*(0:(ng-1))
  cdf <- cumsum( delta * (y[2:(ng+1)] + y[1:ng]) / 2 )
  return(approx(cdf, g, p, ties=mean)$y)
}


#' Compute the values of n normal PDFs (or their derivatives at m grid points).
#'
#' Grid points are specified in g. Returns an m-by-n matrix.  The (i,j)th element of
#' the matrix is the rth derivative of N(mu(j),sd(j)) at g(i).
#'
#' @param g  Locations at which to evaluate the normal densities.
#' @param r  Derivative degree of the densities to evaluate (0, 1, 2, or 3).
#' @param mu Means/locations for the normal densities.
#' @param sd Standard deviations of the normal densities.
#'
#' @return  An m-by-n matrix as described above.
#'
#' @keywords internal
NormalGridFcn <- function(g,r,mu,sd){

  #Get dimensions for output (m-by-n). Expand sd if necessary.
  m <- length(g)
  n <- length(mu)
  if (length(sd)==1) {
    sd <- rep(sd,n)
  }

  #Create matrix versions of g, mu, and sd (all m-by-n)
  G <- matrix(g,nrow=m,ncol=n)
  MU <- matrix(mu,nrow=m,ncol=n,byrow=TRUE)
  SIG <- matrix(sd,nrow=m,ncol=n,byrow=TRUE)

  #Calculate the output
  if (r==0) {
    value <- 1/(sqrt(2*pi)*SIG)*exp(-(G-MU)^2/(2*SIG^2))
  } else if (r==1) {
    value <- -(G-MU)/(sqrt(2*pi)*SIG^3)*exp(-((G-MU)^2)/(2*SIG^2))
  } else if (r==2) {
    value <- ((G-MU)^2-SIG^2)/(sqrt(2*pi)*SIG^5)*exp(-((G-MU)^2)/(2*SIG^2))
  } else if (r==3) {
    value <- (3*SIG^2*(G-MU) - (G-MU)^3)/(sqrt(2*pi)*SIG^7)*exp(-((G-MU)^2)/(2*SIG^2));
  } else {
    stop("r must be 0, 1, 2, or 3")
  }

  return(value)

}


#' Computes a matrix of Gaussian kernel convolution values given two vectors.
#'
#' If \code{x} and \code{s} are n- and m-vectors, respectively, returns the n-by-m matrix
#' of convolution values using Gaussian kernel with bandwidth \code{h}.  If \code{x} and
#' \code{s} are equal, the spectral shift is done to ensure the matrix is numerically
#' positive definite.
#'
#' @param x  A numeric vector.
#' @param s  A numeric vector.
#' @param h  A positive bandwidth.
#' @param threshold  Threshold value passed to SpectralShift.
#'
#' @return The matrix of convolution values.
#'
#' @keywords internal
ConvolutionMatrix <- function(x, s, h, threshold = 1e-10){

  convolution <- function(a,b){1/(2*sqrt(pi)) * exp(-((a-b)/h)^2/4)}
  L <- outer(x, s, convolution)
  if (isTRUE(all.equal(x,s))) {
    L <- SpectralShift(L, threshold)
  }
  return(L)
}


#' Performs the spectral shift on a matrix to make it numerically positive definite.
#'
#' Matrix \code{L} is assumed to have eigenvalues that are either all positive, or very
#' close to zero.  If any eigenvalues are less than less than \code{threshold}, a positive
#' quantity is added to the diagonal.
#'
#' @param L  A square numeric matrix.
#' @param threshold  The eigenvalue threshold. Default 1E-10.
#'
#' @return The spectral-shifted matrix.
#'
#' @keywords internal
SpectralShift <- function(L, threshold = 1e-10){
  stopifnot(dim(L)[1] == dim(L)[2])
  minval <- min(eigen(L)$values)
  if (minval < threshold) {
    addAmount <- abs(minval) + threshold
    L <- L + addAmount*diag(dim(L)[1])
  }
  return(L)
}


#' A wrapper to call solve.QP.
#'
#' This function takes arguments slightly differently from solve.QP, to make it more
#' convenient for internal use.  It also implements measures to robustify calls to solve.QP:
#' \itemize{
#' \item A rounding hack to prevent a bug in solve.QP that occasionally produces all-NaN
#'  solutions without returning a warning or error.  Rounding has been found to eliminate
#'  \emph{almost} all such bugs.
#' \item A call to lpSolve's lp() to check feasibility before running solve.QP
#' \item solve.QP is called within tryCatch to eliminate unwanted crashes.
#' }
#' The output of this function is a list with elements
#' \itemize{
#' \item \code{flag} is 0 for successful completion, 1 for failure at the LP check stage,
#'    and 2 for failure at the QP stage (usually the "NaN solution" bug).
#' \item \code{QP} is the list returned by \code{solve.QP}.  If the QP was not run due to
#' infeasibility, this element is NULL.
#' }
#'
#' solve.QP defines its quadratic program as minimizing 1/2 * x'Dx - x'd, subject to
#' constraints A'x >= b.  Equality constraints have to be in the first rows of A'.
#'
#' This function minimizes x'Dx - x'd, subject to inequality constraints Ax >= b and
#' Equality constraints Aeq*x = beq.
#'
#' @param D The matrix of the quadratic objective.
#' @param d The vector in the linear term of the quadratic objective.
#' @param A The matrix of inequality constraints.
#' @param b The vector of RHS of the inequalities.
#' @param Aeq The matrix of equality constraints.
#' @param beq The vector of RHS of the equalities.
#'
#' @return A list with elements described above.
#'
#' @importFrom stats runif
#'
#' @keywords internal
QPsolve <- function(D, d, A, b, Aeq, beq) {

  # Initialize objects.  Round the constraint matrix to ameliorate the solve.QP bug.
  out <- list(flag=0, QP=NULL)
  Afull <- round(rbind(Aeq, A), 8)
  bfull <- c(beq, b)
  meq <- length(beq)

  # Solve a LP to check feasibility.  If infeasible, set flag and stop.
  LPsol <- lpSolve::lp(direction="min",
                       objective.in = rep(1,length(d)),
                       const.mat = Afull,
                       const.dir = c(rep("==",meq), rep(">=", length(b))),
                       const.rhs = bfull)
  if (LPsol$status != 0) {
    #**logging**cat("had an LP goof\n")
    out$flag <- 1
    return(out)
  }

  # If the LP went through, try the QP.
  QPsol <- tryCatch(quadprog::solve.QP(Dmat=2*D, dvec=d, Amat=t(Afull), bvec=bfull,
                                       meq=meq, factorized=FALSE),
                    error = function(e) NULL)
  if (is.null(QPsol) || any(is.na(QPsol$solution))) {

    #***** REMEIDATING SOLVE.QP ***********
    # There is a bug in solve.QP that causes it to occasionally return solutions
    # that are all NA values.  It does this without error or warning. Also it sometimes
    # fails claiming infeasibility, when the problem is actually feasible (as can be shown
    # by solving a linear program with the same constraints).  Both of these problems have
    # been confirmed by Berwin Turlach (quadprog author).  He suggested scaling the
    # constraint matrix sometimes works. So far I've tried two scaling approaches, with
    # some success:
    #   a) Just scale the entire set of inequalities by a positive number.
    #   b) Scale each inequality independently with different random numbers.
    # Either of these options seem to work if you try it with enough random numbers.
    # Here I implement option b) to try to get solutions in cases where the first
    # call to solve.QP fails.
    #
    # Note: Dr. Turlach has indicated that an update to quadprog is in the works.
    #       After it comes out I should review this code and delete it if it's no
    #       longer needed.
    counter <- 0
    repeat {
      counter <- counter + 1
      nr <- dim(A)[1]
      nc <- dim(A)[2]
      vals <- runif(nr, min=0.001, max=2)
      Afull <- round(rbind(Aeq, A*matrix(vals,ncol=nc,nrow=nr)), 8)
      bfull <- c(beq, b*vals)
      QPsol <- tryCatch(quadprog::solve.QP(Dmat=2*D, dvec=d, Amat=t(Afull), bvec=bfull,
                                           meq=meq, factorized=FALSE),
                        error = function(e) NULL)
      if (!is.null(QPsol) && !any(is.na(QPsol$solution))) {
        out$QP <- QPsol
        break
      } else {
        if (counter == 25) {
          out$flag <- 2
          break
        }
      }
    }
    #**** END REMEIDATING SOLVE.QP *************

  } else {
    out$QP <- QPsol
  }
  return(out)

}


#' Clean out near-zeros from a probability weights vector and re-normalize
#'
#' @param w A vector of probability weights
#'
#' @return A vector with the near-zeros removed, still summing to 1.
#'
#' @keywords internal
CleanWeights <- function(w) {
  w <- round(w,6)
  w <- w / sum(w)
}



#' Create a vector of kernel centers covering \[LB, UB\].
#'
#' Make the spacing as large as possible without going over \code{width}.
#' If symmetric, ensure we have an even number of centers that are symmetric around \code{PoS}.
#'
#' @param LB The lower bound.
#' @param UB The upper bound.
#' @param width The maximum possible spacing.
#' @param minspace The minimum allowable spacing.
#' @param PoS Point of symmetry (default NULL)
#'
#' @return A vector of kernel centers.
#'
#' @keywords internal
MakeCenters <- function(LB, UB, width, minspace, PoS=NULL) {
  if (!is.null(PoS)) {
    #For symmetry case, lay out the grid such that the middle interval has minimum spacing.
    #This prevents that interval from ever being bisected.
    M <- ceiling((UB-LB)/width) + 1  #Initial no. of kernel centers.
    M <- M + M %% 2                  #make it even
    sR <- seq(from = PoS + width*minspace/2, to = UB, length.out = M/2)
    s <- sort(c(2*PoS - sR, sR))
  } else {
    #Non-symmetric case, just lay out M initial kernel centers evenly spaced.
    M <- ceiling((UB-LB)/width) + 1
    s <- seq(from=LB, to=UB, length.out=M)
  }
  return(s)
}


