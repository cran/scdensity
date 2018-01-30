# FUNCTIONS IN THIS FILE:
#    improve
#    isUnimodal
#    isMonotoneR
#    isMonotoneL
#    isBoundedL
#    isBoundedR


#*****************************************************************************************
#*** improve *****************************************************************************
#*****************************************************************************************

#' Move points closer to a target while maintaining a constraint.
#'
#' \code{improve(startValue, x, confun)} uses a greedy algorithm to move the elements of a
#' user-supplied vector \code{startValue} closer to their target values \code{x}, while
#' continually satisfying the constraint-checking function \code{confun}.
#'
#' The algorithm implemented here is the one in Wolters (2012), "A Greedy Algorithm for Unimodal
#' Kernel Density Estimation by Data Sharpening," \emph{Journal of Statistical Software}, 47(6).
#' It could conceivably be useful as a part of other gradient-free optimization schemes where
#' we have an infeasible point and a feasible one, and we seek a point that is on
#' the constraint boundary near the infeasible one.
#'
#' @param startValue The vector of starting values for the search.  Must satisfy
#' \code{confun(startValue) == TRUE}
#' @param x  The target values.
#' @param confun  The constraint-checking function. \code{confun(y)} must return a Boolean value
#' that is invariant to permutations of its vector argument \code{y}.
#' @param verbose  A logical value indicating whether or not information about iteration progress
#' should be printed to the console.
#' @param maxpasses  The maximum allowable number of sweeps through the data points.  At each pass,
#' every point that is not pinned at the constraint boundary is moved toward its target point in a
#' stepping-out procedure.
#' @param tol Numerical tolerance for constraint checking.  A point is considered to be at the
#' constraint boundary if adding \code{tol} to it causes the constraint to be violated. If \code{tol}
#' is too large, the algorithm will terminate prematurely.  If it is too small, run time will be
#' increased with no discernible benefit in the result.
#'
#' @return A vector of the same length as \code{startValue}, with elements closer to \code{x}.
#' @export
#'
#' @examples
#' #Constrain points to be inside the hypercube with vertices at -1 and +1.
#' #The target point is a vector of independent random standard normal variates.
#' #Start at rep(0,n) and "improve" the solution toward the target.
#' n <- 20
#' incube <- function(x) all(x <= 1 & x >= -1)
#' x0 <- rep(0,n)
#' target <- sort(rnorm(n))
#' xstar <- improve(x0, target, incube, verbose=TRUE)
#' dist <- abs(target - xstar)
#' zapsmall(cbind(target, xstar, dist), 4)

improve<-function(startValue, x, confun, verbose=FALSE, maxpasses=500,
                  tol=diff(range(c(startValue, x)) / 1e5)){

  #=== PRELIMINARIES ================================================

  #--- Input checking ---------------------------
  stopifnot(
    is.vector(startValue, mode='numeric'), !any(is.na(startValue)),
    is.vector(x, mode='numeric'), !any(is.na(x)),
    length(startValue)==length(x),
    length(maxpasses)==1, maxpasses>=1,
    length(tol)==1 && tol>0,
    is.function(confun)
  )
  if (!confun(startValue))
    stop("The initial solution does not satisfy the constraint.")

  #--- Initialize objects -----------------------
  # During the search, we want y to be matched to x such that the jth smallest value of y
  # has the jth smallest value of x as its target.  But we also must prevent cycling. Work
  # with x sorted in ascending order.  Then order(y) gives the sort order to bring y into
  # a matched state.  And actually doing the matching is just a sort(y) operation. Save
  # rank(x) for putting the final y back into matching state with the original x.
  n <- length(startValue)
  idx <- 1:n
  count <- 1
  ranks_of_x <- rank(x, ties.method="first")
  x <- sort(x)
  y <- sort(startValue)
  ordermem <- matrix(0, nrow=maxpasses+1, ncol=n)
  ordermem[1, ] <- order(y)
  unitvec <- function(a,b) (b-a)/sqrt(sum((b-a)^2))  #-unit vector from a to b.

  #--- Find the set of moveable points ----------
  D <- abs(x-y)                        #-Distances between elements of x and y.
  S <- sign(x-y)                       #-Sign of move from y to x.
  home <- D<=tol
  nothome <- idx[!home]
  pinned <- rep(FALSE,n)
  for (i in 1:length(nothome)) {
    tst <- y;
    j <- nothome[i]
    tst[j] <- tst[j] + tol*S[j]
    if (!confun(tst)) {
      pinned[j] <- TRUE
    }
  }
  moveable <- !(home | pinned)
  m <- sum(moveable)

  #=== RUN THE ALGORITHM ============================================
  # Perform repeated passes through the data points.  At each pass, try to move each point
  # toward home by stepping out from its current location.  The step size is a proportion
  # of the distance |y(j) - x(j)|; this proportion starts at 1 (step all the way) and is
  # reduced by factors of 2 whenever no moves can be made.

  #--- Set the sweep order ----------------------
  # Sweep in descending order of distance from home.
  pts <- idx[moveable]
  ix <- order(D[moveable], decreasing = TRUE)
  pts <- pts[ix]                       #-The ordered set of moveables.

  #--- SWEEP UNTIL NO MOVEABLE POINTS -----------
  steps <- 1
  while (m>0 && count<maxpasses) {

    nmoved <- 0                        #-move counter.

    #--- A SINGLE PASS THROUGH THE POINTS -------
    for (i in 1:m) {
      j <- pts[i]
      A <- y[j]                        #-The current moveable point.
      B <- x[j]                        #-The target of the current point.
      D <- abs(B-A)                    #-Distance between the points.
      S <- sign(B-A)                   #-Sign of the move from A to B.
      k <- 0                           #-Counter for number of successful steps.
      stepsize <- max(D/steps, tol)    #-Min poss step size is equal to tol.
      stop <- FALSE
      while (!stop && (k+1)*stepsize<=D) {
        # Start stepping out.  Stop once constraint is violated or B is reached.
        y[j] <- A + (k+1)*stepsize*S
        if (confun(y)) {
          k <- k + 1;
        } else {
          stop <- TRUE
        }
      }
      y[j] <- A + k*stepsize*S         #-Set y[j] to last feasible step.
      if (k>0) {
        nmoved <- nmoved + 1
      }
    }

    #--- Show progress if requested -------------
    if (verbose) {
      cat('Sweep:', count, '  moves:', nmoved, '  moveable:', m, '  steps:', steps, '\n')
    }

    #--- Prepare for the next sweep -------------
    # Increment the sweep counter and the order memory. If points have been moved,
    # re-match the points as necessary.  If not, halve the step size.
    count <- count + 1

    if (nmoved > 0) {
      #--- Re-sort points -----------------------
      thisorder <- order(y)
      hasmatch <- any(apply(ordermem, 1, function(v,w) all(v==w), thisorder))
      if (!hasmatch) {
        y <- sort(y)
        ordermem[count, ] <- thisorder
      }
      #--- Re-find the moveable points ----------
      D <- abs(x-y)
      S <- sign(x-y)
      home <- D<=tol
      nothome <- idx[!home]
      pinned <- rep(FALSE,n)
      for (i in 1:length(nothome)) {
        tst <- y;
        j <- nothome[i]
        tst[j] <- tst[j] + tol*S[j]
        if (!confun(tst)) {
          pinned[j] <- TRUE
        }
      }
      moveable <- !(home | pinned)
      m <- sum(moveable)

      #--- Re-set sweep order -------------------
      pts <- idx[moveable]
      ix <- order(D[moveable], decreasing = TRUE)
      pts <- pts[ix]
    } else {
      steps <- 2*steps
    }

  }

  #=== RETURN OUTPUTS ===============================================
  if (m>0) {
    warning("Search stopped because maxpasses was reached.")
  }
  y[ranks_of_x]

}


#*****************************************************************************************
#*** isUnimodal **************************************************************************
#*****************************************************************************************

#' Check for unimodality of function values.
#'
#' Given a set of function values for increasing abscissa values, we call this unimodal if
#' there are zero or one values that are greater than all of their neighbors. Before
#' checking for modes, the values are scaled to fill [0, 1] and then rounded to four
#' decimal places.  This eliminates unwanted detection of tiny differences as modes.
#'
#' This function is intended to be called from other functions in the scdensity package.
#' It does not implement any argument checking.
#'
#' @param f  A vector of function values for increasing abscissa values.
#'
#' @return A logical value indicating if unimodality is satisfied.
#'
#' @keywords internal
isUnimodal = function(f)  {

  f <- (f - min(f))/(max(f) - min(f))
  f <- round(f, 4)

  d = diff(f)
  d = d[d!=0];
  chgs = sum(diff(sign(d))!=0);       #-Number of sign changes.

  chgs<=1;                            #-Chgs could be 0 for monotone, 1 for unimodal.

}


#*****************************************************************************************
#*** isMonotoneR *************************************************************************
#*****************************************************************************************

#' Check for monotonicity of function values in the right tail.
#'
#' Given a vector of n function values and an index ix, determines whether the function
#' values having indices greater than or equal to ix are non-decreasing or non-increasing.
#' Returns TRUE if they are, FALSE otherwise.
#'
#' As in \code{isUnimodal}, the values are first scaled to fill [0, 1] and then rounded to
#' four decimal places.  This eliminates unwanted detection of tiny differences as modes.
#'
#' This function is intended to be called from other functions in the scdensity package.
#' It does not implement any argument checking.
#'
#' @param f  A vector of function values for increasing abscissa values.
#' @param ix An index giving the cutoff for checking monotonicity.
#'
#' @return A logical value indicating if the constraint is satisfied.
#'
#' @keywords internal
isMonotoneR = function(f, ix)  {

  f <- (f - min(f))/(max(f) - min(f))
  f <- round(f, 4)

  n <- length(f)
  d = diff(f[ix:n])
  d = d[d!=0];
  chgs = sum(diff(sign(d))!=0);       #-Number of sign changes.

  chgs<=0;                            #-Chgs must be 0 for monotone.

}


#*****************************************************************************************
#*** isMonotoneL *************************************************************************
#*****************************************************************************************

#' Check for monotonicity of function values in the left tail.
#'
#' Given a vector of n function values and an index ix, determines whether the function
#' values having indices less than or equal to ix are non-decreasing or non-increasing.
#' Returns TRUE if they are, FALSE otherwise.
#'
#' As in \code{isUnimodal}, the values are first scaled to fill [0, 1] and then rounded to
#' four decimal places.  This eliminates unwanted detection of tiny differences as modes.
#'
#' This function is intended to be called from other functions in the scdensity package.
#' It does not implement any argument checking.
#'
#' @param f  A vector of function values for increasing abscissa values.
#' @param ix An index giving the cutoff for checking monotonicity.
#'
#' @return A logical value indicating if the constraint is satisfied.
#'
#' @keywords internal
isMonotoneL = function(f, ix)  {

  f <- (f - min(f))/(max(f) - min(f))
  f <- round(f, 4)

  n <- length(f)
  d = diff(f[1:ix])
  d = d[d!=0];
  chgs = sum(diff(sign(d))!=0);       #-Number of sign changes.

  chgs<=0;                            #-Chgs must be 0 for monotone.

}


#*****************************************************************************************
#*** isBoundedL *************************************************************************
#*****************************************************************************************

#' Check for zeros at the left side of a vector of function values.
#'
#' Given a vector of n function values and an index ix, check whether values 1:ix are zero.
#' Returns TRUE if they are, FALSE otherwise.
#'
#' As in \code{isUnimodal}, the values are first scaled to fill [0, 1] and then rounded to
#' four decimal places.  Because of this it is still possible to use the "bounded support"
#' constraints with the Gaussian kernel.
#'
#' This function is intended to be called from other functions in the scdensity package.
#' It does not implement any argument checking.
#'
#' @param f  A vector of function values for increasing abscissa values.
#' @param ix An index giving the cutoff for checking for zero.
#'
#' @return A logical value indicating if the constraint is satisfied.
#'
#' @keywords internal
isBoundedL = function(f, ix)  {

  f <- (f - min(f))/(max(f) - min(f))
  f <- round(f, 4)

  all(f[1:ix]==0)

}


#*****************************************************************************************
#*** isBoundedR *************************************************************************
#*****************************************************************************************

#' Check for zeros at the right side of a vector of function values.
#'
#' Given a vector of n function values and an index ix, check whether values with indices
#' greater than or equal to ix are zero. Returns TRUE if they are, FALSE otherwise.
#'
#' As in \code{isUnimodal}, the values are first scaled to fill [0, 1] and then rounded to
#' four decimal places.  Because of this it is still possible to use the "bounded support"
#' constraints with the Gaussian kernel.
#'
#' This function is intended to be called from other functions in the scdensity package.
#' It does not implement any argument checking.
#'
#' @param f  A vector of function values for increasing abscissa values.
#' @param ix An index giving the cutoff for checking for zero.
#'
#' @return A logical value indicating if the constraint is satisfied.
#'
#' @keywords internal
isBoundedR = function(f, ix)  {

  f <- (f - min(f))/(max(f) - min(f))
  f <- round(f, 4)
  n <- length(f)

  all(f[ix:n]==0)

}





