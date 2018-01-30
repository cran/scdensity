# FUNCTIONS IN THIS FILE:
#    InitializeP
#    WeightedKDE
#    BinningStep
#    BuildConstraints
#    EstimationStep

# For maximum flexibility, use a list called P (for "Problem") to hold all objects related
# to the estimation process.  Then pass it around among functions BinningStep, BuildConstraints,
# and EstimationStep that mutate the list to carry out the optimization.


#' Initialize the list of problem-related objects.
#'
#' @param x The data vector.
#' @param h The bandwidth (a positive scalar).
#' @param constraints The vector of constraint strings.
#' @param method Either "weightedKDE" or "adjustedKDE".
#'
#' @return A list in the form needed by WeightedKDE(), with additional elements initialized
#'  to their default values.
#' @keywords internal
InitializeP <- function(x, h, constraints, method, opts) {

  stopifnot(length(h)==1, h > 0,
            is.numeric(x), all(is.finite(x)),
            method %in% c("weightedKDE", "adjustedKDE")
            )

  if ("ncheck" %in% names(opts)) {
    ng <- opts$ncheck
  } else {
    ng <- NULL
  }

  if ("verbose" %in% names(opts)) {
    if (isTRUE(opts$verbose)) {
      verbose <- TRUE
      cat("\n")
    } else {
      verbose <- FALSE
    }
  } else {
    verbose <- FALSE
  }

  list(x = x,
       h = h,
       constraints = constraints,
       method = method,
       minSpace = 0.1,      #Minimum spacing of kernel centers for binning step.
       approxTol = 0.001,   #Tolerance for approximating the unconstrained estimator in binning.
       ng = ng,             #Size of constraint-checking grid.
       PoS = NULL,          #Point of symmetry (only for symmetric constraint)
       modeLoc=NULL,        #Mode location (only for unimodality & bimodality constraints)
       rightPoint=NULL,     #Used for monotoneRightTail constraint.
       leftPoint=NULL,      #Used for monotoneLeftTail constraint.
       pts=NULL,            #Used for constraints on inflection points.
       lowerBound=NULL,     #Used for boundedLeft constraint.
       upperBound=NULL,     #Used for boundedRight constraint.
       boundTol=0.001/diff(range(x)),   #Tolerance for the bounded support constraints.
       LB = NULL,   #Lower bound for g.
       UB = NULL,   #Upper bound for g.
       g = NULL,    #Constraint-checking grid.
       s = NULL,    #Kernel centers after binning (including zero-weighted ones).
       y = NULL,    #Kernel centers passed to estimation step.
       w = NULL,    #Optimal unconstrained weights corresponding to s.
       v = NULL,    #Optimal unconstrained weights corresponding to y.
       Apos = NULL, #Matrix of inequality constraints (positivity).
       Aunit = NULL,  #Matrix of equality constraints (unit integral).
       Aeq = NULL,  #Matrix of equality constraints (derivative(s) at important points).
       Ashape = NULL, #Matrix of inequality constraints (shape restrictions).
       bshape = NULL, #Constant vector for inequalities using Ashape.
       Lbin = NULL,  #Matrix in QP for binning step.
       dbin = NULL,  #vector in QP for binning step.
       Lest = NULL,  #Matrix in QP for estimation step.
       dest = NULL,  #vector in QP for estimation step.
       vhat = NULL,  #Optimal weights for y (shape constrained).
       flag = 0,      #Status flag to be filled
       verbose = verbose  #Flag indicating whether to show progress at the console.
       )
}




#' Function to carry out the weighted or adjusted KDE optimization.
#'
#' This function sets up the problem and finds an optimal shape-constrained estimate for a
#' specified set of important points.
#'
#' @param P A list, as created by InitializeP().
#'
#' @return A mutated version of the input list, with additional elements giving the optimization
#'         output.
#' @keywords internal
WeightedKDE <- function(P) {

  symmetric <- "symmetric" %in% P$constraints

  # If symmetry required, or if it hasn't been done already, set up the constraint-checking grid.
  if (symmetric || is.null(P$g)) {
    P <- BuildConCheckGrid(P)
  }

  # If symmetry required, or if it hasn't been done already, do the binning step.
  if (symmetric || is.null(P$s)) {
    P <- BinningStep(P)
  }

  # Formulate the shape constraints
  P <- BuildConstraints(P)

  # Get the optimal weights
  P <- EstimationStep(P)

  # Return output
  return(P)
}


#' A function to add the constraint-checking grid
#'
#' @param P List with problem details.
#'
#' @return Modified list.
#' @keywords internal
BuildConCheckGrid <- function(P) {
  P$LB <- min(P$x) - 3*P$h
  P$UB <- max(P$x) + 3*P$h
  if(is.null(P$ng)) {
    P$ng <- max(100, ceiling((P$UB-P$LB)/P$h))
  }
  if ("symmetric" %in% P$constraints) {
    halfWidth <- max(abs(c(P$LB, P$UB) - P$PoS))
    P$LB <- P$PoS - halfWidth
    P$UB <- P$PoS + halfWidth
    P$ng <- P$ng + P$ng %% 2   #For symmetry, need even no. of check points.
  }
  P$g <- seq(P$LB - 2*P$h, P$UB + 2*P$h, length.out=P$ng)
  return(P)
}



#' A small function to add the objects needed for the binning QP to the problem list.
#'
#' @param P A list with the problem details.
#' @param s A vector of centers.
#'
#' @return A modified problem list.
#' @keywords internal
AddBinObjects <- function(P,s) {
  nx <- length(P$x)
  M <- length(s)
  P$Lbin <- nx * ConvolutionMatrix(s, s, P$h)
  P$dbin <- 2 * ConvolutionMatrix(s, P$x, P$h) %*% rep(1, nx)
  P$Apos <- diag(M)
  P$Aunit <- matrix(1, nrow=1, ncol=M)
  return(P)
}


#' Carry out the binning step.
#'
#' @param P The list of problem objects.
#'
#' @return The input list, with extra members modified.
#'
#' @importFrom stats density
#' @keywords internal
BinningStep <- function(P) {

  symmetric <- "symmetric" %in% P$constraints

  # Create the initial s vector.
  nx <- length(P$x)
  s <- MakeCenters(P$LB, P$UB, P$h, P$minSpace, P$PoS)
  M <- length(s)

  # Loop to bisect intervals until stopping criteria are met.
  epsilon <- P$approxTol * max(density(P$x,bw=P$h)$y)
  stop = FALSE
  while(!stop) {

    # Set up the QP and solve.
    P <- AddBinObjects(P, s)
    QPsol <- QPsolve(D=P$Lbin, d=P$dbin, A=P$Apos, b=rep(0,M), Aeq=P$Aunit, beq=1)
    if (QPsol$flag != 0) {
      # If had QP problems, just try a fixed-width grid with width about h/2.
      s <- MakeCenters(P$LB, P$UB, P$h/2, P$minSpace, P$PoS)
      QPsol <- QPsolve(D=P$Lbin, d=P$dbin, A=P$Apos, b=rep(0,M), Aeq=P$Aunit, beq=1)
      if (QPsol$flag == 0) {
        w <- QPsol$QP$solution
      } else {
        # If that didn't work, fall back to using x with uniform weights (with warning).
        # For symmetry case, throw an error.
        if (symmetric) {
          stop(paste0("The QP solver had numerical problems at the binning step.\n",
                      "Try changing the bandwidth or the number of constraint-checking points."))
        } else {
          s <- P$x
          w <- rep(1/nx, nx)
          warning(paste0("The QP solver had numerical problems at the binning step.\n",
                         "Using the full data vector as kernel centers."))
        }
      }
      break
    }
    w <- QPsol$QP$solution
    P$binSoln <- QPsol$QP

    # Identify the intervals to bisect. If symmetric, bisect on both sides.
    okgaps <- (s[2:M] - s[1:(M-1)]) >= 2*P$h*P$minSpace  #-Intervals wide enough for bisection.
    midpts <- (s[2:M] + s[1:(M-1)])/2
    f0 <- NormalGridFcn(midpts, 0, P$x, P$h) %*% rep(1/nx, nx)
    fhat <- NormalGridFcn(midpts, 0, s, P$h) %*% w
    errs <- (abs(f0 - fhat) > epsilon)
    pick <- okgaps & errs
    if (symmetric) {
      pickFlipped <- M - which(pick)
      pick[pickFlipped] <- TRUE
    }

    # If any intervals need bisection, do it. Otherwise stop.
    if (any(pick)) {
      s <- sort(c(s, midpts[pick]))
      M <- length(s)
    } else {
      stop <- TRUE
    }
  }

  # Add the results to P and return
  if (symmetric) {
    sR <- s[s > P$PoS]
    sL <- 2*P$PoS - sR
    s <- c(sL, sR)
    w[1:length(sL)] <- w[length(sL):1]
  }
  P$s <- s
  P$w <- CleanWeights(w)
  if (P$method == "adjustedKDE") {
    P$y <- P$s
    P$v <- P$w
  } else {
    P$y <- s[P$w>0]
    P$v <- w[P$w>0]
  }

  if (P$verbose) {
    cat("Binning    | completed with ", length(P$s), " kernel centers, ",
        sum(P$w>0), " nonzero weights.\n")
  }

  return(P)

}



#' Build objects Apos, Aunit, Aeq, Ashape, bshape
#'
#' This function builds the matrices/vectors needed to implement shape constraints in the
#' estimation step.
#'
#' @param P A list of problem obejcts
#'
#' @return The same list, with Apos, Aunit, Aeq, Ashape, bshape filled in.
#' @keywords internal
BuildConstraints <- function(P) {

  symmetric <- "symmetric" %in% P$constraints
  ny <- length(P$y)
  P$bshape <- numeric(0)
  if (symmetric) {
    yL <- P$y[1:(ny/2)]
    yR <- P$y[(ny/2 + 1):ny]
    P$Ashape <- matrix(nrow=0, ncol=ny/2)
    P$Aeq <- matrix(nrow=0, ncol=ny/2)
  } else {
    P$Ashape <- matrix(nrow=0, ncol=ny)
    P$Aeq <- matrix(nrow=0, ncol=ny)
  }

  # Sum constraint (equality constraint)
  if (!symmetric) {
    P$Aunit <- matrix(1, nrow=1, ncol=ny)
  } else {
    P$Aunit <- matrix(1, nrow=1, ncol=ny/2)
  }

  # Nonnegativity constraint (inequality)
  if (!symmetric) {
    P$Apos <- diag(ny)
  } else {
    P$Apos <- diag(ny/2)
  }

  # Unimodality constraint: requires P$modeLoc to be non-null
  if ("unimodal" %in% P$constraints) {
    if (!symmetric) {
      AeqNew <- NormalGridFcn(P$modeLoc, 1, P$y, P$h)
      AshapeNew <- NormalGridFcn(P$g, 1, P$y, P$h)
      pick <-  P$g >= P$modeLoc
      AshapeNew[pick, ] <- -AshapeNew[pick, ]
    } else {
      if (!isTRUE(all.equal(P$modeLoc, P$PoS))) {
        stop("Mode location does not equal point of symmetry.")
      }
      gR <- P$g[(P$ng/2+1):P$ng]
      AeqNew <- NormalGridFcn(P$modeLoc, 1, yL, P$h) + NormalGridFcn(P$modeLoc, 1, yR, P$h)
      AshapeNew <- -(NormalGridFcn(gR, 1, yL, P$h) + NormalGridFcn(gR, 1, yR, P$h))
    }
    P$Aeq <- rbind(P$Aeq, AeqNew)
    P$Ashape <- rbind(P$Ashape, AshapeNew)
    P$bshape <- c(P$bshape, rep(0, dim(AshapeNew)[1]))
  }

  # Monotone right tail constraint: requires rightPoint to be non-null
  if ("monotoneRightTail" %in% P$constraints) {
    gR <- P$g[P$g >= P$rightPoint]
    if (!symmetric) {
      AshapeNew <- NormalGridFcn(gR, 1, P$y, P$h)
    } else {
      if (P$rightPoint < P$PoS) {
        stop("Right tail cannot begin on the left side of the point of symmetry.")
      }
      AshapeNew <- NormalGridFcn(gR, 1, yL, P$h) + NormalGridFcn(gR, 1, yR, P$h)
    }
    P$Ashape <- rbind(P$Ashape, -AshapeNew)
    P$bshape <- c(P$bshape, rep(0, length(gR)))
  }

  # Monotone left tail constraint: requires leftPoint to be non-null
  if ("monotoneLeftTail" %in% P$constraints) {
    gL <- P$g[P$g <= P$leftPoint]
    if (!symmetric) {
      AshapeNew <- NormalGridFcn(gL, 1, P$y, P$h)
    } else {
      if (P$leftPoint > P$PoS) {
        stop("Left tail cannot begin on the right side of the point of symmetry.")
      }
      AshapeNew <- NormalGridFcn(gL, 1, yL, P$h) + NormalGridFcn(gL, 1, yR, P$h)
    }
    P$Ashape <- rbind(P$Ashape, AshapeNew)
    P$bshape <- c(P$bshape, rep(0, length(gL)))
  }

  # Two inflections constraint: requires pts to be non-null
  if ("twoInflections" %in% P$constraints) {
    if (!symmetric) {
      AeqNew <- NormalGridFcn(P$pts, 2, P$y, P$h)
      AshapeNew <- NormalGridFcn(P$g, 2, P$y, P$h)
      pick <- P$g > P$pts[1] & P$g < P$pts[2]
      AshapeNew[pick, ] <- -AshapeNew[pick, ]
    } else {
      if (!isTRUE(all.equal(P$pts[1] + P$pts[2], 2*P$PoS))) {
        stop("Inflection points are not equidistant from the point of symmetry.")
      }
      if (P$pts[2] < P$PoS) {
        stop("Inflection points cannot both be on the same side of the point of symmetry.")
      }
      gR <- P$g[P$g >= P$PoS]
      AeqNew <- NormalGridFcn(P$pts, 2, yL, P$h) + NormalGridFcn(P$pts, 2, yR, P$h)
      AshapeNew <- NormalGridFcn(gR, 2, yL, P$h) + NormalGridFcn(gR, 2, yR, P$h)
      pick <- gR < P$pts[2]
      AshapeNew[pick, ] <- -AshapeNew[pick, ]
    }
    P$Aeq <- rbind(P$Aeq, AeqNew)
    P$Ashape <- rbind(P$Ashape, AshapeNew)
    P$bshape <- c(P$bshape, rep(0, dim(AshapeNew)[1]))
  }

  # Two inflections "plus" constraint: requires pts to be non-null
  if ("twoInflections+" %in% P$constraints) {
    if (!symmetric) {
      AeqNew <- NormalGridFcn(P$pts, 3, P$y, P$h)
      AshapeNew <- NormalGridFcn(P$g, 3, P$y, P$h)
      pick1 <- P$g >= P$pts[1] & P$g <= P$pts[2]
      pick2 <- P$g >= P$pts[3]
      AshapeNew[pick1|pick2, ] <- -AshapeNew[pick1|pick2, ]
    } else {
      if (!isTRUE(all.equal(P$pts[2], P$PoS))) {
        stop("The middle inflection point must equal the point of symmetry.")
      }
      if (!isTRUE(all.equal(P$pts[1] + P$pts[3], 2*P$PoS))) {
        stop("Inflection points are not symmetric around the point of symmetry.")
      }
      if (P$pts[3] < P$PoS) {
        stop("Inflection points cannot all be on the same side of the point of symmetry.")
      }
      gR <- P$g[(P$ng/2+1):P$ng]
      AeqNew <- NormalGridFcn(P$pts, 3, yL, P$h) + NormalGridFcn(P$pts, 3, yR, P$h)
      AshapeNew <- NormalGridFcn(gR, 3, yL, P$h) + NormalGridFcn(gR, 3, yR, P$h)
      pick <- gR > P$pts[3]
      AshapeNew[pick, ] <- -AshapeNew[pick, ]
    }
    P$Aeq <- rbind(P$Aeq, AeqNew)
    P$Ashape <- rbind(P$Ashape, AshapeNew)
    P$bshape <- c(P$bshape, rep(0, dim(AshapeNew)[1]))
  }

  # Bounded left constraint: requires lowerBound and boundTol to be non-null
  if ("boundedLeft" %in% P$constraints) {
    gL <- P$g[P$g <= P$lowerBound]
    if(!symmetric) {
      AshapeNew <- NormalGridFcn(gL, 0, P$y, P$h)
    } else {
      if (P$lowerBound > P$PoS) {
        stop("The lower bound cannot be greater than the point of symmetry.")
      }
      AshapeNew <- NormalGridFcn(gL, 0, yL, P$h) + NormalGridFcn(gL, 0, yR, P$h)
    }
    P$Ashape <- rbind(P$Ashape, -AshapeNew)
    P$bshape <- c(P$bshape, rep(-P$boundTol, dim(AshapeNew)[1]))
  }

  # Bounded right constraint: requires upperBound and boundTol to be non-null
  if ("boundedRight" %in% P$constraints) {
    gR <- P$g[P$g >= P$upperBound]
    if (!symmetric) {
      AshapeNew <- NormalGridFcn(gR, 0, P$y, P$h)
    } else {
      if (P$upperBound < P$PoS) {
        stop("The upper bound cannot be less than the point of symmetry.")
      }
      if ("boundedLeft" %in% P$constraints) {
        if (!isTRUE(all.equal(P$lowerBound + P$upperBound, 2*P$PoS))) {
          stop("Lower and upper bounds must be equidistant from the point of symmetry.")
        }
      } else {
        AshapeNew <- NormalGridFcn(gR, 0, yL, P$h) + NormalGridFcn(gR, 0, yR, P$h)
        P$Ashape <- rbind(P$Ashape, -AshapeNew)
        P$bshape <- c(P$bshape, rep(-P$boundTol, dim(AshapeNew)[1]))
      }
    }
  }

  # Bimodal constraint: requires modeLoc to be non-null
  if ("bimodal" %in% P$constraints) {
    if (!symmetric) {
      AeqNew <- NormalGridFcn(P$modeLoc, 1, P$y, P$h)
      AshapeNew <- NormalGridFcn(P$g, 1, P$y, P$h)
      pick1 <- P$g >= P$modeLoc[1] & P$g <= P$modeLoc[2]
      pick2 <- P$g >= P$modeLoc[3]
      AshapeNew[pick1|pick2, ] <- -AshapeNew[pick1|pick2, ]
    } else {
      if (!isTRUE(all.equal(P$modeLoc[2], P$PoS))) {
        stop("The antimode must be located at the point of symmetry.")
      }
      if (!isTRUE(all.equal(P$modeLoc[1] + P$modeLoc[3], 2*P$PoS))) {
        stop("Mode locations must be eqidistant from the point of symmetry.")
      }
      if (P$modeLoc[3] < P$PoS) {
        stop("Both modes can not be on the same side of the point of symmetry.")
      }
      gR <- P$g[(P$ng/2+1):P$ng]
      AeqNew <- NormalGridFcn(P$modeLoc, 1, yL, P$h) + NormalGridFcn(P$modeLoc, 1, yR, P$h)
      AshapeNew <- NormalGridFcn(gR, 1, yL, P$h) + NormalGridFcn(gR, 1, yR, P$h)
      pick <- gR >= P$modeLoc[3]
      AshapeNew[pick, ] <- -AshapeNew[pick, ]
    }
    P$Aeq <- rbind(P$Aeq, AeqNew)
    P$Ashape <- rbind(P$Ashape, AshapeNew)
    P$bshape <- c(P$bshape, rep(0, dim(AshapeNew)[1]))
  }

  return(P)

}



#' Carry out the shape-constrained estimation
#'
#' @param P The problem list
#'
#' @return Problem list with estimate objects.
#' @keywords internal
EstimationStep <- function(P) {

  ny <- length(P$y)
  if("symmetric" %in% P$constraints) {
    yL <- P$y[1:(ny/2)]
    yR <- P$y[(ny/2 + 1):ny]
    L_RR <- ConvolutionMatrix(yR, yR, P$h)
    L_LR <- ConvolutionMatrix(yL, yR, P$h)
    L_Ly <- ConvolutionMatrix(yL, P$y, P$h)
    L_Ry <- ConvolutionMatrix(yR, P$y, P$h)
    P$Lest <- L_RR + L_LR
    P$dest <- (L_Ly +L_Ry) %*% P$v
    QPsol <- QPsolve(D=P$Lest, d=P$dest,
                     A=rbind(P$Apos, P$Ashape), b=c(rep(0, ny/2), P$bshape),
                     Aeq=rbind(P$Aunit, P$Aeq), beq=c(1/2, rep(0, dim(P$Aeq)[1])) )
    if (QPsol$flag != 0) {
      P$estSoln <- list(solution=NULL, value=.Machine$double.xmax)
      P$vhat <- rep(1/ny, ny)
      P$flag <- QPsol$flag
    } else {
      P$estSoln <- QPsol$QP
      P$vhat <- CleanWeights(rep(P$estSoln$solution, 2))
      P$flag <- 0
    }
  } else {
    Lyy <- ConvolutionMatrix(P$y, P$y, P$h)
    P$Lest <- Lyy
    P$dest <- 2 * Lyy %*% P$v
    QPsol <- QPsolve(D=P$Lest, d=P$dest,
                     A=rbind(P$Apos, P$Ashape), b=c(rep(0, ny), P$bshape),
                     Aeq=rbind(P$Aunit, P$Aeq), beq=c(1, rep(0, dim(P$Aeq)[1])) )
    if (QPsol$flag != 0) {
      P$estSoln <- list(solution=NULL, value=.Machine$double.xmax)
      P$vhat <- rep(1/ny, ny)
      P$flag <- QPsol$flag
    } else {
      P$estSoln <- QPsol$QP
      P$vhat <- CleanWeights(P$estSoln$solution)
      P$flag <- 0
    }
  }
  if (P$verbose) {
    cat("Estimation | ")
  }

  return(P)

}

#' Display estimation results to console in a successful case.
#'
#' @param pts The important points.
#' @param value The objective function value.
#'
#' @keywords internal
displaySuccess <- function(pts, value) {
  cat("important points: ", pts, "  objective value: ", value, "\n")
}

#' Display estimation results to console in an unsuccessful case.
#'
#' @param pts The important points.
#'
#' @keywords internal
displayFailure <- function(pts) {
  cat("important points: ", pts, "  optimization FAILED. Continuing search.\n")
}

#' A function factory for making the search objective function.
#'
#' Used when we need to search for important points.
#' P is the problem list. It should have already gone through BuildConCheckGrid and
#' BinningStep. The returned function must return a value even if WeightedKDE() fails.
#' In case of failure, just assign a large random value to the objective value (to keep
#' the search from stagnating or moving systematically in one direction).
#'
#' @param P The list of problem details.
#'
#' @return The objective function.
#' @importFrom stats runif
#' @keywords internal
makeOF <- function(P) {

  fcn <- function(theta) {}  #-default return value
  symmetric <- "symmetric" %in% P$constraints

  FailureValue <- function(x) 1000 + 10*runif(1)

  # Make the function (unimodality case)
  if ("unimodal" %in% P$constraints) {
    fcn <- function(theta) {
      P$modeLoc <- theta
      if (symmetric) P$PoS <- theta
      Pout <- WeightedKDE(P)
      if (Pout$flag == 0) {
        if (P$verbose) displaySuccess(theta, Pout$estSoln$value)
        return(Pout$estSoln$value)
      } else {
        if (P$verbose) displayFailure(theta)
        return(FailureValue(theta))
      }
    }
  }

  # Make the function (2IP case)
  if ("twoInflections" %in% P$constraints) {
    if (!symmetric) {
      fcn <- function(theta) {
        P$pts <- theta
        Pout <- WeightedKDE(P)
        if (Pout$flag == 0) {
          if (P$verbose) displaySuccess(theta, Pout$estSoln$value)
          return(Pout$estSoln$value)
        } else {
          if (P$verbose) displayFailure(theta)
          return(FailureValue(theta))
        }
      }
    }
    if (symmetric && !is.null(P$PoS)) {
      fcn <- function(theta) {
        P$pts <- c(2*P$PoS - theta, theta)
        Pout <- WeightedKDE(P)
        if (Pout$flag == 0) {
          if (P$verbose) displaySuccess(theta, Pout$estSoln$value)
          return(Pout$estSoln$value)
        } else {
          if (P$verbose) displayFailure(theta)
          return(FailureValue(theta))
        }
      }
    }
    if (symmetric && is.null(P$PoS)) {
      fcn <- function(theta) {
        P$PoS <- theta[1]
        P$pts <- c(2*P$PoS - theta[2], theta[2])
        Pout <- WeightedKDE(P)
        if (Pout$flag == 0) {
          if (P$verbose) displaySuccess(theta, Pout$estSoln$value)
          return(Pout$estSoln$value)
        } else {
          if (P$verbose) displayFailure(theta)
          return(FailureValue(theta))
        }
      }
    }
  }

  # Make the function (2IP+ case)
  if ("twoInflections+" %in% P$constraints) {
    if (!symmetric) {
      fcn <- function(theta) {
        P$pts <- theta
        Pout <- WeightedKDE(P)
        if (Pout$flag == 0) {
          if (P$verbose) displaySuccess(theta, Pout$estSoln$value)
          return(Pout$estSoln$value)
        } else {
          if (P$verbose) displayFailure(theta)
          return(FailureValue(theta))
        }
      }
    }
    if (symmetric && !is.null(P$PoS)) {
      fcn <- function(theta) {
        P$pts <- c(2*P$PoS - theta, P$PoS, theta)
        Pout <- WeightedKDE(P)
        if (Pout$flag == 0) {
          if (P$verbose) displaySuccess(theta, Pout$estSoln$value)
          return(Pout$estSoln$value)
        } else {
          if (P$verbose) displayFailure(theta)
          return(FailureValue(theta))
        }
      }
    }
    if (symmetric && is.null(P$PoS)) {
      fcn <- function(theta) {
        P$PoS <- theta[1]
        P$pts <- c(2*P$PoS - theta[2], P$PoS, theta[2])
        Pout <- WeightedKDE(P)
        if (Pout$flag == 0) {
          if (P$verbose) displaySuccess(theta, Pout$estSoln$value)
          return(Pout$estSoln$value)
        } else {
          if (P$verbose) displayFailure(theta)
          return(FailureValue(theta))
        }
      }
    }
  }

  # Make the function (bimodal case)
  if ("bimodal" %in% P$constraints) {
    if (!symmetric) {
      fcn <- function(theta) {
        P$modeLoc <- theta
        Pout <- WeightedKDE(P)
        if (Pout$flag == 0) {
          if (P$verbose) displaySuccess(theta, Pout$estSoln$value)
          return(Pout$estSoln$value)
        } else {
          if (P$verbose) displayFailure(theta)
          return(FailureValue(theta))
        }
      }
    }
    if (symmetric && !is.null(P$PoS)) {
      fcn <- function(theta) {
        P$modeLoc <- c(2*P$PoS - theta, P$PoS, theta)
        Pout <- WeightedKDE(P)
        if (Pout$flag == 0) {
          if (P$verbose) displaySuccess(theta, Pout$estSoln$value)
          return(Pout$estSoln$value)
        } else {
          if (P$verbose) displayFailure(theta)
          return(FailureValue(theta))
        }
      }
    }
    if (symmetric && is.null(P$PoS)) {
      fcn <- function(theta) {
        P$PoS <- theta[1]
        P$modeLoc <- c(2*P$PoS - theta[2], P$PoS, theta[2])
        Pout <- WeightedKDE(P)
        if (Pout$flag == 0) {
          if (P$verbose) displaySuccess(theta, Pout$estSoln$value)
          return(Pout$estSoln$value)
        } else {
          if (P$verbose) displayFailure(theta)
          return(FailureValue(theta))
        }
      }
    }
  }

  # Make the function (only symmetric case)
  if (symmetric && length(intersect(P$constraints,
      c("unimodal", "twoInflections", "twoInflections+", "bimodal")))==0) {
    fcn <- function(theta) {
      P$PoS <- theta
      Pout <- WeightedKDE(P)
      if (Pout$flag == 0) {
        if (P$verbose) displaySuccess(theta, Pout$estSoln$value)
        return(Pout$estSoln$value)
      } else {
        if (P$verbose) displayFailure(theta)
        return(FailureValue(theta))
      }
    }
  }

  return(fcn)
}

