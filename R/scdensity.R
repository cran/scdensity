# FUNCTIONS IN THIS FILE:
#    scdensity
#    print.scdensity
#    summary.scdensity
#    print.summary.scdensity
#    plot.scdensity


#*****************************************************************************************
#*** scdensity ***************************************************************************
#*****************************************************************************************

#' Shape-constrained kernel density estimation.
#'
#' @description
#' \code{scdensity} computes kernel density estimates that satisfy specified shape
#' restrictions. It is used in the same way as [stats::density()], and takes
#' most of that function's arguments. Its default behavior is to compute a unimodal estimate.
#' Use argument \code{constraint} to choose different shape constraints, \code{method} to
#' choose a different estimation method, and \code{opts} to specify method- and
#' constraint-specific options. The result is a list of S3 class \code{scdensity}, which
#' may be inspected via print, summary, and plot methods.
#'
#' @details
#' All density estimates in this package use the Gaussian kernel.  It is the only common
#' kernel function with three continuous derivatives everywhere.  The \code{adjustedKDE} and
#' \code{weightedKDE} methods require continuous derivatives to ensure numerical stability.
#'
#' The default estimation method, \code{adjustedKDE}, can handle all of the available constraints.  The
#' \code{weightedKDE} method can handle every constraint except \code{symmetric}, while the
#' \code{greedySharpenedKDE} method can handle only \code{unimodal}, \code{monotoneRightTail},
#' \code{monotoneLeftTail}, \code{boundedLeft}, and \code{boundedRight}. The \code{opts} list can
#' also be used to supply method-specific control parameters.  See the "Method details" section
#' for more.
#'
#' Each constraint has a corresponding control parameter that can be supplied as an element of
#' \code{opts}.  The control parameters are described in the following table.  See the "Constraint
#' details" section for definitions of each constraint.
#'
#' \if{html}{\figure{ConstraintsTable.svg}{options: width=650 alt="constraints Table"}}
#' \if{latex}{\figure{ConstraintsTable.pdf}{options: width=5in}}
#'
#' More than one shape constraint can be specified simultaneously.  Certain combinations of constraints
#' (e.g., \code{unimodal} and \code{monotoneRightTail}) are redundant, and will cause a warning. Other
#' combinations (e.g., \code{unimodal} and \code{bimodal}) are incompatible and will cause an error.
#' The figure below summarizes the valid constraint combinations.
#'
#' \if{html}{\figure{ConstraintCombos.svg}{options: width=650 alt="valid constraint combinations"}}
#' \if{latex}{\figure{ConstraintCombos.pdf}{options: width=5in}}
#'
#'
#' # Constraint details
#'
#' All of the constraints other than \code{symmetric} are restrictions on the sign of the estimate, or
#' its derviatives, over certain intervals.  The boundaries of the intervals may be called
#' \emph{important points}. If \code{method="greedySharpenedKDE"}, the important points are determined
#' implicitly during estimation.  For the other methods, the locations of the important points may be
#' supplied in \code{opts}; in most cases they are optional.  If they are not provided, estimation
#' will be run iteratively inside a search routine (\code{\link{SequentialLineMin}}) to find good values,
#' and these values will be returned in the \code{extra} list.
#'
#' Here is a list of the constraints with their definitions and any relevant comments about their
#' usage.
#'
#' * \code{unimodal}: The estimate is nondecreasing to the left of \code{opts$modeLocation}, and
#'   nonincreasing to the right.  If \code{modeLocation} is not supplied, it is found by search.
#' * \code{monotoneRightTail}: The estimate is nonincreasing to the right of the \code{opts$rightTail}
#'   percentile of the unconstrained estimate. \code{rightTail} is a numeric value between 0 and 100.
#'   If it is not supplied, it is set to its default value, 90.
#' * \code{monotoneLeftTail}: The estimate is nondecreasing to the left of the \code{opts$leftTail}
#'   percentile of the unconstrained estimate. \code{leftTail} is a numeric value between 0 and 100. If
#'   it is not supplied, it is set to its default value, 10.
#' * \code{twoInflections}: The estimate has two inflection points, found at
#'   \code{opts$inflectionPoints[1]} and \code{opts$inflectionPoints[2]}. This constraint implies unimodality,
#'   but provides greater smoothness than \code{unimodal}. If \code{inflectionPoints} is not supplied, it is
#'   found by search.
#' * \code{twoInflections+}: The \emph{derivative} of the estimate has three inflection
#'   points, located at \code{opts$inflectionPoints[1]}, \code{opts$inflectionPoints[2]}, and
#'   \code{opts$inflectionPoints[3]}.  This constraint implies \code{twoInflections} but is even smoother.
#'   Most parametric densities with two tails satisfy this constraint.  If \code{inflectionPoints} is not
#'   supplied, it is found by search.
#' * \code{boundedLeft}: The estimate is zero to the left of \code{opts$lowerBound}. The value of
#'   \code{lowerBound} must be specified in \code{opts}. This constraint is implemented only up to a
#'   numerical tolerance. Consequently it is still possible to use it with the Gaussian kernel.
#' * \code{boundedRight}: The estimate is zero to the right of \code{opts$upperBound}. The value of
#'   \code{upperBound} must be specified in \code{opts}. This constraint is
#'   implemented only up to a numerical tolerance.  Consequently it is still possible to use
#'   it with the Gaussian kernel.
#' * \code{symmetric}: The estimate is symmetric around \code{opts$pointOfSymmetry}. If
#'   \code{pointOfSymmetry} is not provided, it is found by search.
#' * \code{bimodal}: The estimate has modes at \code{opts$modeLocation[1]} and \code{opts$modeLocation[3]},
#'   with an antimode (local minimum) at \code{opts$modeLocation[2]}. If \code{modeLocation} is not
#'   specified, it is found by search.
#' 
#'
#' # Method details
#'
#' The \code{adjustedKDE} and \code{weightedKDE} methods are implemented using a common framework
#' where the standard KDE is first approximated by a binning step, after which the constrained estimate
#' is obtained. The \code{greedySharpenedKDE} method uses a different approach.
#'
#'
#' ## adjustedKDE and weightedKDE
#' 
#' The \code{adjustedKDE} method is based on the method of Wolters and Braun (2017).  The method
#' uses the usual unconstrained kernel density estimate as a pilot estimate, and adjusts the shape of
#' this estimate by adding a function to it.  The function is selected to minimally change the
#' shape of the pilot estimate while ensuring the constraints are satisfied. Any of the constraints
#' can be used with this method.
#'
#' The \code{weightedKDE} method is based on the method of Hall and Huang (2002).
#' The method uses a weighted kernel density estimator, with the weights minimally
#' perturbed such that the constraint is satisfied. Any of the constraints except \code{symmetric}
#' may be used with this method.
#'
#' For either of these methods, the following optional arguments can be provided as elements of \code{opts}:
#' * \code{ncheck}:  The number of abscissa points used for constraint checking.  By default,
#'   this is set to \code{max(100, ceiling((diff(range(x)) + 6*h) / h))}, where \code{h} is
#'   the bandwidth. With this default it should be rare to encounter constraint violations large enough
#'   to be visible in a plot.  In the event that constraint violations are observed, re-run the estimation
#'   with a larger value of \code{ncheck}.
#' * \code{verbose}: If \code{TRUE}, progress information will be displayed in the console.
#'   The main use of this is to track the progress of the search for important points. Default is \code{FALSE}.
#'
#' When either of these methods are used, the output list \code{extra} contains elements giving the locations of the
#' important points used in the final estimate (e.g., \code{modeLocation} if the estimate is unimodal or
#' bimodal). Additionally, it containts the following elements:
#' * \code{conCheckGrid}: A vector giving the abscissa values at which the constraints were enforced.
#' * \code{binnedCenters}: A vector giving the locations of the kernel centers determined in the
#'   binning step.
#' * \code{binnedWeights}: The weights corresponding to the binned centers.
#' * \code{finalCenters}: The kernel centers used for the final estimate.
#' * \code{finalWeights}: The weights used for the final estimate.
#'
#' ## greedySharpenedKDE
#' 
#' The \code{greedySharpenedKDE} method is described in Wolters (2012a, 2012b). It uses a data sharpening
#' (shifting the data points) approach.  Starting from an initial solution that satisfies the constraints,
#' a greedy algorithm (implemented in the function \code{\link{improve}}) is used to move the points as close as
#' possible to the observed data while maintaining feasibility.
#'
#' The following optional arguments can be provided as elements of \code{opts}:
#' * \code{startValue} --- A vector of the same length as \code{x}, giving the feasible
#'   initial solution from which the algorithm is started.  If not specified, a vector with
#'   all data points at the location of the unconstrained estimate's highest mode will be used.
#'   Note, it is not guaranteed that the default will satisfy every constraint for every data
#'   set.
#' * \code{verbose}: If \code{TRUE}, information about iteration progress will be printed
#'   to the console. Default is \code{FALSE}.
#' * \code{maxpasses}: Each "pass" through the data points moves each point one-by-one in a greedy fasion.
#'   This option limits the maximum number of passes. Default is 500.
#' * \code{tol}: A numerical tolerance for constraint checking.  See \code{\link{improve}}.
#' * \code{ILS}: An integer greater than zero.  If supplied, the greedy algorithm is run inside an
#'   iterated local search metaheuristic, as described in Wolters (2012b, sec. 3.4). This can improve solution
#'   quality, but requires the greedy search to be run \code{2*ILS} extra times.
#'
#' When this method is used, the output list \code{extra} contains the following elements:
#' * \code{xstar}: The final vector of "sharpened" data points used to generate the
#'   estimate.
#'   
#'
#' # References
#'
#' Hall and Huang (2002), Unimodal Density Estimation Using Kernel Methods, \emph{Statistica Sinica},
#' 12, 965-990.
#'
#' Wolters and Braun (2017), Enforcing Shape Constraints on a Probability Density Estimate Using an Additive
#' Adjustment curve, \emph{Communications in Statistics - Simulation and Computation},
#' 47(3), 672-691.
#'
#' Wolters (2012a), A Greedy Algorithm for Unimodal Kernel Density Estimation by Data Sharpening,
#' \emph{Journal of Statistical Software}, 46(6), 1–26.
#'
#' Wolters (2012b), Methods for Shape-Constrained Kernel Density Estimation. Ph.D. Thesis, University
#' of Western Ontario.
#'
#' @param x A vector of data from which the estimate is to be computed.
#' @param bw The bandwidth.  It is specified as either a numerical value or as one of the
#' character strings \code{"nrd0"},  \code{"nrd"}, \code{"ucv"}, \code{"bcv"}, or
#' \code{"SJ"}, exactly as in [stats::density()].
#' @param constraint A vector of strings giving the operative shape constraints. Elements
#' must partially match different alternatives among \code{"unimodal"},
#' \code{"monotoneRightTail"},\code{"monotoneLeftTail"}, \code{"twoInflections"},
#' \code{"twoInflections+"}, \code{"boundedLeft"}, \code{"boundedRight"},
#' \code{"symmetric"}, and \code{"bimodal"}.
#' @param method A string giving the method of enforcing shape constraints.  It must
#' paritally match one of \code{"adjustedKDE"}, \code{"weightedKDE"}, or
#' \code{"greedySharpenedKDE"}.
#' @param opts A list giving options specific to the chosen constraints and/or method. E.g.
#' use \code{opts = list(modeLocation = 0)} to force the mode to be at zero when the
#' constraint is \code{unimodal}. See below for lists of
#' available options.
#' @param adjust A scaling factor for the bandwidth, just as in [stats::density()].
#' @param n The number of points returned in the density estimate.  Same as in
#' [stats::density()].
#' @param na.rm Logical indicating whether or not to remove missing values from \code{x}.
#' Same as in [stats::density()].
#'
#' @return A list with the following elements:
#' * \code{constraint} The constraint(s) used for estimation.  Might differ from
#'   the constraints supplied to the function if they included redundant constraints.
#' * \code{method} The estimation method used.
#' * \code{f0} A function.  Use \code{f0(v)} to evaluate the unconstrained KDE at the points in
#'   \code{v}.
#' * \code{fhat} A function. Use \code{fhat(v)} to evaluate the constrained KDE at the points in
#'   \code{v}.
#' * \code{data} The data used to generate the estimate.
#' * \code{bw} The bandwidth used.
#' * \code{extra} A list holding additional outputs that are specific to the chosen method.
#'   See the "method details" section. 
#' * \code{x} A vector of abscissa values for plotting the estimate.  Same as in
#'   [stats::density()].
#' * \code{y} A vector of ordinate values for plotting the estimate.  Same as in
#'   [stats::density()].
#' * \code{n} The sample size, not including missing values.  Note, this \code{n} has
#'   no relation to the \code{n} provided in the arguments.
#' * \code{data.name} Deparsed name of the \code{x} argument, used in plotting.
#' * \code{call} The call to the function.
#' * \code{has.na} Always \code{FALSE}.  Included for consistency with
#'   [stats::density()].
#'
#' @importFrom stats approx approxfun bw.SJ bw.bcv bw.nrd bw.nrd0 bw.ucv density rnorm
#'
#' @export
#' @aliases density
#' @seealso \code{\link{plot.scdensity}} plot method, \code{\link{print.scdensity}} print
#' method, and \code{\link{summary.scdensity}} summary method.
#' @examples
#' # Default method gives a unimodal estimate using adjustment curve method.
#' x <- rlnorm(30)
#' scKDE <- scdensity(x)
#' scKDE
#' summary(scKDE)
#' plot(scKDE, detail=2)
#' plot(scKDE, detail=4)
#'
#' # Constrain the first and fourth quartiles to be monotone, using greedy sharpening method.
#' x <- rt(50, df=3)
#' scKDE <- scdensity(x, bw="SJ", adjust=0.5, constraint=c("monotoneL", "monotoneR"),
#'                    opts=list(verbose=TRUE, leftTail=25, rightTail=75), method="greedy")
#' plot(scKDE)
#'
#' # Compare unimodal, twoInflections, and twoInflections+ constraints
#' x <- rnorm(100)
#' h <- 0.5 * bw.SJ(x)
#' fhat1 <- scdensity(x, bw=h, constraint="unimodal")
#' fhat2 <- scdensity(x, bw=h, constraint="twoInflections")
#' fhat3 <- scdensity(x, bw=h, constraint="twoInflections+")
#' plot(density(x, bw=h))
#' lines(fhat1$x, fhat1$y, col="red")
#' lines(fhat2$x, fhat2$y, col="blue")
#' lines(fhat3$x, fhat3$y, col="green", lwd=2)
#'
scdensity <- function(x,
                      bw = "nrd0",
                      constraint = c("unimodal", "monotoneRightTail", "monotoneLeftTail",
                                     "twoInflections", "twoInflections+", "boundedLeft",
                                     "boundedRight", "symmetric", "bimodal"),
                      method = c("adjustedKDE", "weightedKDE", "greedySharpenedKDE"),
                      opts=NULL,
                      adjust = 1,
                      n = 512,
                      na.rm = FALSE) {

  #=== Input checking =====================================

  #--- Check inputs in common with density() ----
  # This code was copied and adapted from density()
  if (!is.numeric(x))
    stop("argument 'x' must be numeric")
  name <- deparse(substitute(x))
  x <- as.vector(x)
  x.na <- is.na(x)
  if (any(x.na)) {
    if (na.rm)
      x <- x[!x.na]
    else
      stop("'x' contains missing values")
  }

  x.finite <- is.finite(x)
  if (any(!x.finite)) {
    x <- x[x.finite]
  }
  nx <- length(x)

  if (is.character(bw)) {
    if (nx < 2)
      stop("need at least 2 points to select a bandwidth automatically")
    bw <- switch(tolower(bw), nrd0 = bw.nrd0(x), nrd = bw.nrd(x), ucv = bw.ucv(x),
                 bcv = bw.bcv(x), sj = , `sj-ste` = bw.SJ(x, method = "ste"),
                 `sj-dpi` = bw.SJ(x, method = "dpi"), stop("unknown bandwidth rule"))
  }
  if (!is.finite(bw))
    stop("non-finite 'bw'")
  bw <- adjust * bw
  if (bw <= 0)
    stop("'bw' is not positive.")

  #--- Set up rescaling ---
  # Scale the problem to a fixed range to ensure predictable numerics.  Also set from, to.
  x.orig <- x
  bw.orig <- bw
  from.orig <- min(x) - 4*bw
  to.orig <- max(x) + 4*bw
  pilotEstimate.orig <- density(x=x.orig, bw=bw.orig, kernel="gaussian",
                                from=from.orig, to=to.orig)
  oldRange <- range(x)
  newRange <- c(-1, 1)
  ScaleDown <- function(y) newRange[1] + diff(newRange)*(y-oldRange[1])/diff(oldRange)
  ScaleUp <- function(y) oldRange[1] + diff(oldRange)*(y-newRange[1])/diff(newRange)
  x <- ScaleDown(x)
  bw <- bw * diff(newRange)/diff(oldRange)
  from <- min(x) - 4 * bw
  to <- max(x) + 4 * bw

  #--- Handle the constraint specification -------
  if ("constraint" %in% names(as.list(match.call()))) {
    userConstraints <- constraint
    constraint <- match.arg(constraint, several.ok = TRUE)
    if (length(constraint) != length(userConstraints))
      warning("Some elements of 'constraint' could not be matched and have been discarded.")
  } else {
    constraint <- match.arg(constraint)
  }

  #--- Check constraints for redundancies or incompatibilities ---
  # Ensure the following relationships are upheld in the constraints:
  #  1) twoInflections+ ==> twoInflections ==> unimodal ==> monotone tails.
  #  2) bimodal ==> monotone tails
  #  3) bimodal is incompatible with unimodal, twoInflections, and twoInflections+
  constraint.user <- constraint
  if ("twoInflections+" %in% constraint) {
    throwOut <- c("twoInflections", "unimodal", "monotoneRightTail", "monotoneLeftTail")
    constraint <- setdiff(constraint, throwOut)
  }
  if ("twoInflections" %in% constraint) {
    throwOut <- c("unimodal", "monotoneRightTail", "monotoneLeftTail")
    constraint <- setdiff(constraint, throwOut)
  }
  if ("unimodal" %in% constraint) {
    throwOut <- c("monotoneRightTail", "monotoneLeftTail")
    constraint <- setdiff(constraint, throwOut)
  }
  if ("bimodal" %in% constraint) {
    throwOut <- c("monotoneRightTail", "monotoneLeftTail")
    constraint <- setdiff(constraint, throwOut)
    notOK <- c("unimodal", "twoInflections", "twoInflections+")
    if (any(notOK %in% constraint))
      stop("One or more elements of 'constraint' are incompatible with 'bimodal'")
  }
  if (length(constraint) < length(constraint.user))
    warning("Redundant members of 'constraint' have been removed")

  #--- Is method compatible with constraint? ----
  method <- match.arg(method)
  weightedOK <- c("unimodal", "monotoneRightTail", "monotoneLeftTail", "twoInflections",
                  "twoInflections+", "boundedLeft", "boundedRight", "bimodal")
  greedyOK <- c("unimodal", "monotoneRightTail", "monotoneLeftTail", "boundedLeft",
                "boundedRight")
  adjustedOK <- c("unimodal", "monotoneRightTail", "monotoneLeftTail", "twoInflections",
                  "twoInflections+", "boundedLeft", "boundedRight", "symmetric", "bimodal")
  if (method == "weightedKDE" && !all(constraint %in% weightedOK))
    stop("The weightedKDE method is not compatible with the chosen constraint(s)")
  if (method == "greedySharpenedKDE" && !all(constraint %in% greedyOK))
    stop("The greedySharpenedKDE method is not compatible with the chosen constraint(s)")
  if (method == "adjustedKDE" && !all(constraint %in% adjustedOK))
    stop("The adjustedKDE method is not compatible with the chosen constraint(s)")


  #--- unpack constraint options from opts -------
  # opts contains options both for the constraints and for the method. Here pull out the
  # constraint-related options into local variables, and delete them from the opts list.
  # Then we can pass opts on to the method function calls.
  if ("unimodal" %in% constraint) {
    if ("modeLocation" %in% names(opts)) {
      modeLocation <- ScaleDown(opts$modeLocation)
      opts$modeLocation <- NULL
      stopifnot(
        length(modeLocation) == 1,
        modeLocation <= max(x),
        modeLocation >= min(x)
      )
    } else {
      modeLocation <- NULL
    }
  }
  if ("monotoneRightTail" %in% constraint) {
    if ("rightTail" %in% names(opts)) {
      rightTail <- opts$rightTail
      opts$rightTail <- NULL
      stopifnot(
        length(rightTail) == 1,
        rightTail > 0,
        rightTail < 100
      )
    } else {
      rightTail <- 90
    }
  }
  if ("monotoneLeftTail" %in% constraint) {
    if ("leftTail" %in% names(opts)) {
      leftTail <- opts$leftTail
      opts$leftTail <- NULL
      stopifnot(
        length(leftTail) == 1,
        leftTail > 0,
        leftTail < 100
      )
    } else {
      leftTail <- 10
    }
    if (exists("rightTail")) {
      stopifnot(leftTail < rightTail)
    }
  }
  if ("twoInflections" %in% constraint) {
    if ("inflectionPoints" %in% names(opts)) {
      inflectionPoints <- ScaleDown(opts$inflectionPoints)
      opts$inflectionPoints <- NULL
      stopifnot(
        length(inflectionPoints) == 2,
        min(inflectionPoints) > min(x),
        max(inflectionPoints) < max(x),
        inflectionPoints[1] != inflectionPoints[2]
      )
      inflectionPoints <- sort(inflectionPoints)
    } else {
      inflectionPoints <- NULL
    }
  }
  if ("twoInflections+" %in% constraint) {
    if ("inflectionPoints" %in% names(opts)) {
      inflectionPoints <- ScaleDown(opts$inflectionPoints)
      opts$inflectionPoints <- NULL
      stopifnot(
        length(inflectionPoints) == 3,
        min(inflectionPoints) > min(x),
        max(inflectionPoints) < max(x),
        length(unique(inflectionPoints)) == 3
      )
      inflectionPoints <- sort(inflectionPoints)
    } else {
      inflectionPoints <- NULL
    }
  }
  if ("boundedLeft" %in% constraint) {
    if ("lowerBound" %in% names(opts)) {
      lowerBound <- ScaleDown(opts$lowerBound)
      opts$lowerBound <- NULL
      stopifnot(
        length(lowerBound) == 1,
        lowerBound < max(x)
      )
    } else {
      stop("Constraint 'boundedLeft' selected but no 'lowerBound' in opts list")
    }
  }
  if ("boundedRight" %in% constraint) {
    if ("upperBound" %in% names(opts)) {
      upperBound <- ScaleDown(opts$upperBound)
      opts$upperBound <- NULL
      stopifnot(
        length(upperBound) == 1,
        upperBound > min(x)
      )
    } else {
      stop("Constraint 'boundedRight' selected but no 'upperBound' in opts list")
    }
    if (exists("lowerBound")) {
      stopifnot(lowerBound < upperBound)
    }
  }
  if ("symmetric" %in% constraint) {
    symmetric <- TRUE
    if ("pointOfSymmetry" %in% names(opts)) {
      pointOfSymmetry <- ScaleDown(opts$pointOfSymmetry)
      opts$pointOfSymmetry <- NULL
      stopifnot(
        length(pointOfSymmetry) == 1,
        pointOfSymmetry > min(x),
        pointOfSymmetry < max(x)
      )
    } else {
      pointOfSymmetry <- NULL
    }
  } else {
    symmetric <- FALSE
  }
  if ("bimodal" %in% constraint) {
    if ("modeLocation" %in% names(opts)) {
      modeLocation <- ScaleDown(opts$modeLocation)
      opts$modeLocation <- NULL
      stopifnot(
        length(modeLocation) == 3,
        min(modeLocation) > min(x),
        max(modeLocation) < max(x),
        length(unique(modeLocation)) == 3
      )
      modeLocation <- sort(modeLocation)
    } else {
      modeLocation <- NULL
    }
  }

  #=== Create the output list =============================
  pilotEstimate <- density(x=x, bw=bw, kernel="gaussian", from=from, to=to)
  out <- list(constraint = constraint,
              method = method,
              f0 = approxfun(pilotEstimate.orig$x, pilotEstimate.orig$y, yleft=0, yright=0),
              fhat = NULL,             #-To be filled.
              data = x.orig,
              bw = bw.orig,
              extra = list(),          #-To be filled.
              x = NULL,                #-To be filled.
              y = NULL,                #-To be filled.
              n = length(x),
              data.name = name,
              call = match.call(),
              has.na = FALSE)
  class(out) <- "scdensity"

  #=== Weighted KDE or Adjusted KDE methods ===============
  if (method %in% c("weightedKDE", "adjustedKDE")) {

    #--- Initialize the problem object P and populate its members ---
    # Include opts in the args so we can put any method-related options into P.
    P <- InitializeP(x=x, h=bw, constraints=constraint, method=method, opts=opts)
    P$pilotModeLoc <- pilotEstimate$x[which.max(pilotEstimate$y)]

    #--- Set up constraints that never require a search ---
    P$constraints <- constraint
    if ("boundedLeft" %in% constraint) {
      P$lowerBound <- lowerBound
    }
    if ("boundedRight" %in% constraint) {
      P$upperBound <- upperBound
    }
    if ("monotoneRightTail" %in% constraint) {
      P$rightPoint <- getQuantile(pilotEstimate$x, pilotEstimate$y, rightTail/100)
      out$extra$rightPoint <- ScaleUp(P$rightPoint)
    }
    if ("monotoneLeftTail" %in% constraint) {
      P$leftPoint <- getQuantile(pilotEstimate$x, pilotEstimate$y, leftTail/100)
      out$extra$leftPoint <- ScaleUp(P$leftPoint)
    }

    #--- Pre-compute concheck grid and binning for non-symmetry cases -----
    if (!symmetric) {
      P <- BuildConCheckGrid(P)
      P <- BinningStep(P)
    }

    #--- Set the tolerance and initial bounds for searches -----
    searchTol <- diff(range(x))/1000
    searchBounds <- numeric(2)
    if (exists("lowerBound")) {
      searchBounds[1] <- lowerBound
    } else {
      searchBounds[1] <- min(x)
    }
    if (exists("upperBound")) {
      searchBounds[2] <- upperBound
    } else {
      searchBounds[2] <- max(x)
    }

    #--- If UNIMODAL, set or find mode location -----
    if ("unimodal" %in% constraint) {

      # No symmetry case.
      if (!symmetric) {
        if (!is.null(modeLocation)) {
          P$modeLoc <- modeLocation
        } else {
          ObjFun <- makeOF(P)
          thetaBest <- SequentialLineMin(ObjFun, searchBounds, v0=P$pilotModeLoc, tol=searchTol)
          P$modeLoc <- thetaBest$minimizer
        }
      }

      # Symmetry case.
      # Note: 0, 1, or 2 of {modeLocation, pointOfSymmetry} could be given.
      if (symmetric) {
        case <- is.null(modeLocation) + is.null(pointOfSymmetry)
        if (case == 0) {
          if (modeLocation != pointOfSymmetry) {
            stop("Specified both modeLocation and pointOfSymmetry, but they are not equal.")
          } else {
            P$modeLoc <- modeLocation
            P$PoS <- pointOfSymmetry
          }
        }
        if (case == 1) {
          theNonNullValue <- max(modeLocation, pointOfSymmetry)
          P$modeLoc <- theNonNullValue
          P$PoS <- theNonNullValue
        }
        if (case == 2) {
          ObjFun <- makeOF(P)
          thetaBest <- SequentialLineMin(ObjFun, searchBounds, v0=P$pilotModeLoc, tol=searchTol)
          P$modeLoc <- thetaBest$minimizer
          P$PoS <- thetaBest$minimizer
        }
        out$extra$pointOfSymmetry <- ScaleUp(P$PoS)
      }

      out$extra$modeLocation <- ScaleUp(P$modeLoc)
    }

    #--- If TWO INFLECTIONS, set or find important points -----
    if ("twoInflections" %in% constraint) {

      # No symmetry case.
      if (!symmetric) {
        if (!is.null(inflectionPoints)) {
          P$pts  <- inflectionPoints
        } else {
          ObjFun <- makeOF(P)
          v0 <- c(getQuantile(pilotEstimate$x, pilotEstimate$y, 0.3),
                  getQuantile(pilotEstimate$x, pilotEstimate$y, 0.7))
          thetaBest <- SequentialLineMin(ObjFun, searchBounds, v0=v0, tol=searchTol)
          P$pts <- thetaBest$minimizer
        }
      }

      # Symmetry case.
      # Note: 0, 1, or 2 of {inflectionPoints, pointOfSymmetry} could be given.
      if (symmetric) {
        ipOK <- !is.null(inflectionPoints)
        posOK <- !is.null(pointOfSymmetry)
        if (ipOK && posOK) {
          if (!isTRUE(all.equal(sum(inflectionPoints), 2*pointOfSymmetry))) {
            stop("Specified inflectionPoints are not symmetric around pointOfSymmetry.")
          } else {
            P$pts <- inflectionPoints
            P$PoS <- pointOfSymmetry
          }
        }
        if (ipOK && !posOK) {
          P$pts <- inflectionPoints
          P$PoS <- sum(inflectionPoints)/2
        }
        if (!ipOK && posOK) {
          P$PoS <- pointOfSymmetry
          ObjFun <- makeOF(P)
          searchBounds[1] <- pointOfSymmetry
          thetaBest <- SequentialLineMin(ObjFun, searchBounds, v0=P$pilotModeLoc, tol=searchTol)
          P$pts <- c(2*P$PoS - thetaBest$minimizer, thetaBest$minimizer)
        }
        if (!ipOK && !posOK) {
          ObjFun <- makeOF(P)
          v0 <- c(getQuantile(pilotEstimate$x, pilotEstimate$y, 0.33),
                  getQuantile(pilotEstimate$x, pilotEstimate$y, 0.66))
          thetaBest <- SequentialLineMin(ObjFun, searchBounds, v0=v0, tol=searchTol)
          P$PoS <- thetaBest$minimizer[1]
          P$pts <- c(2*P$PoS - thetaBest$minimizer[2], thetaBest$minimizer[2])
        }
        out$extra$pointOfSymmetry <- ScaleUp(P$PoS)
      }

      out$extra$inflectionPoints <- ScaleUp(P$pts)
    }


    #--- If TWO INFLECTIONS+, set or find important points -----
    if ("twoInflections+" %in% constraint) {

      # No symmetry case.
      if (!symmetric) {
        if (!is.null(inflectionPoints)) {
          P$pts  <- inflectionPoints
        } else {
          ObjFun <- makeOF(P)
          v0 <- c(getQuantile(pilotEstimate$x, pilotEstimate$y, 0.1),
                  getQuantile(pilotEstimate$x, pilotEstimate$y, 0.5),
                  getQuantile(pilotEstimate$x, pilotEstimate$y, 0.9))
          thetaBest <- SequentialLineMin(ObjFun, searchBounds, v0=v0, tol=searchTol)
          P$pts <- thetaBest$minimizer
        }
      }

      # Symmetry case.
      # Note: 0, 1, or 2 of {inflectionPoints, pointOfSymmetry} could be given.
      if (symmetric) {
        ipOK <- !is.null(inflectionPoints)
        posOK <- !is.null(pointOfSymmetry)
        if (ipOK && posOK) {
          check1 <- isTRUE(all.equal(inflectionPoints[2], pointOfSymmetry))
          check2 <- isTRUE(all.equal(sum(inflectionPoints[c(1,3)]), 2*pointOfSymmetry))
          if (!(check1 & check2)) {
            stop("Specified inflectionPoints are not symmetric around pointOfSymmetry.")
          } else {
            P$pts <- inflectionPoints
            P$PoS <- pointOfSymmetry
          }
        }
        if (ipOK && !posOK) {
          if (!isTRUE(all.equal(sum(inflectionPoints[c(1,3)]), 2*inflectionPoints[2]))) {
            stop("With symmetry, inflection points 1 and 3 must be equidistant from 2.")
          }
          P$pts <- inflectionPoints
          P$PoS <- inflectionPoints[2]
        }
        if (!ipOK && posOK) {
          P$PoS <- pointOfSymmetry
          ObjFun <- makeOF(P)
          searchBounds[1] <- pointOfSymmetry
          thetaBest <- SequentialLineMin(ObjFun, searchBounds, v0=P$pilotModeLoc, tol=searchTol)
          P$pts <- c(2*P$PoS - thetaBest$minimizer, P$PoS, thetaBest$minimizer)
        }
        if (!ipOK && !posOK) {
          ObjFun <- makeOF(P)
          v0 <- c(getQuantile(pilotEstimate$x, pilotEstimate$y, 0.5),
                  getQuantile(pilotEstimate$x, pilotEstimate$y, 0.9))
          thetaBest <- SequentialLineMin(ObjFun, searchBounds, v0=v0, tol=searchTol)
          P$PoS <- thetaBest$minimizer[1]
          P$pts <- c(2*P$PoS - thetaBest$minimizer[2], P$PoS, thetaBest$minimizer[2])
        }
        out$extra$pointOfSymmetry <- ScaleUp(P$PoS)
      }

      out$extra$inflectionPoints <- ScaleUp(P$pts)
    }


    #--- If BIMODAL, set or find important points -----
    if ("bimodal" %in% constraint) {

      # No symmetry case.
      if (!symmetric) {
        if (!is.null(modeLocation)) {
          P$modeLoc  <- modeLocation
        } else {
          ObjFun <- makeOF(P)
          v0 <- c(getQuantile(pilotEstimate$x, pilotEstimate$y, 0.3),
                  getQuantile(pilotEstimate$x, pilotEstimate$y, 0.5),
                  getQuantile(pilotEstimate$x, pilotEstimate$y, 0.7))
          thetaBest <- SequentialLineMin(ObjFun, searchBounds, v0=v0, tol=searchTol)
          P$modeLoc <- thetaBest$minimizer
        }
      }

      # Symmetry case.
      # Note: 0, 1, or 2 of {modeLocation, pointOfSymmetry} could be given.
      if (symmetric) {
        modesOK <- !is.null(modeLocation)
        posOK <- !is.null(pointOfSymmetry)
        if (modesOK && posOK) {
          check1 <- isTRUE(all.equal(modeLocation[2], pointOfSymmetry))
          check2 <- isTRUE(all.equal(sum(modeLocation[c(1,3)]), 2*pointOfSymmetry))
          if (!(check1 & check2)) {
            stop("Specified modeLocation values are not symmetric around pointOfSymmetry")
          } else {
            P$modeLoc <- modeLocation
            P$PoS <- pointOfSymmetry
          }
        }
        if (modesOK && !posOK) {
          if (!isTRUE(all.equal(sum(modeLocation[c(1,3)]), 2*modeLocation[2]))) {
            stop("With symmetry, modeLocation points 1 and 3 must be equidistant from 2.")
          }
          P$modeLoc <- modeLocation
          P$PoS <- modeLocation[2]
        }
        if (!modesOK && posOK) {
          P$PoS <- pointOfSymmetry
          ObjFun <- makeOF(P)
          searchBounds[1] <- pointOfSymmetry
          thetaBest <- SequentialLineMin(ObjFun, searchBounds, v0=P$pilotModeLoc, tol=searchTol)
          P$modeLoc <- c(2*P$PoS - thetaBest$minimizer, P$PoS, thetaBest$minimizer)
        }
        if (!modesOK && !posOK) {
          ObjFun <- makeOF(P)
          v0 <- c(getQuantile(pilotEstimate$x, pilotEstimate$y, 0.5),
                  getQuantile(pilotEstimate$x, pilotEstimate$y, 0.7))
          thetaBest <- SequentialLineMin(ObjFun, searchBounds, v0=v0, tol=searchTol)
          P$PoS <- thetaBest$minimizer[1]
          P$modeLoc <- c(2*P$PoS - thetaBest$minimizer[2], P$PoS, thetaBest$minimizer[2])
        }
        out$extra$pointOfSymmetry <- ScaleUp(P$PoS)
      }

      out$extra$modeLocation <- ScaleUp(P$modeLoc)
    }

    #--- If ONLY SYMMETRIC, set or find point of symmetry -----
    #"only symmetric" means symmetry without unimodal, twoInflections, twoInflections+,
    #or bimodal constraints.
    if (symmetric && length(intersect(constraint,
        c("unimodal", "twoInflections", "twoInflections+", "bimodal")))==0) {
      if (!is.null(pointOfSymmetry)) {
        P$PoS <- pointOfSymmetry
      } else {
        ObjFun <- makeOF(P)
        if (exists("leftTail") && searchBounds[1] < P$leftPoint) {
          searchBounds[1] <- P$leftPoint
        }
        if (exists("rightTail") && P$rightPoint < searchBounds[2]) {
          searchBounds[2] <- P$rightPoint
        }
        PoSBest <- SequentialLineMin(ObjFun, searchBounds, v0=P$pilotModeLoc, tol=searchTol)
        P$PoS <- PoSBest$minimizer
      }

      out$extra$pointOfSymmetry <- ScaleUp(P$PoS)
    }


    #--- Calculate the estimate (given important points) -----
    if (P$verbose) {
      cat("--- Important points identified. Now obtain final estimate. ---\n")
    }
    Pout <- WeightedKDE(P)
    if (P$verbose) {
      cat("Final solution obtained. Objective value ", Pout$estSoln$value, "\n")
    }

    #--- Handle error state (infeasible or numerical problems) -----
    if (Pout$flag == 1 && method=="weightedKDE") {
      stop(paste0("Could not find a feasible solution. Change constraints or important points,\n",
           "or use method = adjustedKDE"))
    }
    if (Pout$flag == 1 && method=="adjustedKDE") {
      stop("Could not find a feasible solution.\nCheck the constraints and/or important points.")
    }
    if (Pout$flag == 2) {
      stop(paste0("solve.QP encountered numerical problems and could not find a solution.\n",
                  "Small changes to the inputs may resolve this."))
    }

    #--- Prepare the output list -----
    constrainedEstimate <- density(x=ScaleUp(Pout$y), bw=Pout$h * diff(oldRange)/diff(newRange),
                                   kernel="gaussian", n=n, weights=Pout$vhat)
    out$x <- constrainedEstimate$x
    out$y <- constrainedEstimate$y
    out$fhat <- approxfun(constrainedEstimate$x, constrainedEstimate$y, yleft=0, yright=0)
    out$extra$conCheckGrid <- ScaleUp(Pout$g)
    out$extra$binnedCenters <- ScaleUp(Pout$s)
    out$extra$binnedWeights <- Pout$w
    out$extra$finalCenters <- ScaleUp(Pout$y)
    out$extra$finalWeights <- Pout$vhat
  }

  #=== greedy sharpening method ========================
  if (method == "greedySharpenedKDE") {

    if ("monotoneRightTail" %in% constraint) {
      q <- getQuantile(pilotEstimate$x, pilotEstimate$y, rightTail/100)
      rightTail.ix <- min(which(pilotEstimate$x >= q))
    }
    if ("monotoneLeftTail" %in% constraint) {
      q <- getQuantile(pilotEstimate$x, pilotEstimate$y, leftTail/100)
      leftTail.ix <- max(which(pilotEstimate$x <= q))
    }
    if ("boundedLeft" %in% constraint) {
      lowerBound.ix <- max(which(pilotEstimate$x < lowerBound))
    }
    if ("boundedRight" %in% constraint) {
      upperBound.ix <- min(which(pilotEstimate$x > upperBound))
    }

    confun <- function(x) {
      fvals <- density(x=x, bw=bw, kernel="gaussian", from=from, to=to)$y
      isuni <- !("unimodal" %in% constraint) || isUnimodal(fvals)
      monoR <- !("monotoneRightTail" %in% constraint) || isMonotoneR(fvals, rightTail.ix)
      monoL <- !("monotoneLeftTail" %in% constraint) || isMonotoneL(fvals, leftTail.ix)
      boundedL <- !("boundedLeft" %in% constraint) || isBoundedL(fvals, lowerBound.ix)
      boundedR <- !("boundedRight" %in% constraint) || isBoundedR(fvals, upperBound.ix)

      return(isuni && monoR && monoL && boundedL && boundedR)
    }

    if ("startValue" %in% names(opts)){
      opts$startValue <- ScaleDown(opts$startValue)
    } else {
      pilotMaxLocation <- pilotEstimate$x[which(pilotEstimate$y==max(pilotEstimate$y))]
      opts$startValue <- rep(pilotMaxLocation, nx)
    }
    if ("tol" %in% names(opts)){
      opts$tol <- opts$tol * diff(newRange)/diff(oldRange)
    }
    if ("ILS" %in% names(opts)) {
      stopifnot(
        length(opts$ILS) == 1,
        opts$ILS > 0
      )
      ILS <- opts$ILS
      opts$ILS <- NULL
    } else {
      ILS <- 0
    }
    argslist <- c(list(x=x, confun=confun), opts)
    xstar <- do.call("improve", argslist)

    # Do ILS procedure if requested.
    if (ILS > 0) {
      if ("verbose" %in% names(opts)) {
        cat("\n--- Begin iterated local searches ---\n")
      }
      for (i in 1:ILS) {
        # perturb the previous solution (it likely won't satisfy the constraint)
        noise <- rnorm(nx) * bw
        perturbed <- xstar + noise
        # run improve() with xstar as target to make it satisfy the constraint
        if ("verbose" %in% names(opts)) {
          cat("--- Iteration ", i, " repair step ---\n")
        }
        argslist$startValue <- xstar
        argslist$x <- perturbed
        repaired <- do.call("improve", argslist)
        # run improve() to get the perturbed, repaired solution closer to x.
        if ("verbose" %in% names(opts)) {
          cat("--- Iteration ", i, " improve step ---\n")
        }
        argslist$startValue <- repaired
        argslist$x <- x
        candidate <- do.call("improve", argslist)
        # keep the candidate solution if it's better than xstar
        if (sum(abs(candidate - x)) < sum(abs(xstar - x))) {
          xstar <- candidate
        }
      }
    }

    constrainedEstimate <- density(x=ScaleUp(xstar), bw=bw.orig, kernel="gaussian",
                                   from=from.orig, to=to.orig, n=n)
    out$extra$xstar <- ScaleUp(xstar)
    out$x <- constrainedEstimate$x
    out$y <- constrainedEstimate$y
    out$fhat <- approxfun(constrainedEstimate$x, constrainedEstimate$y, yleft=0, yright=0)
  }

  #=== Return output ========================

  return(out)

}


#*****************************************************************************************
#*** print.scdensity *********************************************************************
#*****************************************************************************************

#' Print method for class \code{scdensity}.
#'
#' Displays the names of the elements of the scdensity list object
#' and their sizes and types. Includes minimal comments about the most important elements.
#'
#' @param x An object of S3 class \code{scensity}.
#' @param ... Included for consistency with generic functions.
#'
#' @export
print.scdensity <- function(x, ...) {

  x.names <- names(x)
  m <- length(x.names)
  extraLoc <- which(x.names=="extra")
  r <- length(x$extra)
  Names <- rep("", m+r-1)
  Lengths <- rep(0, m+r-1)
  Modes <- rep("", m+r-1)
  Comments <- rep("", m+r-1)

  for (i in 1:(extraLoc-1)) {
    Lengths[i] <- length(x[[i]])
    Modes[i] <- mode(x[[i]])
    Names[i] <- x.names[i]
  }
  for (i in seq_along(x$extra)) {
    Lengths[extraLoc-1+i] <- length(x$extra[[i]])
    Modes[extraLoc-1+i] <- mode(x$extra[[i]])
    Names[extraLoc-1+i] <- paste0("extra$", names(x$extra)[i])
  }
  for (i in (extraLoc+1):m) {
    Lengths[i+r-1] <- length(x[[i]])
    Modes[i+r-1] <- mode(x[[i]])
    Names[i+r-1] <- x.names[i]
  }

  smry <- data.frame(Length=Lengths, Mode=Modes, Comment=Comments, stringsAsFactors=FALSE)
  rownames(smry) <- Names

  smry['constraint', 'Comment'] <- paste(x$constraint, collapse=", ")
  smry['method', 'Comment'] <- x$method
  smry['f0', 'Comment'] <- "unconstrained estimator is f0(x)"
  smry['fhat', 'Comment'] <- "constrained estimator is fhat(x)"
  smry['data', 'Comment'] <- "The data vector"
  smry['bw', 'Comment'] <- paste0("Bandwidth = ", zapsmall(x$bw,3))
  smry['x', 'Comment'] <- "For quick plotting"
  smry['y', 'Comment'] <- "For quick plotting"

  cat("A shape-constrained density estimate.\nA list with the following elements:\n\n")
  print(smry)

}

#*****************************************************************************************
#*** summary.scdensity *******************************************************************
#*****************************************************************************************

#' Summary method for class \code{scdensity}.
#'
#' Collects high-level information about the \code{scdensity} object and some descriptive
#' statistics.
#'
#' @param object An object of S3 class \code{scensity}.
#' @param ... Included for consistency with generic functions.
#'
#' @importFrom stats integrate
#' @export
summary.scdensity <- function(object, ...) {

  out <- list()
  out$ndata <- length(object$data)
  out$constraintString <- paste(object$constraint, collapse=", ")
  out$methodString <- object$method
  LB <- min(object$x) - 4*object$bw
  UB <- max(object$x) + 4*object$bw

  out$TVdist <- integrate(function(y) abs(object$f0(y)-object$fhat(y)), LB, UB,
                          rel.tol=10^-3, stop.on.error=FALSE)$value
  out$mean <- integrate(function(y) y*object$fhat(y), LB, UB,
                        rel.tol=10^-3, stop.on.error=FALSE)$value
  out$sd <- sqrt(integrate(function(y) (y-out$mean)^2*object$fhat(y), LB, UB,
                            rel.tol=10^-3, stop.on.error=FALSE)$value)

  out$percentPoints <- 100*c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
  Q <- getQuantileFunction(object)
  out$percentiles <- Q(out$percentPoints/100)
  names(out$percentiles) <- paste0(out$percentPoints, "%")

  class(out) <- "summary.scdensity"
  out

}

#*****************************************************************************************
#*** print.summary.scdensity *************************************************************
#*****************************************************************************************

#' Prints the information in a \code{summary.scdensity} object to the console.
#'
#' @param x An object of S3 class \code{summary.scdensity}.
#' @param ... Included for consistency with generic functions.
#' @export
print.summary.scdensity <- function(x, ...) {

  cat("A constrained density estimate based on ", x$ndata, " observations.\n", sep="")
  cat("  Constraint: ", x$constraintString, "\n", sep="")
  cat("  Method: ", x$methodString, "\n\n", sep="")

  cat("Total variation distance between constrained\n")
  cat("and unconstrained estimates: ", zapsmall(x$TVdist,3), ".\n\n", sep="")

  cat("The mean of the estimated distribution: ", zapsmall(x$mean,3), "\n", sep="")
  cat("The standard deviation of the estimated distribution: ", zapsmall(x$sd, 3), "\n\n",
      sep="")

  cat("Percentiles of the distribution:\n")
  print(round(x$percentiles, digits = 3))

}


#*****************************************************************************************
#*** plot.scdensity *********************************************************************
#*****************************************************************************************

#' Plot method for class \code{scdensity}.
#'
#' Creates a plot of a shape-constrained kernel density estimate.  The amount of information
#' in the plot is controlled by \code{detail}.
#'
#' @param x An object of S3 class \code{scdensity}.
#' @param detail An integer from 1 to 4, indicating the level of information to include in
#' the plot.  1: plot only the constrained estimate.  2: draw both the constrained and
#' unconstrained estimates on the same plot.  3: add a rug showing the data points.
#' 4: additionally plot a Q-Q plot of the observed data versus the constrained estimate in a
#' second panel (for qualitative assessment of goodness-of-fit).
#' @param main A string passed on to the \code{main} argument of the plot command. If
#' \code{detail == 4}, pass a vector of two strings to specify titles for both subfigures.
#' @param xlab A string passed on to the \code{xlab} argument of the plot command. If
#' \code{detail == 4}, pass a vector of two strings to specify x labels for both subfigures.
#' @param ylab A string passed on to the \code{ylab} argument of the plot command. If
#' \code{detail == 4}, pass a vector of two strings to specify y labels for both subfigures.
#' @param type A vector of up to 3 strings specifying the \code{type} of plot used for 1)
#' the constrained estimate, 2) the unconstrained estimate, and 3) the Q-Q plot.
#' @param lty A vector of up to length 3, specifying the \code{lty} arguments passed to the
#' plot commands for 1) the constrained estimate, 2) the unconstrained estimate, and 3)
#' the Q-Q plot. See the description of \code{lty} in [graphics::par()].
#' @param pch A vector of up to 3 integers specifying the \code{pch} argument passed to the
#' plot commands for 1) the constrained estimate, 2) the unconstrained estimate, and 3)
#' the Q-Q plot.  See [graphics::points()] for the integer codes.
#' @param col A vector of up to 3 strings specifying the \code{col} argument passed to the
#' plot commands for 1) the constrained estimate, 2) the unconstrained estimate, and 3)
#' the Q-Q plot.
#' @param lwd A vector of up to length 3 specifying the \code{lwd} argument passed to the
#' plot commands for 1) the constrained estimate, 2) the unconstrained estimate, and 3)
#' the Q-Q plot.
#' @param zero.line A logical value indicating whether or not a horizontal line should be
#' drawn through zero to aid visualization.
#' @param ... Extra parameters passed to the initial \code{plot} command for each subfigure.
#'
#' @importFrom grDevices gray
#' @importFrom graphics abline layout legend lines plot rug
#'
#' @export
#'
#' @examples
#' # Basic usage:
#' x <- rlnorm(30)
#' scKDE <- scdensity(x)
#' plot(scKDE)
#'
#' # Show only the constrained estimate
#' plot(scKDE, detail=1)
#'
#' # Show the constrained and unconstrained estimates.  Change line color and width.
#' plot(scKDE, detail=2, col=c("red","blue"), lwd=c(3,2))
#'
#' # Show the Q-Q plot, but change that plot's symbol and its size.
#' plot(scKDE, detail=4, pch=c(-1, -1, 3), cex=0.5)
plot.scdensity <- function(x, detail = 4,
                           main = c("Density Estimate", "Q-Q Plot"),
                           xlab = c(x$data.name, "Constrained KDE Quantiles"),
                           ylab = c("Density", "Sample Quantiles"),
                           type = c("l","l","p"), lty=c(1,2,0), pch=c(-1,-1,1),
                           col=c("black",gray(0.4),"black"),
                           lwd=c(2,1,0), zero.line=TRUE, ...) {

  if (!is.element(detail, 1:4)) {
    detail <- 4
    warning("Argument 'detail' is invalid.  Setting detail = 4.")
  }

  if (detail>3) {
    layout(matrix(c(1,1,2),1,3,byrow=TRUE))
  } else {
    layout(1)
  }

  if (detail==1) {
    plot(x$x, x$y, type=type[1], lty=lty[1], pch=pch[1], col=col[1], lwd=lwd[1],
         main = main[1], xlab = xlab[1], ylab = ylab[1], ...)
  }

  if (detail>1) {
    yrng <- c(0, max(c(x$fhat(x$x), x$f0(x$x))))
    plot(x$x, x$y, type=type[1], lty=lty[1], pch=pch[1], col=col[1], lwd=lwd[1],
         ylim=yrng, main = main[1], xlab = xlab[1], ylab = ylab[1], ...)
    lines(x$x, x$f0(x$x), type=type[2], lty=lty[2], pch=pch[2], col=col[2], lwd=lwd[2])
    legend("topright", legend=c("constrained","unconstrained"),
           lty=lty[1:2], pch=pch[1:2], col=col[1:2], lwd=lwd[1:2])
  }

  if (detail>2) {
    rug(jitter(x$data), ticksize=0.02)
  }

  if (zero.line) {
    abline(h=0, col="gray")
  }

  if (detail>3) {
    n <- length(x$data)
    Q <- getQuantileFunction(x)
    q.theo <- Q((1:n - 0.5)/n)
    q.data <- sort(x$data)
    plot(q.theo, q.data, type=type[3], lty=lty[3], pch=pch[3], col=col[3], lwd=lwd[3],
         main = main[2], xlab = xlab[2], ylab = ylab[2], ...)
    qqline <- range(c(q.theo, q.data))
    lines(qqline, qqline, lty=2)
  }

}














