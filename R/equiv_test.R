#' Equivalence and non-inferiority tests for one- and two-sample problems
#'
#' Implementation of uniformly most powerful invariant equivalence tests for one- and two-sample problems (paired and unpaired).
#' Also one-sided alternatives (non-inferiority and non-superiority tests) are supported.
#' Basically a variant of a t-test with (relaxed) null and alternative hypotheses exchanged.
#'
#' \code{equiv.test} is modelled after (and borrows code from) R's \code{t.test()} and is intended to work as similarly as possible.
#'
#' This functions implements uniformly most powerful invariant equivalence tests for one-sample and (paired
#' or unpaired) two-sample problems. Also supported are one-sided versions (so-called non-inferiority or
#' non-superiority tests).
#'
#' All tests are on standardized (differences of) means \eqn{theta}:
#' \deqn{theta = (mu_x - mu) / sigma}
#' for the one-sample case,
#' \deqn{theta = (mu_d - mu) / sigma_d}
#' for the paired two-sample case and
#' \deqn{theta = (mu_x - mu_y - mu) / sigma}
#' for the unpaired test, where \eqn{sigma} is the standard deviation of \eqn{x} and \eqn{y} and
#' \eqn{sigma_d} is the standard deviation of the differences. \eqn{mu} is a shift parameter that can
#' be used to compare against a known value in the one-sample case. \eqn{mu} should usually be zero for
#' two-sample problems.
#'
#' The null and alternative hypotheses in equivalence tests (\code{alternative = "two.sided"}) are
#' \deqn{H_0: theta <= -eps \qquad or \qquad theta >= eps}
#' vs
#' \deqn{H_1: -eps < theta < eps}
#'
#' Currently, only symmetric equivalence intervals \eqn{(-eps, eps)} are supported.
#'
#' In the non-inferority-case (\code{alternative = "greater"}) we test
#' \deqn{H_0: theta <= -eps}
#' vs
#' \deqn{H_1: theta > -eps}
#'
#' In the non-superiority-case (\code{alternative = "less"}) we test
#' \deqn{H_0: theta >= eps}
#' vs
#' \deqn{H_1: theta < eps}
#'
#' If \code{paired} is \code{TRUE} then both \code{x} and \code{y} must be specified and they must be the same length.
#' Missing values are silently removed (in pairs if \code{paired} is \code{TRUE}).
#'
#' The formula interface is only applicable for the two-sample tests.
#'
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param alternative	a character string specifying the alternative hypothesis, must be one of "\code{two.sided}" (default), "\code{greater}" or "\code{less}". You can specify just the initial letter.
#' @param eps a single strictly positive number giving the equivalence limits.
#' @param mu a number indicating the true value of the mean (or difference in means if you are performing a two sample test).
#' @param paired a logical indicating whether you want a paired equivalence test in the two-sample case.
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} is a numeric variable giving the data values and \code{rhs} a factor with two levels giving the corresponding groups.
#' @param data an optional matrix or data frame containing the variables in the formula \code{formula}. By default the variables are taken from \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be used.
#' @param na.action a function which indicates what should happen when the data contain NAs. Defaults to \code{getOption("na.action")}.
#' @param ... further arguments to be passed to or from methods.
#' @return A list with class \code{htest} containing the following components:
#'     \item{statistic}{the value of the t-statistic.}
#'     \item{parameter}{the degrees of freedom for the t-statistic.}
#'     \item{p.value}{the p-value for the test.}
#'     \item{estimate}{the plug-in estimate of the standardized mean (or mean difference), i.e. the empirical mean
#'     (or difference of empirical means) divided by the empirical standard deviation. Note that this estimate is not unbiaded.}
#'     \item{null.value}{non-equivalence limits, i.e. boundaries of null hypothesis}
#'     \item{alternative}{a character string describing the alternative hypothesis.}
#'     \item{method}{a character string indicating what type of equivalence test was performed.}
#'     \item{data.name}{a character string giving the name(s) of the data.}
#' @references Wellek, S. (2010). Testing Statistical Hypotheses of Equivalence and Noniferiority. Second edition.  Boca Raton: Chapman & Hall.
#'   (especially Chapters 5.3 and 6.1).
#' @examples
#'   # compare two feed from chickwts dataset
#'   data("chickwts")
#'   chickwts2 <- chickwts[chickwts$feed %in% c("linseed", "soybean"),]
#'   chickwts2$feed <- droplevels(chickwts2$feed)
#'
#'   # similar but cannot be shown to be equivalent up to 0.5 sigma at 0.05 level^
#'   plot(weight ~ feed, data = chickwts2)
#'   equiv.test(weight ~ feed, data = chickwts2, eps = 0.5)
equiv.test <- function(x, ...) UseMethod("equiv.test")

#' @describeIn equiv.test Default S3 method:
equiv.test.default <-
  function(x, y = NULL, alternative = c("two.sided", "less", "greater"),
           eps = 1, mu = 0, paired = FALSE,
           ...)
  {
    alternative <- match.arg(alternative)

    if(length(eps) > 1 | eps[1] < 0) {
      warning("Only symmetric intervals are supported currently. \n Using absolutes value of first component of eps.")
      eps <- abs(eps[1])
      # we should stop here if eps is effectively zero.
      }

    if(!missing(mu) && (length(mu) != 1 || is.na(mu)))
      stop("'mu' must be a single number")
    # if(!missing(conf.level) &&
    #    (length(conf.level) != 1 || !is.finite(conf.level) ||
    #     conf.level < 0 || conf.level > 1))
    #   stop("'conf.level' must be a single number between 0 and 1")
    if( !is.null(y) ) {
      dname <- paste(deparse(substitute(x)),"and",
                     deparse(substitute(y)))
      if(paired)
        xok <- yok <- complete.cases(x,y)
      else {
        yok <- !is.na(y)
        xok <- !is.na(x)
      }
      y <- y[yok]
    }
    else {
      dname <- deparse(substitute(x))
      if (paired) stop("'y' is missing for paired test")
      xok <- !is.na(x)
      yok <- NULL
    }
    x <- x[xok]
    if (paired) {
      x <- x - y
      y <- NULL
    }
    nx <- length(x)
    mx <- mean(x)
    vx <- var(x)
    if(is.null(y)) {
      if(nx < 2) stop("not enough 'x' observations")
      df <- nx - 1
      ncp <- sqrt(nx) * eps
      stderr <- sqrt(vx / nx)
      if(stderr < 10 *.Machine$double.eps * abs(mx))
        stop("data are essentially constant")
      tstat <- (mx - mu)/stderr
      d <- (mx - mu) / sqrt(vx)
      estimate <- setNames(d, "d")
      if (alternative == "less") { # non-superiority
        pval <- pt(tstat, df, ncp = ncp)
        alternative <- "non-superiority"
        equivint <- c(eps, Inf)
        names(equivint) <- c("lower", "upper")
        method <- if(paired) "Paired non-superiority test" else "One Sample non-superiority test"
      }
      else if (alternative == "greater") { # non-inferiority
        pval <- pt(tstat, df, ncp = -ncp, lower.tail = FALSE)
        alternative <- "non-inferiority"
        equivint <- c(-Inf, -eps)
        names(equivint) <- c("lower", "upper")
        method <- if(paired) "Paired non-inferiority test" else "One Sample non-inferiority test"
      }
      else {  # actual (two-sided) equivalence test, one sample or paired samples
        pval <- pf(tstat^2, df1 = 1, df2 = df, ncp = ncp^2)
        alternative <- "equivalence"
        equivint <- rbind(c(-Inf, -eps),
                          c(eps, Inf))
        dimnames(equivint)[[2]] <- c("lower", "upper")
        method <- if(paired) "Paired equivalence test" else "One Sample equivalence test"
      }
    } else { # here we are in the unpaired two-sample case
      ny <- length(y)
      if(nx < 2) # we need at least two observations each to estimate the variance
        stop("not enough 'x' observations")
      if(ny < 2)
        stop("not enough 'y' observations")
      my <- mean(y)
      vy <- var(y)
      df <- nx + ny - 2
      v <- (nx - 1) * vx + (ny - 1) * vy
      v <- v / df
        stderr <- sqrt(v * (1 / nx + 1 / ny))
      if(stderr < 10 *.Machine$double.eps * max(abs(mx), abs(my)))
        stop("data are essentially constant")
      tstat <- (mx - my - mu)/stderr  # same as formula (6.5) in Wellek (2010)
      d <- (mx - my - mu) / sqrt(v)
      estimate <- setNames(d, "d")
      ncp <- sqrt(nx * ny) * eps / sqrt(nx + ny)
      if (alternative == "less") { # non-superiority
        pval <- pt(tstat, df, ncp = ncp)
        alternative <- "non-superiority"
        equivint <- c(eps, Inf)
        names(equivint) <- c("lower", "upper")
        method <- "Two sample non-superiority test"
      }
      else if (alternative == "greater") { # non-inferiority
        pval <- pt(tstat, df, ncp = -ncp, lower.tail = FALSE)
        alternative <- "non-inferiority"
        equivint <- c(-Inf, -eps)
        names(equivint) <- c("lower", "upper")
        method <- "Two sample non-inferiority test"
      }
      else {  # actual (two-sided) equivalence test
        pval <- pf(tstat^2, df1 = 1, df2 = df, ncp = ncp^2)
        alternative <- "equivalence"
        equivint <- rbind(c(-Inf, -eps),
                          c(eps, Inf))
        dimnames(equivint)[[2]] <- c("lower", "upper")
        method <- "Two sample equivalence test"
      }
    }
    names(tstat) <- "t"
    params <- c(df, ncp)
    names(params) <- c("df", "ncp")
    rval <- list(statistic = tstat, parameter = params, p.value = pval,
                 estimate = estimate, null.value = equivint,
                 alternative = alternative,
                 method = method, data.name = dname)
    class(rval) <- "htest"
    return(rval)
  }

#' @describeIn equiv.test S3 method for class 'formula'
equiv.test.formula <-
  function(formula, data, subset, na.action, ...)
  {
    if(missing(formula)
       || (length(formula) != 3L)
       || (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
      stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
      m$data <- as.data.frame(data)
    ## need stats:: for non-standard evaluation
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if(nlevels(g) != 2L)
      stop("grouping factor must have exactly 2 levels")
    DATA <- setNames(split(mf[[response]], g), c("x", "y"))
    y <- do.call("equiv.test", c(DATA, list(...)))
    y$data.name <- DNAME
    if(length(y$estimate) == 2L)
      names(y$estimate) <- paste("mean in group", levels(g))
    y
  }
