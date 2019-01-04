# equiv.test is modelled after R's t.test() and intended to work as similarly as possible.

equiv.test <- function(x, ...) UseMethod("equiv.test")

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
      x <- x-y
      y <- NULL
    }
    nx <- length(x)
    mx <- mean(x)
    vx <- var(x)
    if(is.null(y)) {
      if(nx < 2) stop("not enough 'x' observations")
      df <- nx-1
      stderr <- sqrt(vx/nx)
      if(stderr < 10 *.Machine$double.eps * abs(mx))
        stop("data are essentially constant")
      tstat <- (mx-mu)/stderr
      method <- if(paired) "Paired equivalence test" else "One Sample equivalence test"
      estimate <-
        setNames(mx, if(paired)"mean of the differences" else "mean of x")
      if (alternative == "less") { # non-superiority, not treated yet
        pval <- pt(tstat, df)
        stop("non-superiority and non-inferiority test not implemented yet")
        }
      else if (alternative == "greater") { # non-inferiority, not treated yet
        pval <- pt(tstat, df, lower.tail = FALSE)
        stop("non-superiority and non-inferiority test not implemented yet")
      }
      else {  # actual (two-sided) equivalence test, on sample or paired samples
        ncp <- nx * eps^2
        pval <- pf(tstat^2, df1 = 1, df2 = df, ncp = ncp)
        alternative <- "equivalence"
        equivint <- c(mu-eps, mu + eps)
        names(equivint) <- c("lower", "upper")
      }
    } else { # here we are in the unpaired two-sample case
      ny <- length(y)
      if(nx < 2) # we need at least two observations each to estimate the variance
        stop("not enough 'x' observations")
      if(ny < 2)
        stop("not enough 'y' observations")
      my <- mean(y)
      vy <- var(y)
      method <- "Two sample equivalence test" # only true if alternative is two.sided
      estimate <- c(mx,my)
      names(estimate) <- c("mean of x","mean of y")
      df <- nx+ny-2
      v <- (nx-1)*vx + (ny-1)*vy
      v <- v/df
        stderr <- sqrt(v*(1/nx+1/ny))
      if(stderr < 10 *.Machine$double.eps * max(abs(mx), abs(my)))
        stop("data are essentially constant")
      tstat <- (mx - my - mu)/stderr
      if (alternative == "less") { # non-superiority, not treated yet
        pval <- pt(tstat, df)
        stop("non-superiority and non-inferiority test not implemented yet")
      }
      else if (alternative == "greater") { # non-inferiority, not treated yet
        pval <- pt(tstat, df, lower.tail = FALSE)
        stop("non-superiority and non-inferiority test not implemented yet")
      }
      else {  # actual (two-sided) equivalence test
        ncp <- nx * ny * eps^2 / (nx + ny)
        pval <- pf(tstat^2, df1 = 1, df2 = df, ncp = ncp)
        alternative <- "equivalence"
        equivint <- c(mu-eps, mu + eps)
        names(equivint) <- c("lower", "upper")
      }
    }
    names(tstat) <- "t"
    names(df) <- "df"
#    names(mu) <- if(paired || !is.null(y)) "difference in means" else "mean"
    rval <- list(statistic = tstat, parameter = df, p.value = pval,
                 estimate = estimate, null.value = equivint,
                 alternative = alternative,
                 method = method, data.name = dname)
    class(rval) <- "htest" # TO DO
    return(rval)
  }

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