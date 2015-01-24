# Nested Ranks Test
#
# Calculate a Mann-Whitney-Wilcoxon (MWW) test using nested ranks.  Data are
# structured into several groups and each group has received two treatment
# levels. The rankings are compared between treatment levels, taking group
# structure and size into account when permuting ranks to construct a null
# distribution which assumes no influence of treatment.  When there is just one
# group, this test is identical to a standard MWW test.
#
# These statistical tools were developed in collaboration with Peter Smouse
# (Rutgers University) and Victoria Sork (UCLA) and were funded by U.S.
# National Science Foundation awards NSF-DEB-0514956 and NSF-DEB-0516529.
#
# Reference:
#
# Thompson PG, Smouse PE, Scofield DG, Sork VL. 2014.  What seeds tell us 
# about birds: a multi-year analysis of acorn woodpecker foraging movements.
# Movement Ecology 2:12. doi:10.1186/2051-3933-2-12.
#
# TODO: make sure year and granary load as factors
# TODO: doc "n.iter == 1 simply returns the observed Z-score"
# TODO: see documentation for wilcox.test for how to document both interfaces
#       at once
# TODO: \seealso{\code{\link{wilcox.test}}
# TODO: note also wilcox.test documentation for data=
# TODO: Complete package documentation using Roxygen2
# TODO: How can we generate *-package.Rd documentation using Roxygen2?
# TODO: Update documentation for each function
# TODO: Check documentation for plot.htest_boot, make sure params interleaved
#       with comments
# TODO: Set up testing structure
# TODO: Can the permutations in nestedRanksTest.default() be sped up?

#' Generic
#'
nestedRanksTest <- function(x, ...) UseMethod("nestedRanksTest")

nestedRanksTest.formula <- function(formula, data, groups = NULL, subset, ...)
{
  # initial version largely copied from stats:::t.test, extensively modified since
    if (missing(formula) || (length(formula) != 3L) || 
        (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame())))  # if data passed in as a matrix
        m$data <- as.data.frame(data)
    # Now, expand the formula, see if we have groups there via '|' or via groups=
    if (is.null(m$groups)) {
        if (length(formula[[3]]) == 3L && formula[[3]][[1]] == as.name("|")) {
            # remove grouping from formula, add it to groups=
            GROUP.NAME <- as.character(formula[[3]][[3]])
            m$groups = formula[[3]][[3]]
            formula[[3]] <- formula[[3]][[2]]
            m$formula <- formula
        } else stop("invalid group specification in formula")
    } else {
        # group structure specifed via groups=, could be name or vector
        if (length(formula[[3]]) == 3L && formula[[3]][[1]] == as.name("|"))
            stop("groups are specified with '|' in formula or groups=, but not both")
        GROUP.NAME <- deparse(substitute(groups))
    }
    m[[1L]] <- quote(stats::model.frame)  # this handles subset for us
    m$... <- NULL
    mf <- eval(m, parent.frame())  # mf contains our variables as y, x, "(groups)"
    if (! nrow(mf))
        stop("no data")
    Y.NAME <- names(mf)[attr(attr(mf, "terms"), "response")]
    X.NAME <- attr(attr(mf, "terms"), "term.labels")
    DNAME <- paste(Y.NAME, "by", X.NAME, "grouped by", GROUP.NAME)
    DATA <- list(y = mf[[Y.NAME]], x = factor(mf[[X.NAME]]), 
                 groups = factor(mf[["(groups)"]]))
    if (nlevels(DATA$x) != 2L)
        stop(X.NAME, " must have exactly 2 levels")
    # TODO: consider renaming these...
    #DATA <- list(y = mf[[response]], x = tmt, groups = mf[["(groups)"]])
    #y <- nestedRanksTest.default(y = TREATMENT.NAME, x = tmt, groups = mf[["(groups)"]], ...)
    y <- do.call("nestedRanksTest.default", c(DATA, list(...)))
    y$data.name <- DNAME
    return(y)
}


# must suffix with '.default' when a proper S3 class ???
#'
#' @keywords htest nonparametric
nestedRanksTest.default <- function(x, y, groups, n.iter = 10000, lightweight = FALSE)
{
    # TODO: deparse substitute with do.call
    if (n.iter < 1)
        stop("n.iter must be greater than or equal to 1")
    X.NAME = deparse(substitute(x))
    Y.NAME = deparse(substitute(y))
    GROUPS.NAME = deparse(substitute(groups))
    DNAME <- paste(Y.NAME, "by", X.NAME, "grouped by", GROUPS.NAME)
    METHOD <- "Nested Ranks Test"
    STATISTIC.NAME <- "Z"
    dat <- data.frame(y = y, x = factor(x), groups = factor(groups))
    nr <- nrow(dat)
    #  remove any entries with NA for y, x, or groups
    dat <- dat[apply(dat, 1, function(x) all(! is.na(x))), ]
    BAD.OBS <- nr - nrow(dat)
    x_levels <- levels(dat$x)
    if (length(x_levels) != 2) 
        stop(X.NAME, " requires exactly two levels for treatment")
    if (any(table(dat$groups, dat$x) == 0))
      stop(X.NAME, " must have values for all groups in both treatment levels")
    groups_df <- .nestedRanksTest_weights(dat$x, dat$groups)
    groups_levels <- row.names(groups_df)
    weights <- setNames(groups_df$weights, groups_levels)
    Z = matrix(0, n.iter, length(groups_levels), dimnames = list(NULL, groups_levels))
    for (i in seq_along(groups_levels)) {
        group.info <- groups_df[i, ]
        y1 <- dat$y[dat$groups == groups_levels[i] & dat$x == x_levels[1]]
        y2 <- dat$y[dat$groups == groups_levels[i] & dat$x == x_levels[2]]
        stopifnot(length(y1) == group.info$n1 && length(y2) == group.info$n2)
        y.vals <- c(y1, y2)
        this.Z <- numeric(n.iter)
        if (n.iter > 1)
            for (j in 1:(n.iter - 1))
                this.Z[j] <- .nestedRanksTest_MWW(sample(y.vals), group.info$n1, group.info$n2)
        this.Z[n.iter] <- .nestedRanksTest_MWW(y.vals, group.info$n1, group.info$n2)
        Z[, i] <- this.Z
    }
    null.distribution <- apply(Z, 1, function(x) sum(x * weights))
    stopifnot(length(null.distribution) == n.iter)
    quantiles <- c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)
    NULL.VALUES <- round(quantile(null.distribution, quantiles), ceiling(log10(n.iter)) + 1)
    STATISTIC <- setNames(null.distribution[n.iter], STATISTIC.NAME)
    PVAL <- sum(null.distribution >= STATISTIC) / n.iter
    RVAL <- list(statistic = STATISTIC,
                 p.value = PVAL,
                 alternative = paste(STATISTIC.NAME, "lies above bootstrapped null values"),
                 method = METHOD,
                 data.name = DNAME,
                 bad.obs = BAD.OBS,
                 null.values = NULL.VALUES,
                 n.iter = n.iter,
                 weights = weights)
    if (! lightweight)
        RVAL$null.distribution = null.distribution
    class(RVAL) <- c("htest_boot", "htest")
    return(RVAL)
}

#' Calculates Mann-Whitney-Wilcoxon Z-score
#'
#' \code{.nestedRanksTest_MWW} calculates the Mann-Whitney-Wilcoxon
#' Z-score for the ranks of \code{x} divided into two treatment
#' levels, the first \code{n1} of \code{x} and the final \code{n2}
#' of \code{x}.
#'
#' @param x    Values to be ranked for the test.  Its length must
#'             be equal to the sum of \code{n1} and \code{n2}.
#' @param n1   The first \code{n1} values in \code{x} belong to the
#'             first treatment level.
#' @param n2   The final \code{n2} values in \code{x} belong to the
#'             second treatment level.
#' 
#' @return The calculated Z-score
#'
#' @seealso \code{\link{nestedRanksTest}}, \code{\link{wilcox.test}}
.nestedRanksTest_MWW <- function(x, n1, n2) {
    stopifnot(length(x) == n1 + n2)
    r <- rank(x)
    r1 <- r[1:n1]
    r2 <- r[(n1 + 1):(n1 + n2)]
    R1 <- sum(r1)
    R2 <- sum(r2)
    n1.n2 <- n1 * n2
    U1 <- R1 - (((n1 + 1) * n1) / 2)
    U2 <- R2 - (((n2 + 1) * n2) / 2)
    Z.12 <- (U2 - U1) / n1.n2
    return(Z.12)
}


#' Calculates weights for \code{nestedRanksTest} based on group sizes
#'
#' \code{.ntestedRanksTest_weights} calculates weights for 
#' \code{nestedRanksTest} based on group sizes.  The number of group members
#' in each of the two treatment levels is determined (\code{n1} and \code{n2})
#' together with their product (\code{n1.n2}), and the group-specific weight
#' is calculated by dividing \code{n1.n2} by the sum of \code{n1.n2} for all
#' groups.
#'
#' @param  x      Treatments, coerced to factor.  Must contain two levels.
#' @param  groups Groups, coerced to factor, with elements in the same order
#'                as for \code{x}.
#'
#' @return \code{data.frame} containing weights and other information for
#'         each group: columns \code{group}, a factor of group names, also 
#'         used for row names; \code{n1}, \code{n2}, and \code{n1.n2} for 
#'         integer group sizes in the first and second treatment levels and 
#'         their product; and numeric \code{weights} for the calculated 
#'         weights.
#'
#' @seealso \code{\link{nestedRanksTest}}
#'
.nestedRanksTest_weights <- function(x, groups) {
    x_levels <- levels(x <- factor(x))
    w <- table(groups <- factor(groups), x)
    w <- data.frame(groups = row.names(w), n1 = w[, x_levels[1]], 
                    n2 = w[, x_levels[2]])
    row.names(w) <- w$groups
    w <- transform(w, groups = factor(groups), n1.n2 = n1 * n2)
    w$weights <- w$n1.n2 / sum(w$n1.n2)
    return(w)
}

#' Print result of \code{nestedRanksTest}.
#'
#' \code{print.htest_boot} prints the return value of 
#' \code{\link{nestedRanksTest}}, #' a list of class \code{"htest_boot"} 
#' which extends class \code{"htest"} by including group weights, the 
#' number of bootstrap iterations, and the complete null distribution.
#' The latter is not printed by this function; it may be visualised with
#' \code{\link{plot.htest_boot}}.
#'
#' @param  x      Value returned by \code{nestedRanksTest}
#' @param  \dots  Additional arguments passed to \code{print.htest}
#'
#' @return The value of x is returned invisibly.
#'
#' @examples
#  data(woodpecker_multiyear)
#' print(nestedRanksTest(Dist ~ Year | Granary, data = woodpecker_multiyear, 
#'                       subset = Species == "agrifolia"))
#'
#' @seealso \code{\link{nestedRanksTest}} for the test description, 
#'   \code{\link{plot.htest_boot}} for a graphical plot of test
#'   results, and \code{\link{print.htest}} for the print method of
#'   the base class.
#'
print.htest_boot <- function(x, ...) {
    NextMethod(x, ...)
    cat("bootstrap iterations =", x$n.iter, "\ngroup weights:\n")
    print(x$weights)
    invisible(x)
}

#' Diagnostic plot of result from nestedRanksTest.
#'
#' \code{plot.htest_boot} creates a diagnostic plot using the return value
#' of \code{\link{nestedRanksTest}}.
#'
#' @param  x       Value returned by \code{nestedRanksTest}
#' @param  breaks  The number of breaks to use when plotting the distribution,
#'                 the default is calculated from \code{n.iter} of the call to
#'                 \code{nestedRanksTest}.
#'
#' The following control the appearance of the null distribution, and are
#' passed on to \code{\link{hist}}.  Default values are available above.
#'
#' @param  col
#' @param  border
#' @param  main
#' @param  xlab
#' @param  ylab
#'
#' The following control the appearance of the verticle line indicating the
#' observed value, and are passed on to \code{\link{abline}}.  Default values 
#' are available above.
#'
#' @param  p.col   Passed as \code{col}
#' @param  p.lty   Passed as \code{lty}
#' @param  p.lwd   Passed as \code{lwd}
#'
#' @param  \dots   Additional arguments passed to \code{hist} and 
#'                 \code{abline} for plotting
#'
#' @examples
#  data(woodpecker_multiyear)
#' plot(nestedRanksTest(Dist ~ Year | Granary, data = woodpecker_multiyear, 
#'                      subset = Species == "lobata"))
#'
#' @seealso \code{\link{nestedRanksTest}} for test description and 
#'     \code{\link{print.htest_boot}} for printing test results. 
#'
plot.htest_boot <- function(x, breaks, col = "lightblue", border = NA, 
                            p.col = "red", p.lty = 2, p.lwd = 2,
                            main = paste0(x$method, ", ", x$data.name, "\n",
                                          names(x$statistic), " = ", 
                                          round(x$statistic, ceiling(log10(x$n.iter))), 
                                          ", P = ", x$p.value),
                            xlab = "Distribution of Z-scores", 
                            ylab = paste0("Frequency (out of ", x$n.iter, ")"),
                            ...) {
    if (missing(breaks))
        breaks <- min(x$n.iter / 50, 100)
    hist(x$null.distribution, breaks = breaks, col = col, border = border, 
        main = main, xlab = xlab, ylab = ylab, ...)
    abline(v = x$statistic, col = p.col, lty = p.lty, lwd = p.lwd, ...)
}

