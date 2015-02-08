#' Mann-Whitney-Wilcoxon ranks test when data are in groups.
#'
#' Calculate a Mann-Whitney-Wilcoxon test for a difference between treatment
#' levels using nested ranks.  This test can be used when observations are
#' structured into several groups and each group has received both treatment
#' levels.  The p-value is determined via bootstrapping.  This test is intended
#' to be analogous to a mixed-model extension of the \code{\link{wilcox.test}},
#' for which treatment is a fixed effect and group membership is a random
#' effect.
#'
#' The main function is \code{\link{nestedRanksTest}}, which includes a
#' formula interface implementing the familiar \code{"|"} syntax for
#' specifying group membership on the right-hand side of the formula.
#' The value returned is a list of class \code{'htest_boot'}, which
#' extends class \code{'htest'}.  \code{print} and \code{plot} methods
#' are provided to print and visualise results.
#'
#' These statistical tools were developed in collaboration with Peter E.
#' Smouse (Rutgers University) and Victoria L. Sork (UCLA) and were funded by
#' U.S. National Science Foundation awards NSF-DEB-0514956 and
#' NSF-DEB-0516529.
#'
#' @references
#' Thompson, P. G., Smouse, P. E., Scofield, D. G. and Sork, V. L. (2014)
#' What seeds tell us about birds: a multi-year analysis of acorn woodpecker
#' foraging movements.  \emph{Movement Ecology} 2:12.
#' \url{http://www.movementecologyjournal.com/content/2/1/12}
#'
#' \url{https://github.com/douglasgscofield/nestedRanksTest}
#'
#' @docType package
#'
#' @name nestedRanksTest-package
#'
NULL



#' Mann-Whitney-Wilcoxon ranks test when data are in groups.
#'
#' The statistic for the nested ranks test is a Z-score calculated by
#' comparing ranks between treatment levels, with contributions of each group
#' to the final Z-score weighted by group size.  The p-value is determined by
#' comparing the observed Z-score against a distribution of Z-scores
#' calculated by bootstrapping ranks assuming no influence of treatment while
#' respecting group sizes. When there is just one group, this test is
#' essentially identical to a standard Mann-Whitney-Wilcoxon test.  This test
#' is intended to be a mixed-model extension of the \code{\link{wilcox.test}},
#' for which treatment is a fixed effect and group membership is a random
#' effect.
#'
#' @note Cases for which any of \code{x}, \code{y} or \code{groups} are
#'       \code{NA} are removed.
#'
#' @note The generation of a null distribution can take some time.  For
#'       example, if any use of \code{nestedRanksTest} in the examples were
#'       run with the default \code{n.iter = 10000}, completion would require
#'       a few seconds.
#'
#' @param x       A (non-empty) vector of treatments for each \code{y},
#'                coerced to factor.  Must contain exactly two levels.
#' @param y       A (non-empty) numeric vector of data values.
#' @param groups  A (non-empty) vector specifying group membership for each
#'                \code{y}, coerced to a factor.  There must be at least one
#'                \code{y} in each group for each treatment level.
#' @param formula A formula of the form \code{lhs ~ rhs} or
#'                \code{lhs ~ rhs | groups}, where \code{lhs}
#'                is a numeric variable giving the data values, \code{rhs}
#'                is a variable obeying conditions for \code{x}, and
#'                \code{groups} is a variable obeying conditions for
#'                \code{groups}.  If \code{"| groups"} is not included
#'                in the formula, group membership must be specified with
#'                the \code{groups} argument.
#' @param data    An optional matrix or data frame (or similar: see
#'                \code{\link{model.frame}}) containing the variables in the
#'                formula \code{formula}.  By default the variables are taken
#'                from \code{environment(formula)}.
#' @param subset  An optional vector specifying a subset of observations to
#'                be used.
#' @param n.iter  Number of bootstrap iterations to perform.  The value of
#'                the final iteration is provided by the observed Z-score.
#'                Using \code{n.iter = 1} simply returns the observed Z-score.
#' @param lightweight  If \code{TRUE}, the vector of individual values of
#'                the null distribution is excluded from the return value of
#'                class \code{'htest_boot'}.  By default the null distribution
#'                is included. If \code{n.iter} is large, specifying
#'                \code{TRUE} for this option can save space, but note that
#'                calling \code{plot} on the return value will produce an
#'                error if so.
#' @param ...     Further arguments to be passed to or from methods.
#'
#' @return A list with class \code{'htest_boot'} based on class
#'         \code{'htest'} containing the following components.  Components
#'         marked with \code{"*"} are additions to \code{'htest'}.
#' \tabular{ll}{
#'     \code{statistic}           \tab the value of the observed Z-score.\cr
#'     \code{p.value}             \tab the p-value for the test.\cr
#'     \code{alternative}         \tab a character string describing the
#'                                     alternative hypothesis.\cr
#'     \code{method}              \tab a character string indicating the nested
#'                                     ranks test performed.\cr
#'     \code{data.name}           \tab a character string giving the name(s) of
#'                                     the data.\cr
#'     \code{bad.obs}             \tab the number of observations in the data
#'                                     excluded because of \code{NA} values.\cr
#'     \code{null.values}         \tab quantiles of the null distribution used
#'                                     for calculating the p-value.\cr
#'     \code{n.iter*}             \tab the number of bootstrap iterations used
#'                                     for generating the null distribution.\cr
#'     \code{weights*}            \tab the weights for groups, calculated by
#'                                     \code{nestedRanksTest_weights}.\cr
#'     \code{null.distribution*}  \tab null distribution of Z-scores, with
#'                                     \code{statistic} the last value.\cr
#' }
#' The length of \code{null.distribution} equals \code{n.iter}.  Note that
#' \code{null.distribution} will not be present if the
#' \code{lightweight = TRUE} option was given to \code{nestedRanksTest}.
#'
#'
#' @examples
#' require(graphics)
#'
#' data(woodpecker_multiyear)
#'
#' ## S3 method for class 'formula'
#'
#' ## n.iter set to 1000 to shorten completion time
#'
#' ## group in formula
#' nestedRanksTest(Distance ~ Year | Granary, n.iter = 1000,
#'                 data = woodpecker_multiyear,
#'                 subset = Species == "agrifolia")
#' ## group in 'groups='
#' nestedRanksTest(Distance ~ Year, groups = Granary, n.iter = 1000,
#'                 data = woodpecker_multiyear,
#'                 subset = Species == "lobata")
#'
#' ## Default S3 method
#'
#' dat.a <- subset(woodpecker_multiyear, Species == "agrifolia")
#' ## arguments in default order
#' nestedRanksTest(dat.a$Year, dat.a$Distance, dat.a$Granary, n.iter = 1000)
#' ## named arguments used in 'formula' order
#' res <- with(subset(woodpecker_multiyear, Species == "lobata"),
#'            nestedRanksTest(y = Distance, x = Year, groups = Granary,
#'                            n.iter = 1000))
#' plot(res)
#'
#' @seealso \code{\link{wilcox.test}}, \code{\link{print.htest_boot}},
#'          \code{\link{plot.htest_boot}}
#'
#' @references
#' Thompson, P. G., Smouse, P. E., Scofield, D. G. and Sork, V. L. (2014)
#' What seeds tell us about birds: a multi-year analysis of acorn woodpecker
#' foraging movements.  \emph{Movement Ecology} 2:12.
#' \url{http://www.movementecologyjournal.com/content/2/1/12}
#'
#' \url{https://github.com/douglasgscofield/nestedRanksTest}
#'
#' @keywords htest nonparametric
#'
#' @export
#'
#' @name nestedRanksTest
NULL

# Don't document this, just the methods
nestedRanksTest <- function(x, ...) UseMethod("nestedRanksTest")



#' @rdname nestedRanksTest
#'
#' @export
#'
nestedRanksTest.formula <- function(formula, data, groups = NULL, subset, ...) {
    # initial inspiration drawn from stats:::t.test.formula, trying to follow
    # its basic features, structure and naming conventions for both formula and
    # default methods
    if (missing(formula) || (length(formula) != 3L) ||
        (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    # Expand formula, see if groups via '|' or groups=
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
            stop("groups are specified with '|' in formula or with groups=",
                 " argument, but not both")
        GROUP.NAME <- deparse(substitute(groups))
    }
    m[[1L]] <- quote(stats::model.frame)  # handles subset= for us
    m$... <- NULL
    mf <- eval(m, parent.frame())  # columns in returned mf: y x "(groups)"
    if (! nrow(mf))
        stop("no data")
    Y.NAME <- names(mf)[attr(attr(mf, "terms"), "response")]
    X.NAME <- attr(attr(mf, "terms"), "term.labels")
    DNAME <- paste(Y.NAME, "by", X.NAME, "grouped by", GROUP.NAME)
    DATA <- list(y = mf[[Y.NAME]], x = factor(mf[[X.NAME]]),
                 groups = factor(mf[["(groups)"]]))
    if (nlevels(DATA$x) != 2L)
        stop(X.NAME, " must have exactly 2 levels")
    y <- do.call("nestedRanksTest.default", c(DATA, list(...)))
    y$data.name <- DNAME
    return(y)
}




#' @rdname nestedRanksTest
#'
#' @export
#'
nestedRanksTest.default <- function(x, y, groups, n.iter = 10000,
                                    lightweight = FALSE, ...) {
    if (n.iter < 1)
        stop("n.iter must be greater than or equal to 1")
    if (missing(groups))
        stop("'groups' missing")
    X.NAME = deparse(substitute(x))
    Y.NAME = deparse(substitute(y))
    GROUPS.NAME = deparse(substitute(groups))
    DNAME <- paste(Y.NAME, "by", X.NAME, "grouped by", GROUPS.NAME)
    METHOD <- "Nested Ranks Test"
    STATISTIC.NAME <- "Z"
    dat <- data.frame(y = y, x = factor(x), groups = factor(groups))
    nr <- nrow(dat)
    # remove any entries with NA for y, x, or groups, and refactor
    dat <- dat[apply(dat, 1, function(x) all(! is.na(x))), ]
    dat <- transform(dat, x = factor(x), groups = factor(groups))
    BAD.OBS <- nr - nrow(dat)
    # In error messages, we have to use 'x' here because if this
    # is called from the formula method, deparse(substitute(x)) gives
    # us values rather than the name, and I haven't figured out a
    # clean way around that.
    if (nlevels(dat$x) != 2)
        stop("'x' requires exactly two treatment levels")
    if (any(table(dat$groups, dat$x) == 0))
        stop("'x' must have values for all groups in both treatment levels")
    groups_df <- nestedRanksTest_weights(dat$x, dat$groups)
    groups_levels <- row.names(groups_df)
    weights <- setNames(groups_df$weights, groups_levels)
    Z = matrix(0, n.iter, length(groups_levels),
               dimnames = list(NULL, groups_levels))
    x_levels <- levels(dat$x)
    for (i in seq_along(groups_levels)) {
        group.info <- groups_df[i, ]
        y1 <- dat$y[dat$groups == groups_levels[i] & dat$x == x_levels[1]]
        y2 <- dat$y[dat$groups == groups_levels[i] & dat$x == x_levels[2]]
        stopifnot(length(y1) == group.info$n1 && length(y2) == group.info$n2)
        y.vals <- c(y1, y2)
        this.Z <- numeric(n.iter)
        if (n.iter > 1)
            for (j in 1:(n.iter - 1))
                this.Z[j] <- nestedRanksTest_Z(sample(y.vals, replace = FALSE),
                                               group.info$n1, group.info$n2)
        this.Z[n.iter] <- nestedRanksTest_Z(y.vals, group.info$n1,
                                            group.info$n2)
        Z[, i] <- this.Z
    }
    null.distribution <- apply(Z, 1, function(x) sum(x * weights))
    stopifnot(length(null.distribution) == n.iter)
    quantiles <- c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)
    NULL.VALUES <- round(quantile(null.distribution, quantiles),
                         max(ceiling(log10(n.iter)) + 1, 3))
    STATISTIC <- setNames(null.distribution[n.iter], STATISTIC.NAME)
    PVAL <- sum(null.distribution >= STATISTIC) / n.iter
    RVAL <- list(statistic   = STATISTIC,
                 p.value     = PVAL,
                 alternative = paste(STATISTIC.NAME,
                                     "lies above bootstrapped null values"),
                 method      = METHOD,
                 data.name   = DNAME,
                 bad.obs     = BAD.OBS,
                 null.values = NULL.VALUES,
                 n.iter      = n.iter,
                 weights     = weights)
    if (! lightweight)
        RVAL$null.distribution = null.distribution
    class(RVAL) <- c('htest_boot', 'htest')
    return(RVAL)
}




#' Calculates Z-score from ranks.
#'
#' \code{nestedRanksTest_Z} is used by \code{nestedRanksTest} to
#' calculate the Z-score for the ranks of responses \code{y} divided
#' into two treatment levels.
#'
#' Values across both treatments are ranked using the base R function
#' \code{rank} with \code{ties.method = "average"}, which assigns
#' tied values their average rank.  The Mann-Whitney-Wilcoxon test
#' statistic is computed from these ranks.  Because the value of the
#' statistic is sample-size dependent (between \code{-n1*n2} and
#' \code{n1*n2}), it is scaled to be \code{[-1,+1]} by dividing by
#' \code{n1*n2}.
#'
#' The bottleneck for bootstrapping is calculation of ranks, so the most
#' straightforward way to speed up \code{nestedRanksTest} would come from
#' speeding up \code{rank}.  Because of the checks performed prior to
#' calling this routine, it should be sufficient to use a stripped-down
#' function that simply does the equivalent of making an \code{.Internal}
#' call, which is not allowed within package code.  As of this writing, this
#' is sufficient:
#'
#' \code{rank_new <- function (x) .Internal(rank(x, length(x), "average"))}
#'
#' For the example data this is 8-9 times faster than the base R \code{rank},
#' because it avoids error-checking overhead.  For longer vectors, the
#' advantage decreases such that at 10000 elements it is 20-30\%.
#'
#' @param y    Values to be ranked for the test.  Its length must
#'             be equal to the sum of \code{n1} and \code{n2}.
#' @param n1   The first \code{n1} values in \code{y} belong to the
#'             first treatment level.
#' @param n2   The final \code{n2} values in \code{y} belong to the
#'             second treatment level.
#'
#' @return The calculated Z-score
#'
#' @seealso \code{\link{nestedRanksTest}}, \code{\link{wilcox.test}}
#'
#' @export
#'
nestedRanksTest_Z <- function(y, n1, n2) {
    r <- rank(y, ties.method = "average")
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



#' Calculates weights for \code{nestedRanksTest} based on group sizes.
#'
#' \code{ntestedRanksTest_weights} is used by \code{nestedRanksTest} to
#' calculate group weights based on group sizes.  The number of group members
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
#' @export
#'
nestedRanksTest_weights <- function(x, groups) {
    x_levels <- levels(x <- factor(x))
    w <- table(groups <- factor(groups), x)
    w <- data.frame(groups = row.names(w), n1 = w[, x_levels[1]],
                    n2 = w[, x_levels[2]])
    row.names(w) <- w$groups
    w$groups <- factor(w$groups)
    w$n1.n2 <- w$n1 * w$n2
    w$weights <- w$n1.n2 / sum(w$n1.n2)
    return(w)
}



#' Print result of \code{nestedRanksTest}.
#'
#' \code{print.htest_boot} prints the return value of
#' \code{\link{nestedRanksTest}}, a list of class \code{'htest_boot'}
#' which extends class \code{'htest'} by including group weights, the
#' number of bootstrap iterations, and the complete null distribution.
#' The latter is not printed by this function; it may be visualised with
#' \code{\link{plot.htest_boot}}.
#'
#' @param  x      Value of class \code{'htest_boot'} as returned by
#'                \code{nestedRanksTest}.
#' @param  \dots  Additional arguments passed to \code{print.htest}.
#'
#' @return The value of x is returned invisibly.
#'
#' @examples
#' data(woodpecker_multiyear)
#' ## n.iter set to 1000 to shorten completion time
#' res <- nestedRanksTest(Distance ~ Year | Granary, n.iter = 1000,
#'                        data = woodpecker_multiyear,
#'                        subset = Species == "agrifolia")
#' class(res)
#' print(res)
#'
#' @seealso \code{\link{nestedRanksTest}} for the test description,
#'   \code{\link{plot.htest_boot}} for a graphical plot of test
#'   results, and \code{\link{print.htest}} for the print method of
#'   the base class.
#'
#' @export
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
#' of \code{\link{nestedRanksTest}}, a list of class \code{'htest_boot'}.
#' The plot contains a histogram of the null distribution generated by
#' bootstrapping plotted with \code{\link{hist}}, and a verticle line
#' indicating the observed value plotted with \code{\link{abline}}.
#'
#' If there is no null distribution included in the class, because the
#' options \code{lightweight = TRUE} or \code{n.iter = 1} were given to
#' \code{nestedRanksTest}, this function produces an error.
#'
#' @param  x       Value of class \code{'htest_boot'} as returned by
#'                 \code{nestedRanksTest}.
#' @param  breaks  The number of breaks to use when plotting the distribution,
#'                 the default is calculated from \code{n.iter} of the call to
#'                 \code{nestedRanksTest}.
#' @param  col     Fill color for histogram bars, passed to \code{hist}.
#' @param  border  Border color for histogram bars, passed to \code{hist}.
#' @param  main    Main title, passed to \code{hist}.
#' @param  xlab    X-axis label, passed to \code{hist}.
#' @param  ylab    Y-axis label, passed to \code{hist}.
#' @param  p.col   Observed value line colour, passed to \code{abline}.
#' @param  p.lty   Observed value line type, passed to \code{abline}.
#' @param  p.lwd   Observed value line width, passed to \code{abline}.
#' @param  \dots   Additional arguments passed to \code{hist} and
#'                 \code{abline} for plotting.
#'
#' @return None
#'
#' @examples
#' require(graphics)
#'
#' data(woodpecker_multiyear)
#'
#' ## n.iter set to 1000 to shorten completion time
#' res.a <- nestedRanksTest(Distance ~ Year | Granary, n.iter = 1000,
#'                          data = woodpecker_multiyear,
#'                          subset = Species == "agrifolia")
#' res.l <- nestedRanksTest(Distance ~ Year | Granary, n.iter = 1000,
#'                          data = woodpecker_multiyear,
#'                          subset = Species == "lobata")
#'
#' opa = par(mfrow = c(2, 1))
#' ## Defaults
#' plot(res.l)
#' ## Modify colours, line type and main title
#' plot(res.a, main = "Quercus agrifolia", col = "lightgreen",
#'      p.col = "brown4", p.lty = 1)
#' par(opa)
#'
#' @seealso \code{\link{nestedRanksTest}} for test description,
#'     \code{\link{print.htest_boot}} for printing test results,
#'     and \code{\link{hist}} and \code{\link{abline}} for plotting options.
#'
#' @keywords hplot
#'
#' @export
#'
plot.htest_boot <- function(x, breaks, col = "lightblue", border = NA,
                            main = paste(sep = "",
                                         x$method, ", ", x$data.name, "\n",
                                         names(x$statistic), " = ",
                                         round(x$statistic,
                                               ceiling(log10(x$n.iter))),
                                         ", P = ", x$p.value),
                            xlab = "Distribution of Z-scores",
                            ylab = paste(sep = "",
                                         "Frequency (out of ", x$n.iter, ")"),
                            p.col = "red", p.lty = 2, p.lwd = 2,
                            ...) {
    if (is.null(x$null.distribution) || x$n.iter == 1)
        stop(deparse(substitute(x)), " does not contain a null distribution")
    if (missing(breaks))
        breaks <- min(ceiling(x$n.iter / 50), 100)
    hist(x$null.distribution, breaks = breaks, col = col, border = border,
        main = main, xlab = xlab, ylab = ylab, ...)
    abline(v = x$statistic, col = p.col, lty = p.lty, lwd = p.lwd, ...)
}



#' Distances acorns of two oak species were carried by acorn woodpeckers in
#' two different years.
#'
#' A dataset containing distances acorns of two oak species were carried by
#' acorn woodpeckers (\emph{Melanerpes formicivorus}) to their granaries, in
#' two different years for each oak species.  Data were collected in oak
#' savanna habitat in central California.  Acorn woodpeckers store acorns in
#' central granaries, and different woodpecker social groups maintain
#' different granaries. The variables are as follows:
#' \itemize{
#'     \item Species, the species of oak for the observed acorn ("lobata" for
#'         \emph{Quercus lobata}, "agrifolia" for \emph{Quercus agrifolia})
#'     \item Year, the year of observation (2002 and 2004 for
#'         \emph{Quercus lobata}, 2006 and 2007 for \emph{Quercus agrifolia})
#'     \item Granary, the woodpecker granary from which the acorn was collected
#'     \item Distance, distance in metres from the acorn source tree to the
#'         granary
#' }
#'
#' @references
#' Thompson, P. G., Smouse, P. E., Scofield, D. G. and Sork, V. L. (2014)
#' What seeds tell us about birds: a multi-year analysis of acorn woodpecker
#' foraging movements.  \emph{Movement Ecology} 2:12.
#' \url{http://www.movementecologyjournal.com/content/2/1/12}
#'
#' @format   Data frame with 534 rows and 4 variables.
#' @author   Douglas G. Scofield \email{douglasgscofield@@gmail.com}
#' @source   Dataset originates from the lab of Victoria L. Sork
#'           \email{vlsork@@ucla.edu} and is used with permission.
#' @docType  data
#' @name     woodpecker_multiyear
#'
NULL

