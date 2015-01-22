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
# TODO: General package documentation, can we use Roxygen?
# TODO: Update documentation for each function, including interface following
#       the Roxygen stuff from dev_tools
# TODO: Figure out how to register the *.formula and *.default methods
# TODO: Do the *MWW and *weights functions need to remain internal?
# TODO: Set up testing structure
# TODO: Can the permutations in nestedRanksTest.default() be sped up?


nestedRanksTest.formula = function(formula, data, groups = NULL, subset, ...)
{
  # initial version largely copied from stats:::t.test, extensively modified since
  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]), 
                                                                  "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))  # if data passed in as a matrix
    m$data <- as.data.frame(data)
  # Now, expand the formula, see if we have groups there via '|' or via groups=
  if (is.null(m$groups)) {
    if (length(formula[[3]]) == 3L && formula[[3]][[1]] == as.name("|")) { # grouping in formula
      # remove grouping from formula, add it to groups=
      GROUP.NAME <- as.character(formula[[3]][[3]])
      m$groups = formula[[3]][[3]]
      TREATMENT.NAME <- as.character(formula[[3]][[2]])
      formula[[3]] <- formula[[3]][[2]]
      m$formula <- formula
    } else stop("invalid group specification in formula")
  } else {  # group structure specifed via groups=, could be name or vector
    if (length(formula[[3]]) == 3L && formula[[3]][[1]] == as.name("|"))
      stop("groups are specified with '|' in formula or groups=, but not both")
    GROUP.NAME <- deparse(substitute(groups))
  }
  m[[1L]] <- quote(stats::model.frame)  # this handles subset for us
  m$... <- NULL
  mf <- eval(m, parent.frame())  # for (y ~ x, groups=G) mf now y , x , "(groups)"
  if (! nrow(mf))
    stop("no data")
  DNAME <- paste(names(mf)[1:2], collapse = " by ")
  DNAME <- paste(DNAME, GROUP.NAME, sep = " grouped by ")
  response <- attr(attr(mf, "terms"), "response")
  TREATMENT.NAME <- attr(attr(mf, "terms"), "term.labels")
  tmt <- factor(mf[[TREATMENT.NAME]])
  if (nlevels(tmt) != 2L)
    stop("treatment factor must have exactly 2 levels")
  # TODO: consider renaming these...
  DATA <- data.frame(y=mf[[response]], x=tmt, groups=mf[["(groups)"]])
  y <- do.call("nestedRanksTest.default", 
               c(list(y=mf[[response]], x=tmt, groups=mf[["(groups)"]]), list(...)))
  y$data.name <- DNAME
  return(y)
}


# must suffix with '.default' when a proper S3 class ???
nestedRanksTest.default = function(y, x, groups, n.iter=10000)
{
  Y.NAME = deparse(substitute(y))
  X.NAME = deparse(substitute(x))
  GROUPS.NAME = deparse(substitute(groups))
  DNAME <- paste(Y.NAME, "by", X.NAME, "grouped by", GROUPS.NAME)
  METHOD = "Nested Ranks Test"
  STATISTIC.NAME = "Z.weighted.obs"
  dat = data.frame(y=y, x=x, groups=groups)
  nr = nrow(dat)
  dat = dat[apply(dat[,1:3], 1, function(x) all(!is.na(x))), 1:3]
  BAD.OBS = nr - nrow(dat)
  vals = dat[, 1]
  tmt = factor(dat[, 2])
  grp = dat[, 3]
  tmt_levels = unique(sort(tmt))
  if (length(tmt_levels) != 2) 
    stop(GROUPS.NAME, "requires exactly two levels for treatment")
  s = split(grp, tmt)
  if (length(unique(sort(grp))) != length(intersect(s[[1]], s[[2]])))
    stop(X.NAME, "must have values for all groups in both treatment levels")
  wt = .nestedRanksTest_weights(dat)
  Grps = unique(sort(as.character(grp)))
  nm.Grps = setNames(make.names(Grps, unique = TRUE), Grps)
  weights = setNames(wt$Rel_Wt, rownames(wt))
  p = list()
  for (g in Grps) {
    tmt1.dat = dat[grp == g & tmt == tmt_levels[1], ]
    tmt2.dat = dat[grp == g & tmt == tmt_levels[2], ]
    n1 = nrow(tmt1.dat)
    n2 = nrow(tmt2.dat)
    this.vals = c(tmt1.dat[, 1], tmt2.dat[, 1])
    this.z = .nestedRanksTest_MWW(this.vals, n1, n2)
    Z = numeric(n.iter)
    if (n.iter > 1) {
      for (i in 1:(n.iter - 1)) {
        d = sample(this.vals)
        Z[i] = .nestedRanksTest_MWW(d, n1, n2)
      }
    }
    Z[n.iter] = this.z
    p[[ nm.Grps[g] ]] = Z
  }
  p = as.data.frame(p)
  if (! all(names(p) == names(weights)))
    stop("weights out of order")
  null.distribution <- apply(p, 1, function(x) sum(x * weights))
  quantiles <- c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1)
  NULL.VALUES <- quantile(null.distribution, quantiles)
  STATISTIC = setNames(null.distribution[n.iter], STATISTIC.NAME)
  PVAL = sum(null.distribution >= STATISTIC) / n.iter
  RVAL = list(statistic = STATISTIC,
              p.value = PVAL,
              alternative = paste(STATISTIC.NAME, "lies above bootstrapped values"),
              method = METHOD,
              data.name = DNAME,
              bad.obs = BAD.OBS,
              null.values = NULL.VALUES,
              n.iter = n.iter,
              weights = weights,
              null.distribution = null.distribution)
  class(RVAL) = c("htest_boot", "htest")
  return(RVAL)
}

# .nestedRanksTest_MWW calculates the Mann-Whitney-Wilcoxon Z-score
#    x  : values
#    n1 : the first n1 of x belong to the first level
#    n2 : the final n2 of x belong to the second level
# The return value is the calculated Z-score
.nestedRanksTest_MWW = function(x, n1, n2) {
  r = rank(x)
  r1 = r[1:n1]
  r2 = r[(n1+1):(n1+n2)]
  R1 = sum(r1)
  R2 = sum(r2)
  n1.n2 = n1 * n2
  U1 = (R1 - (((n1+1)*n1)/2))
  U2 = (R2 - (((n2+1)*n2)/2))
  Z.12 = (U2 - U1) / n1.n2
  return(Z.12)
}


# MWW.weights calculates sample-size weights
.nestedRanksTest_weights = function(dat) {
  vals = dat[, 1]
  tmt = dat[, 2]
  grp = dat[, 3]
  tmt_levels = unique(sort(tmt))
  Grps = unique(sort(as.character(grp)))
  w = data.frame()
  for (g in Grps) {
    #gdat = subset(dat, grp == g)
    #tmt1.dat = subset(gdat, treatment == tmt_levels[1])
    #tmt2.dat = subset(gdat, treatment == tmt_levels[2])
    tmt1.dat = dat[grp == g & tmt == tmt_levels[1], ]
    tmt2.dat = dat[grp == g & tmt == tmt_levels[2], ]
    n1 = nrow(tmt1.dat)
    n2 = nrow(tmt2.dat)
    n1.n2 = n1 * n2
    w = rbind(w, list(group = g,
                      n1    = n1,
                      n2    = n2,
                      n1.n2 = n1.n2))
  }
  w = transform(w, Rel_Wt = n1.n2 / sum(n1.n2))
  rownames(w) = make.names(w$group, unique = TRUE)
  return(w)
}

#
print.htest_boot <- function(x, ...) {
  NextMethod(x, ...)
  cat("bootstrap iterations =", x$n.iter, "\ngroup weights:\n")
  print(x$weights, ...)
  invisible(x)
}

# plot.htest_boot creates a utility plot of class htest_boot results
plot.htest_boot = function(x, breaks, col = "lightblue", border = NA, 
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

