# Nested Ranks Test
#
# Calculate a Mann-Whitney-Wilcoxon (MWW) test using nested ranks.  Data are
# structured into several groups and each group has received two treatment
# levels. The rankings are compared between treatment levels, taking group
# structure and size into account when permuting ranks to construct a null
# distribution which assumes no influence of treatment.  When there is just one
# group, this test is identical to a standard MWW test.
#
# Reference:
#
# Thompson PG, Smouse PE, Scofield DG, Sork VL. 2014.  What seeds tell us 
# about birds: a multi-year analysis of acorn woodpecker foraging movements.
# Movement Ecology 2:12. doi:10.1186/2051- 3933-2-12, data at 
# doi:10.5061/dryad.64jk6


options(stringsAsFactors=FALSE)

# MWW.nested.test performs the nested ranks test
#    dat    : data frame containing group, treatment and value columns
#    dat    : column 1: value
#           : column 2: treatment
#           : column 3: group
#    n.iter : number of iterations for permutation (n.iter - 1 are random, n.iter-th is data)
# The result of the test is printed, and if the return value is
# assigned, it is a data.frame containing the complete set of permuted values
# for all granaries with the test answers attached as attributes with the
# results of the test attached as attributes.  The data.frame can be passed to
# plot.MWW.nested.test() to plot the test results.

# TODO: can I extend this with a nestedRanksTest class allowing for printing and plotting
# TODO: end print.nestedRanksTest() with invisible(x)
# http://tolstoy.newcastle.edu.au/R/devel/05/03/0094.html
# https://github.com/Rapporter/rapportools/blob/master/R/htest.R
# http://www1.maths.lth.se/matstat/bioinformatics/software/R/library/ctest/html/print.pairwise.htest.html

nestedRanksTest.formula = function(formula, data, groups = NULL, subset, ...)
{
  # initial version largely copied from stats:::t.test
  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]), 
                                                                  "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  # TODO: trace these called, what happens to m
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))  # if data passed in as a matrix
    m$data <- as.data.frame(data)
  # Now, expand the formula, see if we have groups there via '|' or via groups=
  if (is.null(m$groups)) {
    if (length(formula[[3]]) == 3L && formula[[3]][[1]] == as.name("|")) { # grouping in formula
      # modify formula, add groups=
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
  m[[1L]] <- quote(stats::model.frame)  # this could do a lot, including na.action and subset
  m$... <- NULL
  mf <- eval(m, parent.frame())  # for (y ~ x, groups=G) mf now y , x , "(groups)"
  if (! nrow(mf))
    stop("no data")
  DNAME <- paste(names(mf)[1:2], collapse = " by ")  # only response and treatment
  DNAME <- paste(DNAME, GROUP.NAME, sep = " grouped by ")
  response <- attr(attr(mf, "terms"), "response")
  TREATMENT.NAME <- attr(attr(mf, "terms"), "term.labels")
  tmt <- factor(mf[[TREATMENT.NAME]])
  if (nlevels(tmt) != 2L)
    stop("treatment factor must have exactly 2 levels")
  # TODO: consider renaming these...
  DATA <- data.frame(y=mf[[response]], x=tmt, groups=mf[["(groups)"]])
  # TODO: nestedRanksTest doesn't yet exist
  y <- do.call("nestedRanksTest.default", 
               c(list(y=mf[[response]], x=tmt, groups=mf[["(groups)"]]), list(...)))
  #y <- do.call("MWW.nested.test", c(list(dat = DATA), list(...)))
  y$data.name <- DNAME
  y
}


# must suffix with '.default' when a proper S3 class
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
  wt = MWW.weights(dat)
  Grps = unique(sort(as.character(grp)))
  nm.Grps = setNames(make.names(Grps, unique = TRUE), Grps)
  # fill in weight for each granary
  weights = setNames(wt$Rel_Wt, rownames(wt))
  # print(weights)
  # compute permutation for each group individually
  # TODO: speed this up...
  p = list()
  for (g in Grps) {
    tmt1.dat = dat[grp == g & tmt == tmt_levels[1], ]
    tmt2.dat = dat[grp == g & tmt == tmt_levels[2], ]
    n1 = nrow(tmt1.dat)
    n2 = nrow(tmt2.dat)
    this.vals = c(tmt1.dat[, 1], tmt2.dat[, 1])
    this.z = MWW(this.vals, n1, n2)
    Z = numeric(n.iter)
    if (n.iter > 1) {
      for (i in 1:(n.iter - 1)) {
        d = sample(this.vals)
        Z[i] = MWW(d, n1, n2)
      }
    }
    Z[n.iter] = this.z
    p[[ nm.Grps[g] ]] = Z
  }
  p = as.data.frame(p)
  if (! all(names(p) == names(weights)))
    stop("weights out of order")
  NULL.DISTRIBUTION = apply(p, 1, function(x) sum(x * weights))
  N.ITER = n.iter
  WEIGHTS = weights
  STATISTIC = setNames(NULL.DISTRIBUTION[n.iter], STATISTIC.NAME)
  PVAL = sum(NULL.DISTRIBUTION >= STATISTIC) / N.ITER
  RVAL = list(statistic = STATISTIC,
              p.value = PVAL,
              alternative = paste(STATISTIC.NAME, "lies above bootstrapped values"),
              method = METHOD,
              data.name = DNAME,
              bad.obs = BAD.OBS,
              weights = WEIGHTS,
              n.iter = N.ITER,
              null.distribution = NULL.DISTRIBUTION)
  class(RVAL) = "htest"
  return(RVAL)
}

# MWW calculates the Mann-Whitney-Wilcoxon Z-score
#    x  : values
#    n1 : the first n1 of x belong to the first level
#    n2 : the final n2 of x belong to the second level
# The return value is the calculated Z-score
MWW = function(x, n1, n2) {
  r = rank(x)
  r1 = r[1:n1]
  r2 = r[(n1+1):(n1+n2)]
  R1 = sum(r1)
  R2 = sum(r2)
  n1.n2 = n1 * n2
  U1 = (R1 - (((n1+1)*n1)/2))
  U2 = (R2 - (((n2+1)*n2)/2))
  Z.12 = (U2 - U1) / n1.n2
  ####
  Z.12
}


# MWW.weights calculates sample-size weights
MWW.weights = function(dat) {
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
  ####
  w
}

# plot.MWW.nested.test creates a utility plot of the test result
plot.MWW.nested.test = function(test.dat) {
  title = deparse(substitute(test.dat))
  bks = if (max(test.dat$Z.weighted) > 1 || min(test.dat$Z.weighted) < -1)
          200
        else
          seq(-1,1,0.05)
  hist(test.dat$Z.weighted,
       breaks=bks,
       col="lightblue",
       border=NA,
       main=paste0(title, " weighted between-tmt_levels Z-score\nacross all groups, P = ", 
                   attr(test.dat, "P.obs")),
       xlab="Weighted between-year Z-score",
       ylab=paste0("Frequency (out of ", attr(test.dat, "n.iter"), ")"))
  abline(v=attr(test.dat, "Z.weighted.obs"), col="red", lty=2, lwd=2)
}

