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

# TODO: find R function to make strings into proper R names
# TODO: does a two-sided test make sense here?  I don't think so...
# TODO: can I extend this with a nestedRanksTest class allowing for printing and plotting
# TODO: end print.nestedRanksTest() with invisible(x)
# http://tolstoy.newcastle.edu.au/R/devel/05/03/0094.html
# https://github.com/Rapporter/rapportools/blob/master/R/htest.R
# http://www1.maths.lth.se/matstat/bioinformatics/software/R/library/ctest/html/print.pairwise.htest.html

nestedRanksTest.formula = function(formula, data, groups = NULL, ...)
{
  # code largely copied from stats:::t.test
  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]), 
                                                                  "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  # TODO: trace these called, what happens to m
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")  # make sure this is correct
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  tmt <- factor(mf[[-response]])
  if (nlevels(tmt) != 2L)
    stop("treatment factor must have exactly 2 levels")
  ## handle grouping: groups = NULL if unspecified
  # TODO: deal with groups, this is preliminary, from latticeParseFormula()
  # If groups is a variable, then its value is used for groups, otherwise the
  # code below turns the formula into the group variable
  groupVar = deparse(substitute(groups))  # start by assuming groups is a vector
  if (inherits(groups, "formula")) {
    groupVar <- as.character(groups)[2]
    groups <- eval(parse(text = groupVar), data, environment(groups))
  }
  # TODO: parse out grouping if exists in formula
  if (length(formula[[3]]) == 3) {  # grouping specified in formula
    if (formula[[3]][[1]] = as.name("|")) {  # indeed
      if (! is.null(groups))
        stop("cannot specify grouping variable in both formula and groups argument")
      tmt <- as.character(formula[[3]][[2]])
      groupVar <- as.character(formula[[3]][[3]])
      groups <- eval(parse(text = groupVar), data, environment(formula))
    } else stop("invalid group specification in formula")
  }
  if (is.null(groups))
    stop("group structure must be specified in the formula or groups argument")
  DNAME <- paste(DNAME, groupVar, sep = " grouped by ")
  # TODO: must modify DNAME to incorporate groups
  # TODO: consider renaming these...
  DATA <- setNames(split(mf[[response]], tmt), c("x", "y"))
  # TODO: nestedRanksTest doesn't yet exist
  y <- do.call("nestedRanksTest", c(DATA, list(...)))
  y$data.name <- DNAME
  # TODO: are we returning group means?  is this appropriate?  wilcox.test lacks it
  if (length(y$estimate) == 2L)
    names(y$estimate) <- paste("mean in treatment", levels(tmt))
  y
}

MWW.nested.test = function(dat, n.iter=10000)
{
  DNAME = deparse(substitute(dat))
  METHOD = "Nested Ranks Test"
  tmt_levels = unique(sort(dat$treatment))
  if (length(tmt_levels) != 2) 
    stop(DNAME, "requires exactly two levels for treatment")
  # implement a function version
  nr = nrow(dat)
  dat = dat[apply(dat[,1:3], 1, function(x) all(!is.na(x))), 1:3]
  BAD.OBS = nr - nrow(dat)
  vals = dat[, 1]
  tmt = dat[, 2]
  grp = dat[, 3]
  s = split(grp, tmt)
  if (length(unique(sort(grp))) != length(intersect(s[[1]], s[[2]])))
    stop(DNAME, "must have values for all groups in both treatment levels")
  wt = MWW.weights(dat)
  Grps = unique(sort(as.character(grp)))
  # fill in weight for each granary
  weights = wt$Rel_Wt
  # TODO: find R function to make strings into proper R names
  names(weights) = paste0("g", wt$group)  # ensure that groups have proper R names
  # print(weights)
  # compute permutation for each group individually
  p = list()
  for (g in Grps) {
    gdat = subset(dat, group == g)
    tmt1.dat = subset(gdat, treatment == tmt_levels[1])
    tmt2.dat = subset(gdat, treatment == tmt_levels[2])
    n1 = nrow(tmt1.dat)
    n2 = nrow(tmt2.dat)
    vals = c(tmt1.dat$value, tmt2.dat$value)
    this.z = MWW(vals, n1, n2)
    Z = numeric(n.iter)
    if (n.iter > 1) {
      for (i in 1:(n.iter - 1)) {
        d = sample(vals)
        Z[i] = MWW(d, n1, n2)
      }
    }
    Z[n.iter] = this.z
    p[[paste0("g",g)]] = Z
  }
  ans = as.data.frame(p)
  if (! all(names(ans) == names(weights)))
    stop("weights out of order")
  Z.weighted = apply(ans, 1, function(.x) sum(.x * weights))
  ans$Z.weighted = Z.weighted
  attr(ans,"weights") = weights
  #attr(ans,"Z.weighted.obs") = Z.weighted[n.iter]
  STATISTIC = Z.weighted[n.iter]
  names(STATISTIC) = "Z.weighted.obs"
  #attr(ans,"P.obs") = sum(Z.weighted >= Z.weighted[n.iter]) / n.iter
  PVAL = sum(Z.weighted >= Z.weighted[n.iter]) / n.iter
  attr(ans,"n.iter") = n.iter
  cat(" Z.weighted.obs =", attr(ans, "Z.weighted.obs"), 
      " n.iter =", attr(ans, "n.iter"), 
      " P.obs =", attr(ans, "P.obs"),
      "\n")
  RVAL = list(statistic = STATISTIC,
              p.value = PVAL,
              alternative = "Z.weighted.obs lies above the values derived from bootstrapping",
              method = METHOD,
              data.name = DNAME,
              bad.obs = BAD.OBS)
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
  tmt_levels = unique(sort(dat$treatment))
  Grps = as.character(sort(as.integer(unique(dat$group))))
  w = data.frame()
  for (g in Grps) {
    gdat = subset(dat, group == g)
    tmt1.dat = subset(gdat, treatment == tmt_levels[1])
    tmt2.dat = subset(gdat, treatment == tmt_levels[2])
    n1 = nrow(tmt1.dat)
    n2 = nrow(tmt2.dat)
    n1.n2 = n1 * n2
    w = rbind(w, list(group=g,
                      n1=n1,
                      n2=n2,
                      n1.n2=n1.n2))
  }
  w = transform(w, Rel_Wt=n1.n2/sum(n1.n2))
  rownames(w) = paste0("g", w$group)
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

