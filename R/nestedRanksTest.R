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
# Thompson PG, Smouse PE, Scofield DG, Sork VL.  Accepted, Movement Ecology.
# What seeds tell us about birds: a multi-year analysis of acorn woodpecker
# foraging movements.
#
# Data available in Data Dryad at <http://doi.org/10.5061/dryad.64jk6>.


options(stringsAsFactors=FALSE)

# MWW.nested.test performs the nested ranks test
#    dat    : data frame containing group, treatment and value columns
#    n.iter : number of iterations for permutation (n.iter - 1 are random, n.iter-th is data)
#    title  : title to use when reporting test results
# The result of the test is printed, and if the return value is
# assigned, it is a data.frame containing the complete set of permuted values
# for all granaries with the test answers attached as attributes with the
# results of the test attached as attributes.  The data.frame can be passed to
# plot.MWW.nested.test() to plot the test results.
MWW.nested.test = function(dat, n.iter=10000, title=NULL) {
  dat.name = deparse(substitute(dat))
  if (is.null(title)) title = dat.name
  y = unique(sort(dat$treatment))
  if (length(y) != 2) 
    stop(dat.name, "requires exactly two levels for treatment")
  s = split(dat$group, dat$treatment)
  if (length(unique(sort(dat$group))) != length(intersect(s[[1]], s[[2]])))
    stop(dat.name, "must have values for all groups in both treatment levels")
  wt = MWW.weights(dat)
  Grps = as.character(sort(as.integer(unique(dat$group))))
  # fill in weight for each granary
  weights = wt$Rel_Wt
  names(weights) = paste0("g",wt$group)  # ensure that groups have proper R names
  # print(weights)
  # compute permutation for each granary individually
  p = list()
  for (g in Grps) {
    gdat = subset(dat, group == g)
    y1.dat = subset(gdat, treatment == y[1])
    y2.dat = subset(gdat, treatment == y[2])
    n1 = nrow(y1.dat)
    n2 = nrow(y2.dat)
    vals = c(y1.dat$value, y2.dat$value)
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
  attr(ans,"Z.weighted.obs") = Z.weighted[n.iter]
  attr(ans,"P.obs") = sum(Z.weighted >= Z.weighted[n.iter]) / n.iter
  attr(ans,"n.iter") = n.iter
  cat(title, 
      " Z.weighted.obs =", attr(ans, "Z.weighted.obs"), 
      " n.iter =", attr(ans, "n.iter"), 
      " P.obs =", attr(ans, "P.obs"),
      "\n")
  ####
  invisible(ans)
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
  y = unique(sort(dat$treatment))
  Grps = as.character(sort(as.integer(unique(dat$group))))
  w = data.frame()
  for (g in Grps) {
    gdat = subset(dat, group == g)
    y1.dat = subset(gdat, treatment == y[1])
    y2.dat = subset(gdat, treatment == y[2])
    n1 = nrow(y1.dat)
    n2 = nrow(y2.dat)
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
plot.MWW.nested.test = function(test.dat, title=NULL) {
  if (is.null(title)) title = deparse(substitute(test.dat))
  bks = if (max(test.dat$Z.weighted) > 1 || min(test.dat$Z.weighted) < -1)
          200
        else
          seq(-1,1,0.05)
  png(file=paste0(title, "_MWW_nested_test.png"), width=500, height=400)
  hist(test.dat$Z.weighted,
       breaks=bks,
       col="lightblue",
       border=NA,
       main=paste0(title, " weighted between-y Z-score\nacross all groups, P = ", 
                   attr(test.dat, "P.obs")),
       xlab="Weighted between-year Z-score",
       ylab=paste0("Frequency (out of ", attr(test.dat, "n.iter"), ")"))
  abline(v=attr(test.dat, "Z.weighted.obs"), col="red", lty=2, lwd=2)
  dev.off()
}

