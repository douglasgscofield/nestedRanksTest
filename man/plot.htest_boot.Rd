% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/htest_boot-class.R
\name{plot.htest_boot}
\alias{plot.htest_boot}
\title{Diagnostic plot of result held in htest_boot object}
\usage{
\method{plot}{htest_boot}(x, breaks, col = "lightblue", border = NA,
  digits = getOption("digits"), main = paste(sep = "", x$method, ", ",
  x$data.name, "\\n", names(x$statistic), " = ", format(signif(x$statistic,
  max(1L, digits - 2L))), ", P = ", format.pval(x$p.value, digits = max(1L,
  digits - 3L))), xlab = "Distribution of Z-scores", ylab = paste(sep = "",
  "Frequency (out of ", x$n.iter, ")"), p.col = "red", p.lty = 2,
  p.lwd = 2, ...)
}
\arguments{
\item{x}{Value of class \code{'htest_boot'}}

\item{breaks}{The number of breaks to use when plotting the distribution,
the default is calculated from \code{n.iter} of the call to
\code{nestedRanksTest}}

\item{col}{Fill color for histogram bars, passed to \code{hist}.}

\item{border}{Border color for histogram bars, passed to \code{hist}}

\item{digits}{Number of digits to use for statistic and p-value, the
default is taken from the \code{"digits"} option}

\item{main}{Main title, passed to \code{hist}}

\item{xlab}{X-axis label, passed to \code{hist}}

\item{ylab}{Y-axis label, passed to \code{hist}}

\item{p.col}{Observed value line colour, passed to \code{abline}}

\item{p.lty}{Observed value line type, passed to \code{abline}}

\item{p.lwd}{Observed value line width, passed to \code{abline}}

\item{\dots}{Additional arguments passed to \code{hist} and
\code{abline} for plotting}
}
\value{
None
}
\description{
\code{plot.htest_boot} creates a diagnostic plot of a return value
held in a list of class \code{'htest_boot'}, as returned by e.g., 
\code{\link{nestedRanksTest}}.  The plot contains a histogram of the
null distribution generated by' bootstrapping plotted with
\code{\link{hist}}, and a verticle line indicating the observed value
plotted with \code{\link{abline}}.
}
\details{
If there is no null distribution included in the class, e.g., because 
the' options \code{lightweight = TRUE} or \code{n.iter = 1} were given to
\code{nestedRanksTest}, this function produces an error.
}
\examples{
require(graphics)

data(woodpecker_multiyear)

## n.iter set to 1000 to shorten completion time
res.a <- nestedRanksTest(Distance ~ Year | Granary, n.iter = 1000,
                         data = woodpecker_multiyear,
                         subset = Species == "agrifolia")
res.l <- nestedRanksTest(Distance ~ Year | Granary, n.iter = 1000,
                         data = woodpecker_multiyear,
                         subset = Species == "lobata")

opa = par(mfrow = c(2, 1))
## Defaults
plot(res.l)
## Modify colours, line type and main title
plot(res.a, main = "Quercus agrifolia", col = "lightgreen",
     p.col = "brown4", p.lty = 1)
par(opa)

}
\seealso{
\code{\link{print.htest_boot}} for printing test results,
  \code{\link{hist}} and \code{\link{abline}} for plotting options,
  and \code{\link{nestedRanksTest}} for one test that
  returns an object of class \code{'htest_boot'}
}
\keyword{hplot}

