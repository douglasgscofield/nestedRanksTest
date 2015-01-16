\name{nestedRanksTest}
\alias{nestedRanksTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform Mann-Whitney-Wilcoxon test when data are nested in groups
}
\description{
Describe it.
}
\usage{
nestedRanksTest(dat, n.iter = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
\code{data.frame} containing group, treatment and value columns
}
  \item{n.iter}{
Number of iterations for permutation test (\code{n.iter - 1} are random, the \code{n.iter}-th is the data)
}
}
\details{
\code{nestedRanksTest()}  ...
Calling \code{nestedRanksTest()} ...
}
\value{
The result of the test is printed, and a \code{data.frame} containing the complete set of permuted values for all granaries with the test answers attached as attributes with the results of the test attached as attributes.  The \code{data.frame} can be passed to \code{plot.MWW.nested.test()} to plot the test results.
}
\references{
Thompson, P.G., P.E. Smouse, D.G. Scofield and V.L. Sork.  2014.  What seeds tell us about birds: a multi-year analysis of acorn woodpecker foraging movements.  \emph{Movement Ecology} 2, 12.  doi:10.1186/2051-3933-2-12, data at doi:10.5061/dryad.64jk6
}
\author{
Douglas G. Scofield
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(example_data)
example_data
attributes(example_data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ file }
\keyword{ manip }
\keyword{ attribute }
