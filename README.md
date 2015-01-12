nestedRanksTest
===============

Calculate a Mann-Whitney-Wilcoxon (MWW) test using nested ranks.  Data are structured into several groups and each group has received two treatment levels. The rankings are compared between treatment levels, taking group structure and size into account.  A null distribution is constructed by permuting ranks, assuming no influence of treatment.  When there is just one group, this test is identical to a standard MWW test. 

The function `MWW.nested.test(dat, n.iter=10000, title=NULL)` takes three arguments:

* `dat` is a data frame containing three columns.
	1. `$group`, an alphanumeric group identifier
	2. `$treatment`, the treatment level, there must be two treatment levels
	3. `$value`, a numerical value to be ranked between treatments.  All groups must have values in both treatment levels.
* `n.iter` is the number of permutations to use to create the null distribution.  The number of simulated distributions is `n.iter - 1`, with the observed data providing the `n.iter`-th value
* `title` is a title to use when reporting test results, if not provided it is taken from the name of the input data frame

The test results are printed to the console, and are also invisibly returned in a data frame.  The data frame contains the weighted Z-scores across all groups for each iteration, with the observed Z-scores as the last row.  The data frame also contains the following attached attributes:

* `weights` which are the group-size weights used when combining Z-scores
* `Z.weighted.obs` is the Z-score for the observed data
* `P.obs` is the *P*-value as calculated from the empirical distribution
* `n.iter` as provided to the function

The returned data frame can be passed to the utility function `plot.MWW.ntested.test(test.dat, title=NULL)` to plot the simulated distribution and the observed value.


Publication
-----------

We described this test in Thompson _et al._ 2014 [_Movement Ecology_ 2:12](http://dx.doi.org/10.1186/2051-3933-2-12).  We used the MWW nested ranks test to compare dispersal distances between years for a collection of acorn woodpecker granaries, but the test is more generally useful for any MWW-type comparison between two treatment levels when the data are divided into two or more discrete groups.  Each group must be present in both treatment levels, but the number of measurements within each group need not be the same for each level.

The published data are available from Data Dryad at <http://doi.org/10.5061/dryad.64jk6>.


Reference
---------

Thompson PG, Smouse PE, Scofield DG, Sork VL.  2014.  What seeds tell us about birds: a multi-year analysis of acorn woodpecker foraging movements.  [_Movement Ecology_ 2:12](http://dx.doi.org/10.1186/2051-3933-2-12), [data](http://doi.org/10.5061/dryad.64jk6).

