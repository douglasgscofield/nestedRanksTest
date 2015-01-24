nestedRanksTest
===============

Tne `nestedRanksTest` package provides functions for performing a
Mann-Whitney-Wilcoxon-type nonparametric test for a difference between
treatment levels using nested ranks, and also provides functions for displaying
results of the test.  The nested ranks test may be used when observations are
structured into several groups and each group has received both treatment
levels.  The p-value is determined via bootstrapping.

The development version is hosted here and can be installed via:

```R
> install.packages("devtools")
> devtools::install_github("douglasgscofield/nestedRanksTest")
```

Help is available via `help("nestedRanksTest")`.

* * *
These statistical tools were developed in collaboration with Peter E.
Smouse (Rutgers University) and Victoria L. Sork (UCLA) and were funded by
U.S. National Science Foundation awards NSF-DEB-0514956 and
NSF-DEB-0516529.
* * *


To Use It
---------

The principle function is `nestedRanksTest()`, with two interfaces.  The
formula interface allows specification of quantitative measures, treatments,
and group membership using formula syntax.

```R
data(woodpecker_multiyear)
result <- nestedRanksTest(Distance ~ Year | Granary, data = woodpecker_multiyear,
                          subset = Species == "agrifolia")
print(result)
```

~~~~

	Nested Ranks Test

data:  Distance by Year grouped by Granary
Z = 0.277, p-value = 1e-04
alternative hypothesis: Z lies above bootstrapped null values
null values:
      0%       1%       5%      10%      25%      50%      75%      90%      95%      99%     100% 
-0.29492 -0.15583 -0.11059 -0.08705 -0.04554 -0.00124  0.04430  0.08488  0.10936  0.15335  0.27695 

bootstrap iterations = 10000 
group weights:
        10         31        140        151        152        938        942 
0.05204461 0.04646840 0.02478315 0.14560099 0.30359356 0.29120198 0.13630731 
~~~~

```R
plot(result)
```
![nestedRanksTest plot example](example_plot.png "nestedRanksTest plot example")


The default interface uses arguments for specifying the variables.

```R
result <- nestedRanksTest(Year, Distance, Granary, ...)
```

Details
-------

The statistic for the nested ranks test is a Z-score calculated by comparing
ranks between treatment levels, with contributions of each group to the final
Z-score weighted by group size.  The p-value is determined by comparing the
observed Z-score against a distribution of Z-scores calculated by bootstrapping
ranks assuming no influence of treatment while respecting group sizes. When
there is just one group, this test is essentially identical to a standard
Mann-Whitney-Wilcoxon test.

We first described this test in [Thompson _et al._ 2014 _Movement Ecology_
2:12](http://www.movementecologyjournal.com/content/2/1/12).   We examined
year-to-year differences in acorn movement by acorn woodpeckers in an oak
savannah in central California.  Acorn woodpeckers are highly social, and each
social group maintains its own granary in which it stores acorns.  The
woodpeckers are highly territorial and occupy relatively stable territories,
each containing a number of mature oak trees at varying distances from the
granary.  Because each granary has its own neighbourhood of oak trees, it would
not be a very precise test to look for overall differences in acorn movement
within the entire savannah, on the other hand looking for differences within
each granary individually would have required greater than practical sample
sizes to have sufficient statistical power.  We chose instead to test whether
the distance ranks differ between years within each granary, and combine
results across granaries for an aggregate between-year test.

Note: The generation of a null distribution can take some time; for example
completing the default `n.iter = 10000` may require a few seconds.

Return value
------------

`nestedRanksTest()` returns a list of class `"htest_boot"` containing the following components:

Component |  Contents
----------|----------
`statistic` | the value of the observed Z-score.
`p.value` | the p-value for the test.
`weights` | the weights for groups, calculated by `nestedRanksTest_weights`.
`n.iter` | the number of bootstrap iterations used for generating the null distribution.
`null.values` | quantiles of the null distribution used for calculating the p-value.
`null.distribution` | (missing when `lightweight = TRUE`) a vector of length `n.iter` containing all values of the null distribution of Z-scores, with `statistic` as the last value.
`alternative` | a character string describing the alternative hypothesis.
`method` | a character string indicating the nested ranks test performed.
`data.name` | a character string giving the name(s) of the data..
`bad.obs` | the number of observations in the data excluded because of `NA` values.

Dataset
-------

The package also includes a dataset, `woodpecker_multiyear`, which includes the
data on woodpecker acorn movement underlying Figure 2 in [Thompson _et al._
2014 _Movement Ecology_ 2:12](http://www.movementecologyjournal.com/content/2/1/12).

References
----------

Thompson PG, Smouse PE, Scofield DG, Sork VL.  2014.  What seeds tell us about
birds: a multi-year analysis of acorn woodpecker foraging movements.  _Movement
Ecology_ **2**: 12, available as Open Access at
<http://www.movementecologyjournal.com/content/2/1/12>

<https://github.com/douglasgscofield/nestedRanksTest>

