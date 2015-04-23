# For next version

* Add class description to `R/htest_boot-class.R`
* The forthcoming R includes `getOption("digits")`, add this check to `print.htest_boot`, see R development code at <https://github.com/wch/r-source/commit/7e90613c381a011a2f1a625bab5dfb01b63e83ec#diff-2fb78b80e16e1dbbb7e60b6721702ea1>
* Is there more to be done generalising htest_boot?  Should it be a separate package?
* Faster alternative to base R `rank` as that is the bootstrapping bottleneck
* Is travis-ci.org a solution for continual integration?
* Include link to CRAN package passing stats, along with other tags, as is visible in e.g., https://github.com/HenrikBengtsson/matrixStats
