# For 1.0

* Document version dependencies: `knitr (>= 0.5)` (for rmarkdown in the vignette), anything else in the vignette?
* Removed `paste0()`, which was introduced in R 2.15.0.  Are there other unknown dependencies on later R functions?
* Attempt to build for Windows and linux, is travis-ci.org a solution?
* Add figure to vignette explaining structure in the dataset.
* Can the bootstrapping be sped up?  The place to start might be to get rid of `for (i in seq_along(groups_levels)) {`.
* Keywords etc. to get it in a CRAN view for this kind of test
* Submit to CRAN
