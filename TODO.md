# For 1.0

* use Authors@R in DESCRIPTION
* Document version dependencies: `knitr (>= 0.5)` (for rmarkdown in the vignette), anything else in the vignette?
* Removed `paste0()`, which was introduced in R 2.15.0.  Are there other unknown dependencies on later R functions?
* Attempt to build for Windows and linux, is travis-ci.org a solution?
* Add figure to vignette explaining structure in the dataset.
* Speed up bootstrapping by writing a C function using the user-visible `R_orderVector` described at http://cran.r-project.org/doc/manuals/r-release/R-exts.html#Utility-functions
* Keywords etc. to get it in a CRAN view for this kind of test
* Submit to CRAN
