# nestedRanksTest 1.0

* Initial release, available via `devtools::install_github()` but not yet 
  submitted to CRAN
* Provides S3 generic `nestedRanksTest()`, with function and default methods
* Provides `print` and `plot` methods for class `"htest_boot"`, which extends 
  class `"htest"` by including information on the generation and content of the 
  underlying bootstrapped null distribution used for generating p-values
* Provides a dataset `woodpecker_multiyear`
* TODO: modify Makefile to build both html and pdf vignettes
* TODO: manage References better in the vignette, include internal links
* TODO: complete the description of the interface
* TODO: complete test suite by adding error checks
