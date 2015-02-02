# nestedRanksTest 1.0

* Initial release, available via `devtools::install_github()` but not yet 
  submitted to CRAN
* Provides S3 generic `nestedRanksTest()`, with function and default methods
* Provides `print` and `plot` methods for class `"htest_boot"`, which extends 
  class `"htest"` by including information on the generation and content of the 
  underlying bootstrapped null distribution used for generating p-values
* Provides dataset `woodpecker_multiyear`, distances acorn woodpeckers moved
  acorns of two different oak species from source tree to storage granary during
  two years per oak species
