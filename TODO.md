# For 1.0

* Can the perhaps build the README from a subset of the vignette?  Vignette needs some editing if so.
* Add figure to vignette explaining structure in the dataset.
* Complete the description of the interface (weighting, Z-score, bootstrap)
* How to actually run/view a vignette?  Add this to readme.
* More tests
* Can the bootstrapping be sped up?  The place to start might be to get rid of `for (i in seq_along(groups_levels)) {`.
* Can we add a trigger that does `make doc` before any `git commit`, to automatically keep the `man/*.Rd` files up to date?
* Keywords etc. to get it in a CRAN view for this kind of test
* Submit to CRAN
