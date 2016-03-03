## QMRA

[![Travis-CI Build Status](https://travis-ci.org/brechtdv/QMRA.svg?branch=master)](https://travis-ci.org/brechtdv/QMRA)

The `QMRA` package provides Maximum-likelihood and Bayesian parametric methods for exposure and dose-response assessment.

The easiest way to install the development version of the `QMRA` package is via the `devtools` package:

```r
devtools::install_github("brechtdv/QMRA")
library(QMRA)
```

The Bayesian functions in the `QMRA` package call on **`JAGS` (Just Another Gibbs Sampler)**, which therefore has to be available on the user's system. `JAGS` can be downloaded from http://mcmc-jags.sourceforge.net/.

**Mac users** may need to install `Xcode` from the App store to install the `QMRA` package. If the **`rjags` package** cannot be loaded, try reinstalling it and select install from source when prompted.
