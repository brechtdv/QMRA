## QMRA

[![Travis-CI Build Status](https://travis-ci.org/brechtdv/QMRA.svg?branch=master)](https://travis-ci.org/brechtdv/QMRA)

The QMRA package provides Maximum-likelihood and Bayesian parametric methods for exposure and dose-response assessment.

The easiest way to install the development version of the `QMRA` package is via the `devtools` package:

```r
devtools::install_github("brechtdv/QMRA")
library(QMRA)
```

IMPORTANT: the Bayesian functions in the QMRA package call on JAGS (Just Another Gibbs Sampler), which therefore has to be available on the user's system. JAGS can be downloaded from http://mcmc-jags.sourceforge.net/.
