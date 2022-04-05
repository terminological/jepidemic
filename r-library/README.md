# jepidemic

R and Java based tools for epidemic modelling

```R

J = jepidemic::JavaApi$new()
estim = J$CoriEstimator$new(r0Mean = 5,r0SD = 5,maxWindow = 14)
library(EpiEstim)
data(Flu2009)
estim$withInfectivityProfile(vector = Flu2009$si_distr)

```

