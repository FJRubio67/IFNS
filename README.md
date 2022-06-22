# IFNS

## Individual frailty excess hazard models in cancer epidemiology

This repository contains a real data application of the clustered survival models proposed in the paper:

> Rubio, F.J., Putter, H. and Belot, A. (2022). Individual frailty excess hazard models in cancer epidemiology. Under review.

The models are fitted using the R package `IFNS`. To install the `IFNS` R package use:

```
library(devtools)
install_github("FJRubio67/IFNS")

library(IFNS)
```

An example of the use of this package using real data can be found at:

[Simulacrum data - lung cancer: Individual Frailty Model](https://rpubs.com/FJRubio/IFNSSimulacrum)

See also: [GHSurv](https://github.com/FJRubio67/GHSurv) and [LBANS](https://github.com/FJRubio67/LBANS)
