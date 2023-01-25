# IFNS R package

## Individual frailty excess hazard models in cancer epidemiology

This repository contains a real data application of the individual frailty excess hazard models proposed in the paper:

> Rubio, F.J., Putter, H. and Belot, A. (2022). Individual frailty excess hazard models in cancer epidemiology. Statistics in Medicine, in press. https://doi.org/10.1002/sim.9657 . [ [Preprint] ](https://drive.google.com/file/d/16Jc6T4EOgIAoSJa0IJM-kN8hVAV9cZDG/view)

The models are fitted using the R package `IFNS`. To install the `IFNS` R package use:

```
library(devtools)
install_github("FJRubio67/IFNS")

library(IFNS)
```

An example of the use of this package using simulated data can be found at:

[Simulacrum data - lung cancer: Individual Frailty Model](https://rpubs.com/FJRubio/IFNSSimulacrum)

The data set used in this example is available at [LBANS](https://github.com/FJRubio67/LBANS).

See also: [GHSurv](https://github.com/FJRubio67/GHSurv) and [LBANS](https://github.com/FJRubio67/LBANS)
