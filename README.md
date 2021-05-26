# HonestDiD

Robust inference in differences-in-differences and event study designs using methods developed in [Rambachan and Roth (2019)](https://scholar.harvard.edu/jroth/publications/Roth_JMP_Honest_Parallel_Trends).

The [HonestDiD vignette](doc/HonestDiD_Example.pdf) provides a brief decsription of the package and 
an illustration to show users how to use the package to conduct sensitivity analysis on the parallel trends assumption 
in difference-in-differences and event study designs. The package [manual](doc/manual.pdf) provides detailed documentation of the package functions.

This software package is based upon work supported by the National Science
Foundation Graduate Research Fellowship under Grant DGE1745303 (Rambachan) and Grant DGE1144152 (Roth). 

## Installation

The package may be installed by using the function `install_github()` from the `remotes` package:

```
install.packages("remotes") # if remotes package not installed

# Turn off warning-error-conversion, because the tiniest warning stops installation
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")

# install from github
remotes::install_github("asheshrambachan/HonestDiD")
```
