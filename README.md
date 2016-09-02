# simplerspec 

Simplerspec aims to facilitate spectra and additional data handling and model development for spectroscopy applications such as FT-IR soil spectroscopy. Different helper functions are designed to create a 
data and modeling workflow. Data inputs and outputs are stored in `R` objects with specific data structures. The following steps are covered in the current beta version of the package:

1. Read spectral data from text files ((`.csv`); an implementation for reading OPUS binary files is planned)
2. Average spectra for replicate scans
3. Detect and remove outlier spectra based on robust PCA
4. Resample spectra to new wavenumber intervals
5. Perform pre-processing of spectra
6. Join chemial and spectral data sets
7. Perform PLS regression modeling
8. Predict chemical properties from a list of calibrated models and new soil spectra

# Installation

The newest version of the package is available on this GitHub repository. Note that the package is still under development. If you find bugs you are highly welcome to report your issues (write an [email](mailto:philipp.baumann@gmx.ch) or create an [issue](https://github.com/philipp-baumann/simplerspec/issues)) You can install `simplerspec` using the devtools package. Currently, there seems to be still an issue that `install_github()` does not automatically install all packages that are listed under "imports" (see [here](https://github.com/hadley/devtools/issues/1265)). In case you have error obtain error messages that packages can't be found, install the following packages:

```R
# List of packages to be installed
list_packages <- c("ggplot2", "plyr", "data.table", "reshape2",
  "mvoutlier", "hexView", "Rcpp", "hyperSpec", "prospectr",
  "dplyr", "caret")
# Install packages from CRAN
install.packages(list_packages, dependencies = TRUE)
```
Then run:

```R
# install.packages("devtools")
# Install the simplerspec package from the github repository
# (https://github.com/philipp-baumann/simplerspec)
devtools::install_github("philipp-baumann/simplerspec")
```



