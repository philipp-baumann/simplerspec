# simplerspec 

Simplerspec aims to facilitate spectra and additional data handling and model development for spectroscopy applications such as FT-IR soil spectroscopy. Different helper functions are designed to create a 
data and modeling workflow. Data inputs and outputs are stored in `R` objects with specific data structures. The following steps are covered in the current beta version of the package:

1. Read spectral data from text files (`.csv`; an implementation for reading OPUS binary files is planned)
2. Average spectra for replicate scans
3. Detect and remove outlier spectra based on robust PCA
4. Resample spectra to new wavenumber intervals
5. Perform pre-processing of spectra
6. Join chemial and spectral data sets
7. Perform calibration sampling and PLS regression modeling
8. Predict chemical properties from a list of calibrated models and new soil spectra

# Installation

The newest version of the package is available on this GitHub repository. Note that the package is still under development. If you find bugs you are highly welcome to report your issues (write an [email](mailto:philipp.baumann@gmx.ch) or create an [issue](https://github.com/philipp-baumann/simplerspec/issues)) You can install `simplerspec` using the devtools package. Currently, there seems to be still an issue that `install_github()` does not automatically install all packages that are listed under "imports" (see [here](https://github.com/hadley/devtools/issues/1265)). In case you obtain error messages that packages can't be found, install the following packages:

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

# Key concepts and data analysis workflow

The functions are built to work in a pipeline and cover commonly used procedures for spectral model development. Many R packages are available to do tasks in spectral modeling such as pre-processing of spectral data. The motivation to create this package was:

1. Avoid repetition of code in model developement (common source of errors)
2. Provide a reproducible data analysis workflow for FT-IR spectroscopy
3. R packges are an ideal way to organize and share R code
4. Make soil FT-IR spectroscopy modeling accessible to people that have basic R knowledge

This package builds mainly upon functions from the following packages:

* `prospectr `: 

Consistent and reproducible data and metadata management is a important prerequisite for spectral model development. Therefore, different outputs should be stored as R objects in a consistent way using R data structures. Simplerspec functions use lists as R data structures because they allow to store complex, hierarchical objects in a flexible way. Lists can e.g. contain other lists, vectors, data.frames, or matrices.

# Example workflow

In a fist step, the spectra (one file per spectrum and replicate scan) are read from
the text (`.txt`) files. Currently, an export macro within the Bruker OPUS software
is used to convert OPUS binary files to spectra in the form of a text file.
The argument `path` specifies the the folder where all spectral files to be loaded
into R are located. The files contain two columns that are comma-separated. The first
is the wavenumber and the second is the absorbance value.

```R
# Read spectra in text format (Alpha spectrometer) -----------------------------
# Currently 
soilspec_in <- read_spectra(
  path = "data/spectra/alpha_txt"
)
```

Pipes can make R code more readable and fit to the stepwise data processing
in the context of developing spectral models. The pipe operator (%>%, called "then") is a new operator in R that was introduced
with the magrittr package. It facilitates readability of code
and avoids to type intermediate objects. The basic behaviour of
the pipe operator is
that the object on the left hand side is passed as the first argument
to the function on the right hand side.
More details can be found [here](https://github.com/smbache/magrittr).

The model development process can be quickly coded as the example below illustrates:

```R
################################################################################
## Part 1: Read and pre-process spectra, Read chemical data, and join
## spectral and chemical data sets
################################################################################

# Average, remove outlier, resample, then pre-process spectra ------------------
soilspec <- soilspec_in %>% average_spectra() %>%
  remove_outliers(remove = FALSE) %>%
  resample_spectra(wn_lower = 510, wn_upper = 3988, wn_interval = 2) %>%
  do_pretreatment(select = "MIR0")
  
# Read chemical data from csv (comma-separated values) file --------------------
soilchem <- read.csv(
  file = "out/data/soilchem_yamsys.csv"
)

# Join chemical and spectra data -----------------------------------------------
spec_chem <- join_chem_spec(dat_chem = soilchem, dat_spec = soilspec)

################################################################################
## Part 2: Run PLS regression models for different soil variables
################################################################################

# Example Partial Least Squares (PLS) regression model for total Carbon (C)
pls_C <- pls_ken_stone(
  spec_chem = spec_chem[!is.na(spec_chem$C), ],
  ratio_val = 1/3,
  variable = C,
  validation = TRUE
)
```

# Details on functions

## `read.spectra()` function


