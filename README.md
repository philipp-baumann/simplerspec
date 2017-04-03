# simplerspec 

The simplerspec package aims to facilitate spectra and additional data handling and model development for spectroscopy applications such as FT-IR soil spectroscopy. Different helper functions are designed to create a 
data and modeling workflow. Data inputs and outputs are stored in `R` objects with specific data structures. The following steps are covered in the current version of the package:

1. Read spectral data from text files (`.csv`) and Bruker OPUS binary files
2. Average spectra for replicate scans
3. Detect and remove outlier spectra based on robust PCA
4. Resample spectra to new wavenumber intervals
5. Perform pre-processing of spectra
6. Join chemial and spectral data sets
7. Perform calibration sampling and Partial Least Squares (PLS) regression modeling
8. Predict chemical properties from a list of calibrated models and new soil spectra

# Installation

The newest version of the package is available on this GitHub repository. Note that the package is still under development. If you find bugs you are highly welcome to report your issues (write me an [email](mailto:philipp.baumann@gmx.ch) or create an [issue](https://github.com/philipp-baumann/simplerspec/issues)). You can install `simplerspec` using the devtools package. Currently, there seems to be still an issue that `install_github()` does not automatically install all packages that are listed under "imports" (see [here](https://github.com/hadley/devtools/issues/1265)). In case you obtain error messages that packages can't be found, install the following packages:

```R
# List of packages to be installed
list_packages <- c("ggplot2", "plyr", "data.table", "reshape2",
  "mvoutlier", "hexView", "Rcpp", "hyperSpec", "prospectr",
  "dplyr", "caret", "tidyverse")
# Install packages from CRAN
install.packages(list_packages, dependencies = TRUE)
```
Then run:

```R
# Uncomment and run the below line if you have not yet installed
# the devtools package
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
5. Provide a package interface that keeps data with various structures for spectral modeling related in R objects

This package builds mainly upon functions from the following R packages:

* `prospectr `: Various utilities for pre-processing and sample selection based on spectroscopic data. An introduction to the package with examples can be found [here](http://antoinestevens.github.io/prospectr/).
* `plyr` and `dplyr `: Fast data manipulation tools with an unified interface. See [here](https://github.com/hadley/dplyr) for details.
* `ggplot2 `: Alternative plotting system for R, based on the grammar of graphics. See [here](http://ggplot2.org/).
* `caret `: Classification and regression training. A set of functions that attempt to streamline the process for creating predictive models. See [here](http://topepo.github.io/caret/index.html) for details.

Consistent and reproducible data and metadata management is a important prerequisite for spectral model development. Therefore, different outputs should be stored as R objects in a consistent way using R data structures. Simplerspec functions uses tibble data frames as principal data structures because they allow to store lists within the well-known data frame structures. Lists are flexible and can e.g. contain other lists, vectors, data.frames, or matrices.

# Example workflow

In a fist step, the spectra (one file per spectrum and replicate scan) are read from the text (`.txt`) files. Currently, an export macro within the Bruker OPUS
software is used to convert OPUS binary files to spectra in the form of a text 
file. The argument `path` specifies the the folder where all spectral files to 
be loaded into R are located. The files contain two columns that are 
comma-separated. The first is the wavenumber and the second is the absorbance 
value.

```R
# Load simplerspec package for spectral model development wrapper functions
require(simplerspec)
# Load tidyverse package: loads packages frequently used for data manipulation,
# data tidying, import, and plotting
require(tidyverse)

################################################################################
## Part 1: Read and pre-process spectra, Read chemical data, and join
## spectral and chemical data sets
################################################################################

## Read spectra in list ========================================================

# List of OPUS binary spectra files
lf <- list.files("data/spectra/soilspec_eth_bin/", full.names = T)

# Read spectra from files into R list
spc_list <- read_opus(
  fnames = lf_eth,
  in_format = c("binary"),
  out_format = "list"
)
```

Pipes can make R code more readable and fit to the stepwise data processing
in the context of developing spectral models. The pipe operator (`%>%`, called "then") is a new operator in R that was introduced
with the magrittr package. It facilitates readability of code
and avoids to type intermediate objects. The basic behaviour of
the pipe operator is
that the object on the left hand side is passed as the first argument
to the function on the right hand side. When loading the tidyverse package, the
pipe operator is attached to the current R session
More details can be found [here](https://github.com/smbache/magrittr).

The model development process can be quickly coded as the example below illustrates:

```R
## Spectral data processing pipe ===============================================

soilspec_tbl <- spc_list %>%
  # Gather list into tibble data frame
  gather_spc() %>% 
  # Resample spectra to new wavenumber interval
  resample_spc(wn_lower = 500, wn_upper = 3996, wn_interval = 2) %>%
  # Average replicate scans per sample
  average_spc() %>%
  # Preprocess spectra using Savitzky-Golay first derivative with a window size
  # of 21 points
  preprocess_spc(select = "sg_1_w21")
  
## Read chemical reference data and join with spectral data ====================

# Read chemical reference analysis data
soilchem_tbl <- read_csv(file = "data/soilchem/soilchem_yamsys.csv")

# Join spectra tibble and chemical reference analysis tibble
spec_chem <- join_spc_chem(
  spc_tbl = soilspec_tbl_eth , chem_tbl = soilchem_tbl, by = "sample_id")


################################################################################
## Part 2: Run PLS regression models for different soil variables
################################################################################

# Example Partial Least Squares (PLS) Regression model for total Carbon (C)
pls_C <- pls_ken_stone(
  spec_chem = spec_chem,
  ratio_val = 1/3,
  variable = C,
  validation = TRUE,
  pc = 6,
  pls_ncomp_max = 6
)
```

# Details on functions, arguments, and input and output data structures

## `read_opus()` function

# Credits

I would like to thank the following people for the inspiration by concepts, code and packages:

* Antoine Stevens and Leonardo Ramirez-Lopez for their contributions to the [prospectr package](https://cran.r-project.org/web/packages/prospectr/index.html) and the
*Guide to Diffuse Reflectance Spectroscopy & Multivariate Calibration*
* Andrew Sila, Tomislav Hengl, and Thomas Terhoeven-Urselmans for the `read.opus()`
function from the [soil.spec](https://cran.r-project.org/web/packages/soil.spec/index.html) package developed at ICRAF.
* [Hadley Wickham](http://hadley.nz/) for his work and concepts on data science within R


