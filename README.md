<a href="https://github.com/philipp-baumann/simplerspec/blob/master/simplerspec-logo.png"><img src="https://github.com/philipp-baumann/simplerspec/blob/master/simplerspec-logo.png" width="300"/></a>

Travis-CI build status (Debian Linux) [![Travis Build Status](https://travis-ci.org/philipp-baumann/simplerspec.svg?branch=master)](https://travis-ci.org/philipp-baumann/simplerspec)

AppVeyor build status (Windows) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/philipp-baumann/simplerspec?branch=master&svg=true)](https://ci.appveyor.com/project/philipp-baumann/simplerspec)

# Short description

The simplerspec package aims to facilitate spectra and additional data handling and model development for spectroscopy applications such as infrared soil spectroscopy. Different helper functions are designed to create a 
data and modeling workflow. Data inputs and outputs are stored in common S3 `R` objects (`lists` and `data frames`), using in addition `data.table` and `tibble` extensions. The following features are covered in the current version of the package:

1. `read_opus_univ()`: Read spectra and metadata from Bruker OPUS binary files into R list
2. `gather_spc()`: Gather spectra and metadata from list into a tibble object (list-columns)
4. `resample_spc()`: Resample spectra to new wavenumber intervals
2. `average_spc()`: Average spectra for replicate scans
5. `preprocess_spc()`: Perform pre-processing of spectra
6. `join_chem_spc()`: Join chemical and spectral data sets by `sample_id`
7. `plot_spc_ext()`: Extended spectral plotting; e.g. group spectra using different
panels or color spectra based on chemical reference values to explore trends.
8. `fit_pls()`: Perform model tuning and evaluation based on Partial Least Squares (PLS) regression
9. `select_ref_spc()`: Select a set of reference samples to measured by
traditional analysis methods when no a priori sample data except spectra are 
available (based on Kennard-Stones sampling)
10. `predict_from_spc()`: Predict multiple chemical properties from a list of calibrated models and new soil spectra
11. `assess_multimodels()`: Assess model performance given multiple pairs of predicted and measured variables.

# Projects using simplerspec

* [Spectral platform for soil samples of the Democratic Republic of Congo](https://sae-interactive-data.ethz.ch/simplerspec.drc/)

# Cheatsheet

<a href="https://github.com/philipp-baumann/spc-proc-concepts/blob/master/img/simplerspec_cheatsheet_crop.pdf"><img src="https://github.com/philipp-baumann/spc-proc-concepts/blob/master/img/simplerspec_cheatsheet.png" width="630"/></a>

# Installation

The newest version of the package is available on this GitHub repository. Note that the package is still under development. If you find bugs you are highly welcome to report issues (write me an [email](mailto:philipp.baumann@usys.ethz.ch) or create an [issue](https://github.com/philipp-baumann/simplerspec/issues)). You can install `simplerspec` using the devtools package.

```R
# Uncomment and run the below line if you have not yet installed
# the devtools package
# install.packages("devtools")
# Install the simplerspec package from the github repository
# (https://github.com/philipp-baumann/simplerspec)
devtools::install_github("philipp-baumann/simplerspec")
```

## Special installation note for Windows 8 and R versions 3.3 and 3.4

For some Windows versions with recent R versions (3.3 and 3.4), there 
might be an error message that the `Rcpp` package can not be installed because
there is no precompiled binary (packaging up) of the `Rcpp` package available on CRAN. Because the `Rcpp` package contains C++ code, the package needs compilation.
The compiler is supplied in the R tools (contains GCC 4.9.3 and Mingw-W64 V3).
First, you need to download and install the latest R tools version from [here](https://cran.r-project.org/bin/windows/Rtools/). Then, you need to 
install `Rcpp` from source provided on CRAN by 

```R
# install.packages("Rcpp", type = "source")
```


# Motivation and key concepts

The functions are built to work in a pipeline and cover commonly used procedures for spectral model development. Many R packages are available to do tasks in spectral modeling such as pre-processing of spectral data. The motivation to create this package was:

1. Avoid repetition of code in model development (common source of errors).
2. Provide a reproducible data analysis workflow for FT-IR spectroscopy.
3. R packages are an ideal way to organize and share R code.
4. Make soil FT-IR spectroscopy modeling accessible to people that have basic R knowledge.
5. Provide an integrated data-model framework that features tidy data structures designed for both user-friendly printing and efficient data processing.

This package builds mainly upon functions from the following R packages:

* `prospectr `: Various utilities for pre-processing and sample selection based on spectroscopic data. An introduction to the package with examples can be found [here](http://antoinestevens.github.io/prospectr/).
* `plyr` and `dplyr `: Fast data manipulation tools with an unified interface. See [here](https://github.com/hadley/dplyr) for details.
* `ggplot2 `: Alternative plotting system for R, based on the grammar of graphics. See [here](http://ggplot2.org/).
* `caret `: Classification and regression training. A set of functions that attempt to streamline the process for creating predictive models. See [here](http://topepo.github.io/caret/index.html) for details.

Consistent and reproducible data and metadata management is an important prerequisite for spectral model development. Therefore, simplerspec functions are based on storing spectral data and related data in R data structures which keep related data in rows. Every row representing an observation contains data related to a single spectral measurement. Simplerspec functions uses tibble data frames as principal data structures because they allow to store lists within the well-known data frame structures. Lists are flexible data structures and can e.g. contain other lists, vectors, data.frames, or matrices.

List-columns features provided within the tibble framework are an excellent base to work with functional programming tools in R, which allows to efficiently write code. 
Simplerspec internally uses popular functional programming extension tools provided
by the `purrr` package for processing and transforming spectra. 
For learning more, I would recommend
[this nice purrr list-column tutorial](https://jennybc.github.io/purrr-tutorial/ls13_list-columns.html) 
provided by Jenny Brian. Further, simplerspec well integrates with the 
data processing API provided by the dplyr package, which makes spectroscopic
analysis tidy and easy to understand.

# Example workflow

Bruker FTIR spectrometers produce binary files in the OPUS format that can contain different types of spectra and many parameters such as instrument type and settings that were used at the time of data acquisition and internal processing (e.g. Fourier transform operations). Basically, the entire set of setup measurement parameters, selected spectra, supplementary metadata such as the time of measurement are written into OPUS binary files. In contrast to simple text files that contain only plain text with a defined character encoding, binary files can contain any type of data represented as sequences of bytes (a single byte is sequence of 8 bits and 1 bit either represents 0 or 1).

Simplerspec comes with reader function `read_opus_univ()` that is intended to be a universal Bruker OPUS file reader that extracts spectra and key metadata from files. Usually, one is mostly interested to extract the final absorbance spectra (shown as *AB* in the OPUS viewer software).

```R
# Load simplerspec package for spectral model development wrapper functions
library(simplerspec)
# Load tidyverse packages: loads packages frequently used for data manipulation,
# data tidying, import, and plotting
library(tidyverse)

################################################################################
## Part 1: Read and pre-process spectra, read chemical data, and join
## spectral and chemical data sets
################################################################################

## Read spectra in list ========================================================

# List of OPUS binary spectra files
lf <- dir("data/spectra/soilspec_eth_bin", full.names = TRUE)

# Read spectra from files into R list
spc_list <- read_opus_univ(fnames = lf, extract = c("spc"))
# Returns messages:
#> Extracted spectra data from file: <BF_lo_01_soil_cal.0>
#> Extracted spectra data from file: <BF_lo_01_soil_cal.1>
#> Extracted spectra data from file: <BF_lo_01_soil_cal.2>
#> Extracted spectra data from file: <BF_lo_02_soil_cal.0>
#> ...
```

Pipes can make R code more readable and allows step-wise data processing
when developing spectral models. The pipe operator (`%>%`, called "then") is a new operator in R that was introduced
with the magrittr package. It facilitates readability of code
and avoids to type intermediate objects. The basic behavior of
the pipe operator is
that the object on the left hand side is passed as the first argument
to the function on the right hand side. When loading the tidyverse packages, the
pipe operator is attached to the current R session.
More details can be found [here](https://github.com/smbache/magrittr).

The model development process can be quickly coded as the example below illustrates:

```R
## Spectral data processing pipe ===============================================

soilspec_tbl <- spc_list %>%
  # Gather list of spectra data into tibble data frame
  gather_spc() %>% 
  # Resample spectra to new wavenumber interval
  resample_spc(wn_lower = 500, wn_upper = 3996, wn_interval = 2) %>%
  # Average replicate scans per sample_id
  average_spc() %>%
  # Preprocess spectra using Savitzky-Golay first derivative with a window size
  # of 21 points
  preprocess_spc(select = "sg_1_w21")
  
soilspec_tbl
# A tibble: 284 x 11
#>                                 unique_id             file_id         sample_id
#>                                     <chr>               <chr>             <chr>
#> 1 BF_lo_01_soil_cal.0_2015-11-06 14:34:10 BF_lo_01_soil_cal.0 BF_lo_01_soil_cal
#> 2 BF_lo_01_soil_cal.1_2015-11-06 14:38:14 BF_lo_01_soil_cal.1 BF_lo_01_soil_cal
#> 3 BF_lo_01_soil_cal.2_2015-11-06 14:40:55 BF_lo_01_soil_cal.2 BF_lo_01_soil_cal
#> 4 BF_lo_02_soil_cal.0_2015-11-06 17:27:55 BF_lo_02_soil_cal.0 BF_lo_02_soil_cal
#> 5 BF_lo_02_soil_cal.1_2015-11-06 17:30:19 BF_lo_02_soil_cal.1 BF_lo_02_soil_cal
#> 6 BF_lo_02_soil_cal.2_2015-11-06 17:32:47 BF_lo_02_soil_cal.2 BF_lo_02_soil_cal
#> 7 BF_lo_03_soil_cal.0_2015-11-09 11:32:55 BF_lo_03_soil_cal.0 BF_lo_03_soil_cal
#> 8 BF_lo_03_soil_cal.1_2015-11-09 11:35:26 BF_lo_03_soil_cal.1 BF_lo_03_soil_cal
#> 9 BF_lo_03_soil_cal.2_2015-11-09 11:38:08 BF_lo_03_soil_cal.2 BF_lo_03_soil_cal
#> 10 BF_lo_04_soil_cal.0_2015-11-06 10:36:13 BF_lo_04_soil_cal.0 BF_lo_04_soil_cal
#> # ... with 274 more rows, and 8 more variables: spc <list>, wavenumbers <list>,
#> #   metadata <list>, spc_rs <list>, wavenumbers_rs <list>, spc_mean <list>,
#> #   spc_pre <list>, xvalues_pre <list>
  

## Read chemical reference data and join with spectral data ====================

# Read chemical reference analysis data
soilchem_tbl <- read_csv(file = "data/soilchem/soilchem_yamsys.csv")
#> Parsed with column specification:
#> cols(
#>   .default = col_double(),
#>   sample_ID = col_character(),
#>   country = col_character(),
#>   site = col_character(),
#>   material = col_character(),
#>   site_comb = col_character()
#> )
#> See spec(...) for full column specifications.

# Join spectra tibble and chemical reference analysis tibble
spec_chem <- join_spc_chem(
  spc_tbl = soilspec_tbl, chem_tbl = soilchem_tbl, by = "sample_id")
#> Joining, by = "sample_id"

################################################################################
## Part 2: Run PLS regression models for different soil variables
################################################################################

# Example Partial Least Squares (PLS) Regression model for total Carbon (C)
# Use repeated k-fold cross-validation to tune the model (choose optimal 
# number of PLS components) and estimate model performance on hold-out 
# predictions of the finally chosen model (model assessment).
# This allows to use the entire set for both model building and evaluation;
# recommended for small data sets
pls_C <- fit_pls(
  # remove rows with NA in the data
  spec_chem = spec_chem[!is.na(spec_chem$C), ],
  response = C,
  evaluation_method = "resampling",
  tuning_method = "resampling",
  resampling_method = "rep_kfold_cv",
  pls_ncomp_max = 7 # maximum number of PLS components tested during tuning
) 
```

# Package help

After successfully installing simplerspec, you can use the R build-in help
using `?simplerspec::<fun_name>`

# Credits

I would like to thank the following people for the inspiration by concepts, code and packages:

* Antoine Stevens and Leonardo Ramirez-Lopez for their contributions to the [prospectr package](https://cran.r-project.org/web/packages/prospectr/index.html) and the
*Guide to Diffuse Reflectance Spectroscopy & Multivariate Calibration*
* Andrew Sila, Tomislav Hengl, and Thomas Terhoeven-Urselmans for the `read.opus()`
function from the [soil.spec](https://cran.r-project.org/web/packages/soil.spec/index.html) package developed at ICRAF.
* [Hadley Wickham](http://hadley.nz/) for his work and concepts on data science within R
* [Max Kuhn](https://github.com/topepo) for the creation of the caret package and for his excellent teaching
materials on [applied predictive modeling](http://appliedpredictivemodeling.com/blog/)
