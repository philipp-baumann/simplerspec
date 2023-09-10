<!-- NEWS.md is maintained by https://cynkra.github.io/fledge, do not edit -->

# simplerspec 0.2.0

- Add new example data set `soilspec_yamsys`.
- `gather_spc()`: Consolidate the documentation with details on how data are matched from list and gathered into a spectra tibble.

# simplerspec 0.1.0.9001

* `resample_spc()` now supports flexible spectra and x-axis types as inputs. Its
  interface has been carefully augmented without breaking previous 
  functionality ([#9](https://github.com/philipp-baumann/simplerspec/issues/9)
  * New argument `column_in` specifies the string or name (unquoting support)
    of the input column that contains the list of spectra. The following 
    spectrum types, which are automatically matched against the list-column that
    contains the corresponding x-unit value vectors, are currently supported:
    `spc` (raw or unprocessed spectra), `spc_rs` (resampled spectra),
    `spc_mean` (mean spectra), `spc_nocomp` (spectra prior atmospheric 
    compensation), `sc_sm` (single channel sample spectra), `sc_rf` (single 
    channel reference spectra), `spc_pre` (preprocessed spectra).
  * New argument `interpol_method` specifying the interpolation method is 
    introduced. Default is `"linear"` to achieve identical results with both
    prospectr v0.1.0 and v0.2.0. The current CRAN prospectr v0.2.0 has changed
    the default of `interpol` to `"spline"`. The previous `resample_spc()` 
    unfortunatelty did not explicitly state the method internally, and relied 
    on the default instead. The measures taken ensure downward compatibility of 
    `resample_spc()` with previous versions of prospectr and simplerspec.
  * The arguments gain more defensive checks inside the function (supplied types
    and presence of objects in spectra).
  * The function components and the help are updated accordingly. Clearer 
    vocabulary to describe the functionality and more consistent terminology for
    physical quantities and R objects are used.
  
* Add UTF-8 support to DESCRIPTION because roxygen2 version 7.1.0 requires it.


# simplerspec 0.1.0.9000

* Start using Kirill's `{fledge}` for tracking and communicating the simplerspec
  development process in `NEWS`.


# simplerspec 0.1.0

* Added a `NEWS.md` file to track changes to the package

# simplerspec 0.1.0.1

* `read_opus_bin_univ()`: Add support for Bruker files that have undefined `PLF` value (`:= NULL`)
