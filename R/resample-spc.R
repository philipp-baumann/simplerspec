#' @title Resample spectra to new x-axis interval
#' @description Resamples (interpolates) different spectra types with
#' corresponding x-axis values that are both stored in list-columns of a spectra
#' tibble. A spectra tibble hosts spectra, x-unit vectors, metadata, and
#' further linked data with standardized naming conventions. Data input for
#' resampling can for example be generated with `simplerspec::gather_spc()`.
#' Resampling is a key harmonizing step to process and later model spectra
#' measured at different resolutions and spectral ranges (i.e., different
#' spectrometer devices and/or measurement settings).
#' @param spc_tbl Spectra data embedded in a tibble object (classes
#' `"tbl_df", "tbl", "data.frame"`). The spectra tibble needs to contain at
#' least of one of the the spectra columns `spc`, `spc_rs`, `spc_mean`,
#' `spc_nocomp`, `sc_sm`, `sc_rf`, or `spc_pre` (list-colums with spectral
#' `data.table`s), and `wavenumbers` or `wavelengths` (list-column with vectors
#' of x-axis values corresponding to each spectrum). The section *"Matching
#' spectrum type and corresponding x-axis type"* describes the spectra types
#' and corresponding x-axis types.
#' @param column_in Character vector of length 1L or symbol/name
#' specifying the name of list-column that contains the spectra to be resampled.
#' @param x_unit Character vector of length 1L specifying the measurement unit
#' of the x-axis values (list-column) of the input spectra in `spc_tbl`.
#' Possible values are `"wavenumber"` (default) or `"wavelength"`. Wavenumber
#' is a convenient unit of frequency in the mid-infrared spectral range,
#' where wavelength is often used as spatial period for the visible and
#' near-infrared range.
#' @param wn_lower Numeric value of lowest wavenumber. This argument will only
#' be used if `x_unit = "wavenumber"`. The value serves as starting value for
#' the new wavenumber sequence that the spectra will be resampled upon. Default
#' value is 500 (i.e., in reciprocal centimeters).
#' @param wn_upper Numeric value of highest wavenumber. This argument will only
#' be used if `x_unit = "wavenumber`. The value will be used as last value of
#' the new wavenumber sequence that the spectra will be resampled upon. Default
#' value is 4000 (i.e., in reciprocal centimeters).
#' @param wn_interval Numeric value of the wavenumber increment for the new
#' wavenumber sequence that the spectra will be resampled upon. Default value
#' is 2 (i.e., in reciprocal centimeters).
#' @param wl_lower Numeric value of lowest wavelength. This argument will only
#' be used if `x_unit = "wavelength"`. The value serves as starting value of
#' the new wavenumber sequence that the spectra will be resampled upon.
#' Default value is 350 (i.e. in nanometers).
#' @param wl_upper Numeric value of highest wavelength. This argument will only
#' be used if `x_unit = "wavelength"`. The value will be used as last value of
#' the new wavenumber sequence that the spectra will be resampled upon. Default
#' value is 2500 (i.e., in nanometers).
#' @param wl_interval Numeric value of the wavelength increment for the new
#' wavenumber sequence that the spectra will be resampled upon. This argument
#' will only be used if `x_unit = "wavelength"`. Default value is 1 (i.e., in
#' nanometers).
#' @param interpol_method Character of `"linear"` (default) or `"spline"` with
#' the interpolation method. `"spline"` uses a cubic spline to interpolate the
#' input spectra at given x-axis values to new equispaced x-axis intervals.
#' @return A spectra tibble (`spc_tbl`) containing two added list-columns:
#' * `spc_rs:` Resampled spectra as list of `data.table`s
#' * `wavenumbers_rs` or `wavelengths_rs`: Resampled x-axis values as list of
#'    numeric vectors
#' @section Matching spectrum type and corresponding x-axis type:
#' The combinations of input spectrum types (`column_in`) and
#' corresponding x-axis types are generated from a simple lookup list. The
#' following key-value(s) pairs can be matched at given key, which is the column
#' name from `column_in` containing the spectra.
#' * `"spc"` : `"wavenumbers"` or `"wavelengths"` (raw spectra)
#' * `"spc_rs"` : `"wavenumbers_rs"` or `"wavelengths_rs"`) (resampled spectra)
#' * `"spc_mean"` : `"wavenumbers_rs"` or `"wavelengths_rs"` (mean spectra)
#' * `"spc_nocomp"` `"wavenumbers"` or `"wavelengths"` (spectra prior
#'   atmospheric compensation)
#' * `"sc_sm" : c("wavenumbers_sc_sm", "wavelengths_sc_sm")` (single channel
#'   sample spectra)
#' * `"sc_rf" : c("wavenumbers_sc_rf", "wavelengths_sc_rf")` (single channel
#'   reference spectra)
#' * `"spc_pre" : "xvalues_pre"` (preprocessed spectra)
#' @export
resample_spc <- function(spc_tbl,
                         column_in = "spc",
                         x_unit = c("wavenumber", "wavelength"),
                         wn_lower = 500, wn_upper = 4000, wn_interval = 2,
                         wl_lower = 350, wl_upper = 2500, wl_interval = 1,
                         interpol_method = c("linear", "spline")) {
  # Capture user input as expressions (can be both of type character or symbol),
  # also called quoting; convert quosures to characters for later arg matching
  column_in <- rlang::enquo(column_in)
  column_in_chr <- rlang::quo_name(column_in)

  stopifnot(
    is.character(x_unit) && length(x_unit) > 0,
    is.numeric(wn_lower), is.numeric(wn_upper), is.numeric(wn_interval),
    is.numeric(wl_lower), is.numeric(wl_upper), is.numeric(wl_interval)
  )

  # Lookup list to match spectrum types and corresponding x-axis types
  spc_xaxis_types <- list(
    "spc" = c("wavenumbers", "wavelengths"), # raw/unprocessed
    "spc_rs" = c("wavenumbers_rs", "wavelengths_rs"), # resampled
    "spc_mean" = c("wavenumbers_rs", "wavelengths_rs"), # mean
    "spc_nocomp" = c("wavenumbers", "wavelengths"), # no atm. compensation
    "sc_sm" = c("wavenumbers_sc_sm", "wavelengths_sc_sm"), # single channel sample
    "sc_rf" = c("wavenumbers_sc_rf", "wavelengths_sc_rf"), # single channel reference
    "spc_pre" = rep("xvalues_pre", 2) # preprocessed
  )
  spctypes <- names(spc_xaxis_types)
  column_spc <- match.arg(column_in_chr, spctypes)

  x_unit <- match.arg(x_unit)
  switch(x_unit,
         wavenumber = {x_unit_int <- 1L},
         wavelength = {x_unit_int <- 2L})

  interpol_method <- match.arg(interpol_method)

  # Final selection of `x_unit` column name string from user input and lookup
  x_unit_sel <- spc_xaxis_types[[column_spc]][x_unit_int]

  # Both columns with X-values and input spectra need to be present in `spc_tbl`
  colnm <- colnames(spc_tbl)
  stopifnot(x_unit_sel %in% colnm, column_spc %in% colnm)

  # Extract list-column containing spectra
  spc_in_list <- dplyr::pull(spc_tbl, !!column_in)

  # Extract list-column containing x-axis values
  xvalues_in_list <- dplyr::pull(spc_tbl, !!x_unit_sel)

  # Automatically check the arrangement of the input x-Unit values;
  # often, it is convenient to have have a descending ordner of spectral columns
  # if the physical quantity of the x-axis is wavenumbers
  xvalue_order_chr <- purrr::map_chr(xvalues_in_list, seq_order)

  if (length(unique(xvalue_order_chr)) > 1L) {
    stop(
      glue::glue(
        "The column `{x_unit_sel}` which contains the list of X-values
        has both elements of ascending and descending order.
        * To resolve, you can split `spc_tbl` in a list of `spc_tbl`s
          with identical X-value vectors based on `group_by_col_hash()`,
          and apply `resample_spc()` separately to each list element.
        * Alternatively, you could fix the order of x-axis values
          for all input spectra and X-value vectors to all ascending or
          descending"),
      call. = FALSE)
  }
  xvalue_order <- xvalue_order_chr[1L]

  # Generate sequence of new x-axis values
  switch(x_unit_int,
         `1L` = {
           xvalues_out <- seq(from = wn_lower, to = wn_upper, by = wn_interval)
           x_unit_type_rs <- "wavenumbers_rs"
          },
         `2L` = {
           xvalues_out <- seq(from = wl_lower, to = wl_upper, by = wl_interval)
           x_unit_type_rs <- "wavelengths_rs"
         })

  if (xvalue_order == "descending") xvalues_out <- rev(xvalues_out)

  # Repeat sequence of new (resampled) x-axis values in list (for every obs.)
  xvalues_out_list <- rep(list(xvalues_out), nrow(spc_tbl))
  names(xvalues_out_list) <- names(spc_in_list)

  # Resample all spectra extracted from list-column `column_in` using prospectr
  spc_rs <- lapply(
    seq_along(spc_in_list),
    function(i) {
      data.table::data.table(
        prospectr::resample(
          X = spc_in_list[[i]], # spectral data.table to resample
          wav = xvalues_in_list[[i]], # old x-values vector
          new.wav = xvalues_out_list[[i]], # new x-values vector
          interpol = interpol_method
        )
      )
    }
  )
  names(spc_rs) <- names(spc_in_list)

  spc_tbl_out <-
    spc_tbl %>%
    tibble::add_column(
      spc_rs = spc_rs,
      !!x_unit_type_rs := xvalues_out_list
    )
  return(spc_tbl_out)
}

# Helper
seq_order <- function(x) ifelse(x[1L] < x[length(x)], "ascending", "descending")
