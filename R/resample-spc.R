#' @title Resample spectra from spectral data object (tibble)
#' @description Resamples spectra stored in tibble column after
#' gathering spectra by \code{simplerspec::gather_spc()}.
#' @param spc_tbl Spectra data in a tibble object (classes "tibble_df", "tbl"
#' and "data.frame"). The spectra tibble is expected to contain at least
#' the columns 'spc' (list-column with spectral matrices stored in a list) and
#' `wavenumbers` or `wavelengths` (list-column that contains list of x-axis values
#' corresponding to each spectrum in `spc` as wavenumber or wavelength in
#' numeric values).
#' @param x_unit Character vector of the x-axis unit of the spectra in
#' within \code{spc_tbl}. Default is \code{"wavenumber"}. A further possible
#' argument is \code{"wavelength"}, which is usually the unit reported for
#' VIS/NIR spectrometers such as ASD field spectrometers.
#' @param wn_lower Numeric (integer or float) value of lowest wavenumber.
#' This argument will only be used if \code{x_unit = "wavenumber"}. The value
#' will be used as starting value when creating the new wavenumber sequence
#' that the spectra will be resampled upon. Default value is 510.
#' @param wn_upper Numeric (integer or float) value of highest wavenumber.
#' This argument will only be used if \code{x_unit = "wavenumber"}. The value
#' will be used as last value when creating the new wavenumber sequence
#' that the spectra will be resampled upon. Default value is 3988.
#' @param wn_interval Numeric value (integer or float) value of the wavenumber
#' increment when creating the new wavenumber sequence that the spectra will be
#' resampled upon. Default value is 2.
#' @param wl_lower Numeric (integer or float) value of lowest wavelength.
#' This argument will only be used if \code{x_unit = "wavelength"}. The value
#' will be used as starting value when creating the new wavenumber sequence
#' that the spectra will be resampled upon. Default value is 350.
#' @param wl_upper Numeric (integer or float) value of highest wavelength.
#' This argument will only be used if \code{x_unit = "wavelength"}. The value
#' will be used as last value when creating the new wavenumber sequence
#' that the spectra will be resampled upon. Default value is 2500.
#' @param wl_interval Numeric value (integer or float) value of the wavelength
#' increment when creating the new wavenumber sequence that the spectra will be
#' resampled upon. This argument will only be used if \code{x_unit = "wavelength"}.
#' Default value is 1.
#' @return A spectra tibble containing the following additional columns
#' added to the input spectra tibble (`spc_tbl`): \code{wavenumbers_rs} or
#' \code{wavelengths_rs} as list-columns containing the resampled wavenumbers
#' or wavelengths, \code{spc_rs} list-column of resampled spectra returned
#' as list of data.tables (class "data.table" and "data.frame").
#' @export
resample_spc <- function(spc_tbl, x_unit = "wavenumber",
    wn_lower = 510, wn_upper = 3988, wn_interval = 2,
    wl_lower = 350, wl_upper = 2500, wl_interval = 1) {

  if (x_unit == "wavenumber" && "wavenumbers" %in% names(spc_tbl)) {
    # Create sequence of new wavenumbers
    wn_seq <- rev(seq(from = wn_lower, wn_upper, by = wn_interval))

    # Collect sequence of wavenumbers in list
    wavenumbers_rs <- rep(list(wn_seq), nrow(spc_tbl))
    names(wavenumbers_rs) <- names(spc_tbl$spc)

    # Resample all spectra in list column spc using prospectr
    spc_rs <- lapply(seq_along(spc_tbl$spc), function(i) {
      data.table::data.table(
        prospectr::resample(
          X = spc_tbl$spc[[i]], # spectral matrix to resample
          wav = spc_tbl$wavenumbers[[i]], # old wavenumbers
          new.wav = wn_seq # new wavenumbers
        )
      )
    }
    )
    names(spc_rs) <- names(spc_tbl$spc)

    # Add list of resampled spectra matrices to tibble spc_tbl
    spc_tbl %>% tibble::add_column(
      spc_rs = spc_rs,
      wavenumbers_rs = wavenumbers_rs
    )

  } else if (x_unit == "wavelength" && "wavelengths" %in% names(spc_tbl)) {
    # Create sequence of new wavelengths
    wl_seq <- seq(from = wl_lower, wl_upper, by = wl_interval)
    # Convert from wavelength (in nm) to wavenumbers (in cm^-1)
    wn_seq <- 10000000 / wl_seq

    # Collect sequence of wavelengths and wavenumbers in list
    wavelengths_rs <- rep(list(wl_seq), nrow(spc_tbl))
    wavenumbers_rs <- rep(list(wn_seq), nrow(spc_tbl))
    names(wavelengths_rs) <- names(spc_tbl$spc)
    names(wavenumbers_rs) <- names(spc_tbl$spc)

    # Resample all spectra in list column spc using prospectr
    spc_rs <- lapply(seq_along(spc_tbl$spc), function(i) {
      data.table::data.table(
        prospectr::resample(
          X = spc_tbl$spc[[i]], # spectral matrix to resample
          wav = spc_tbl$wavelengths[[i]], # old wavenumbers
          new.wav = wl_seq # new wavenumbers
        )
      )
    }
    )
    names(spc_rs) <- names(spc_tbl$spc)

    # Add list of resampled spectra matrices to tibble spc_tbl
    spc_tbl %>% tibble::add_column(
      spc_rs = spc_rs,
      wavelengths_rs = wavelengths_rs,
      wavenumbers_rs = wavenumbers_rs
    )
  } else {
    stop("Either columns 'wavenumbers' and 'wavelengths' are missing in the
      spectra tibble <spc_tbl> or argument x_unit has not been set to
      'wavenumber' or 'wavelength.")
  }

}
