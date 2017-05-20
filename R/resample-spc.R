#' @title Resample spectra from spectral data object (tibble)
#' @description Resamples spectra stored in tibble column after
#' gathering spectra by \code{simplerspec::gather_spc()}.
#' @export
resample_spc <- function(spc_tbl, wn_lower = 510, wn_upper = 3988,
    wn_interval = 2, x_unit = "wavenumber",
    wl_lower = 350, wl_upper = 2500, wl_interval = 1) {

  if (x_unit == "wavenumber" & "wavenumbers" %in% names(spc_tbl)) {
    # Create sequence of new wavenumbers
    wn_seq <- rev(seq(from = wn_lower, wn_upper, by = wn_interval))

    # Collect sequence of wavenumbers in list
    wavenumbers_rs <- rep(list(wn_seq), nrow(spc_tbl))
    names(wavenumbers_rs) <- names(spc_tbl$spc)

    # Resample all spectra in list column spc using prospectr
    spc_rs <- lapply(names(spc_tbl$spc), function(nm) {
      data.table::data.table(
        prospectr::resample(
          X = spc_tbl$spc[[nm]], # spectral matrix to resample
          wav = spc_tbl$wavenumbers[[nm]], # old wavenumbers
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

  } else if (x_unit = "wavelength" & "wavelengths" %in% names(spc_tbl)) {
    # Create sequence of new wavelengths
    wl_seq <- rev(seq(from = wl_lower, wl_upper, by = wl_interval))

    # Collect sequence of wavelengths in list
    wavelengths_rs <- rep(list(wl_seq), nrow(spc_tbl))
    names(wavelengths_rs) <- names(spc_tbl$spc)

    # Resample all spectra in list column spc using prospectr
    spc_rs <- lapply(names(spc_tbl$spc), function(nm) {
      data.table::data.table(
        prospectr::resample(
          X = spc_tbl$spc[[nm]], # spectral matrix to resample
          wav = spc_tbl$wavelengths[[nm]], # old wavenumbers
          new.wav = wl_seq # new wavenumbers
        )
      )
    }
    )
    names(spc_rs) <- names(spc_tbl$spc)

    # Add list of resampled spectra matrices to tibble spc_tbl
    spc_tbl %>% tibble::add_column(
      spc_rs = spc_rs,
      wavelengths_rs = wavelengths_rs
    )
  } else {
    stop("Either columns 'wavenumbers' and 'wavelengths' are missing in the
      spectra tibble <spc_tbl> or argument x_unit has not been set to
      'wavenumber' or 'wavelength.")
  }

}
