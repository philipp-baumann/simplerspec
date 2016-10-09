#' @title Resample spectra stored to new
#' @description
#' @param list_spectra List of spectra and metadata
#' @param wn_lower Numerical value for lowest  wavenumber in sampling interval
#' @param wn_upper Numerical value for highest wavenumber in sampling interval
#' @export
resample_spectra <- function(
  list_spectra, wn_lower = 510, wn_upper = 3988, wn_interval = 2)
  {
  # Create sequence of new wavenumbers
  wn_seq <- rev(seq(from = wn_lower, wn_upper, by = wn_interval))
  list_spectra$MIR_rs <- prospectr::resample(
    X = list_spectra$MIR_mean, # spectral matrix to resample
    wav = as.numeric(colnames(list_spectra$MIR_mean)), # old wavenumbers
    new.wav = wn_seq # new wavenumbers
    )
  return(list_spectra)
}
