#' @title Preprocess spectra
#' @description Use commonly used preprocessing algorithms on
#' the spectra.
#' @param list_spectra List that contains averaged spectra
#' in the list element called \code{MIR_mean}
#' @param select Character string that specifies the predefined
#' pretreatment options. Possible arguments are:
#' \code{select = "MIR0"} for Savitzky Golay smoothing filter
#' without derivative, \code{select = "MIR1"} for Savitky Golay
#' with first derivative, \code{select = "MIR2"} for Savitzky
#' Golay with second derivative, \code{select = "MIR0_snv"}
#' for Standard Normal Variate after Savitzky Golay without
#' derivative, and \code{select = "MIRb"} for
#' baseline correction.
#' @usage do_pretreatment(list_spectra, select)
#' @return list_spectra: List that contains preprocessed
#' spectra in element \code{MIR0}
#' @import hyperSpec
#' @export
do_pretreatment <- function(list_spectra, select) {
  MIR_mean <- NULL
  MIR_raw <- list_spectra$MIR_mean
  # Filter the data using the Savitzky and Golay smoothing filter
  # with a window size of 5 spectral variables and
  # a polynomial order of 3 (no differentiation)
  # p = polynomial order; plot variance vs polynomial order?
  # w = window size (must be odd)
  # m = m-th derivative of the polynomial coefficients
  # (0 = smoothing)
  MIR0 <- prospectr::savitzkyGolay(X = list_spectra$MIR_rs,
    m = 0, p = 3, w = 9) # smoothing and averaging
  MIR1 <- prospectr::savitzkyGolay(X = list_spectra$MIR_rs,
    m = 1, p = 3, w = 5) # first derivative ***
  # Implement window size of 21, corresponds to ICRAF standard;
  # see e.g. Terhoeven-Urselmans et al. (2010)
  MIR1_w21 <- prospectr::savitzkyGolay(X = list_spectra$MIR_rs,
    m = 1, p = 3, w = 21)
  # First derivative and window size of 13
  MIR1_w13 <- prospectr::savitzkyGolay(X = list_spectra$MIR_rs,
    m = 1, p = 3, w = 13)
  MIR2 <- prospectr::savitzkyGolay(X = list_spectra$MIR_rs,
    m = 2, p = 3, w = 5) # second derivative ***
  # Calculate standard normal variate (SNV) after smoothing
  MIR0_snv <- prospectr::standardNormalVariate(MIR0)
  MIR1_snv <- prospectr::standardNormalVariate(MIR1) # added 2016-08-05
  # Baseline correction
  # Compute baseline but first, create hyperSpec obj
  spc <- new("hyperSpec", spc = as.matrix(list_spectra$MIR_rs),
    wavelength = as.numeric(colnames(list_spectra$MIR_rs)))
  below <- hyperSpec::spc.fit.poly.below(
    fit.to = spc[, , 4000 ~ 900],
    apply.to = spc, npts.min = 20, poly.order = 2)
  spc_corr <- spc - below
  MIRb <- spc_corr[[]]
  pre <- select
  list_spectra$MIR_pre <- get(pre)
  return(list_spectra)
}
