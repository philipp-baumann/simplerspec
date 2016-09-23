# Perform sampling for selection of reference analysis based on spectral PCA ---
#' @title Select a set of reference spectra to be measured by reference analysis
#' methods
#' @description Select a set of calibration spectra to develop spectral models.
#' Samples in this list will be analyzed using laboratory reference methods.
#' @param list_spectra list that contains preprocessed spectra and metadata
#' @param ratio_ref Ratio of desired reference samples to total sample number
#' @param pc Number of principal components (numeric). If pc < 1, the number
#' of principal components kept corresponds to the number of components
#' explaining at least (pc * 100) percent of the total variance.
#' @param print logical expression whether a plot (ggplot2) of sample selection
#' for reference analysis is shown in PCA space
#' @param validation Logical expression whether
#' calibration sampling is performed
#' (\code{TRUE} or \code{FALSE}).
#' @usage ken_stone(spec_chem, ratio_val, pc, print = TRUE,
#' validation = TRUE)
#' @export
select_ref_samples <- function(list_spectra, ratio_ref, pc = 2, print = TRUE) {
  pc_number <- eval(pc, envir = parent.frame())
  sel <- prospectr::kenStone(X = list_spectra$MIR0,
    k = round(ratio_ref * nrow(list_spectra$MIR0)), pc = substitute(pc_number))
  # Select metadata and spectra of reference samples based on row indices
  ref_metadata <- list_spectra$data_meta[sel$model, ]
  ref_spectra <- list_spectra$MIR0[sel$model, ]
  # Select metadata and spectra of prediction samples based on row indices
  pred_metadata <- list_spectra$data_meta[sel$model, ]
  pred_spectra <- list_spectra$MIR0[sel$model, ]
  # Create list of reference samples containing metadata and spectra
  ref_samples <- list(
    ref_metadata = ref_metadata,
    ref_spectra = ref_spectra
  )
  # Create list of prediction samples containing metadata and spectra
  pred_samples <- list(
    pred_metadata = pred_metadata,
    pred_spectra = pred_spectra
  )
  # Return sample list
  list(
    ref_samples = ref_samples,
    pred_samples = pred_samples
  )
}
