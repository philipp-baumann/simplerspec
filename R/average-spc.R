#' @title Average spectra from spectral data object (tibble)
#' @description Averages spectra in tibble column by sample_id after
#' resampling spectra by \code{simplerspec::resample_spc()}.
#' @export
average_spc <- function(spc_tbl) {

  # Avoid R CMD check note: `no visible binding for global variable`
  spc_rs <- sample_id <- NULL

  # Combine rows of all resampled spectra in one data.table
  spc <- data.table::rbindlist(spc_tbl$spc_rs)

  # Add sample_id column to resampled spectra
  spc[, sample_id:=spc_tbl$sample_id]

  # Average spectra, use sample_id as index for grouping
  data.table::setkey(spc, sample_id)
  spc_mean <- spc[, lapply(.SD, mean), by = sample_id]

  # Create vector of sample_id from column sample_id in spc_mean
  sample_id_mean <- spc_mean$sample_id
  # Delete sample_id column in data.table
  spc_mean_noid <- spc_mean[, sample_id := NULL]

  # Create list of averaged spectra, one spectrum is one data.table
  # Method split.data.table is not yet available in data.table 1.9.7
  # Wait for v2.0.0
  spc_mean_list <- setNames(
    split(spc_mean_noid, seq(nrow(spc_mean_noid))),
    sample_id_mean
  )

  # Convert averaged spectra and sample_id to tibble
  spc_mean_tbl <- tibble::tibble(
    sample_id = sample_id_mean, spc_mean = spc_mean_list
  )
  # Join mean spectra tibble spc_tbl_mean to spc_tbl
  spc_tbl <- dplyr::left_join(spc_tbl, spc_mean_tbl)
}
