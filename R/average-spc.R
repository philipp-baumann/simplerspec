#' @title Average spectra in spectral tibble data frame
#' @description Average spectra in tibble by a group column after
#' resampling spectra by \code{simplerspec::resample_spc()}.
#' @param spc_tbl Tibble data.frame containing spectra in list-column
#' \code{spc_rs}. This list-column is created when resampling spectra
#' with \code{resample_spc()}
#' @param by Character vector of length 1; specifies column by which spectra
#' are averaged. Default is \code{"sample_id"}.
#' @return Tibble data frame (class `tbl_df`) with mean spectra appended
#' as list-column named \code{spc_mean}.
#' @import stats
#' @importFrom data.table data.table rbindlist setkey setDT := .SD
#' @importFrom rlang ensym quo_name
#' @export
average_spc <- function(spc_tbl, by = "sample_id", column_in = "spc_rs") {

  # Avoid R CMD check note: `no visible binding for global variable`
  spc_rs <- sample_id <- id <- NULL

  column_in <- rlang::enquo(column_in)

  # Combine rows of all resampled spectra in one data.table
  spc <- data.table::rbindlist(dplyr::pull(spc_tbl, !!column_in))

  # Add sample_id column to resampled spectra
  spc[, id := spc_tbl[, by][[by]]]

  # Average spectra, use sample_id as index for grouping
  data.table::setkey(spc, id)
  spc_mean <- spc[, lapply(.SD, mean), by = id]

  # Create vector of sample_id from column sample_id in spc_mean
  sample_id_mean <- spc_mean[, id]
  # Delete sample_id column in data.table
  spc_mean_noid <- spc_mean[, id := NULL]

  # Create list of averaged spectra, one spectrum is one data.table
  # Use best performing alternative:
  # https://github.com/jennybc/row-oriented-workflows/blob/master/iterate-over-rows.md
  spc_mean_list <- stats::setNames(
    purrr::transpose(spc_mean_noid),
      sample_id_mean
  )

  # Quote the symbol or the string supplied to by argument
  by <- ensym(by)

  # Convert averaged spectra and sample_id to tibble
  spc_mean_tbl <- tibble::tibble(
    !! by := sample_id_mean, spc_mean = spc_mean_list
  )
  # Join mean spectra tibble spc_tbl_mean to spc_tbl
  spc_tbl <- dplyr::left_join(spc_tbl, spc_mean_tbl, by = quo_name(by))
}
