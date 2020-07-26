#' @title Average spectra in list-column by entries in grouping column
#' @description Average spectra in list-column of spectra tibble (`spc_tbl`) by
#' groups given in group column.
#' @param spc_tbl Tibble data frame containing at least the grouping column
#' given in argument `by` and input spectra given in list-column `column_in`.
#' @param by Character vector of length 1L or name/symbol that specifies the
#' column by which groups of spectra are averaged. Default is `"sample_id"`.
#' @param column_in Character vector of length 1L or or name/symbol that
#' specifies the list-column that contains the inputs spectra to be averaged.
#' Default is `"spc_rs"`, which are resampled spectra (i.e., resulting after
#' preceding `resample_spc()` step).
#' @return Spectra tibble data frame (class `"tbl_df"`, `"tbl"`, `"data.frame"`)
#' with a new list-column of column name `"spc_mean"` at the last position,
#' containing mean spectra with identical row replicates within the same
#' `by`-group.
#' @details For memory efficiency and subsequent modeling, consider slicing the
#' extra row copies of `spc_mean` resulting from `average_spc()` for example by
#' * `split(x = spc_tbl, f = spc_tbl$<by>) %>% lapply(., function(x) x x[1, ]) %>% do.call(., rbind)`
#' * `dplyr::group_by(spc_tbl, <by>) %>% dplyr::slice(1L)`
#' @import stats
#' @importFrom data.table data.table rbindlist setkey setDT := .SD
#' @importFrom rlang ensym quo_name
#' @export
average_spc <- function(spc_tbl, by = "sample_id", column_in = "spc_rs") {

  # Avoid R CMD check note: `"...no visible binding for global variable..."`
  spc_rs <- sample_id <- id <- NULL

  # Quote the symbol or the string supplied by the second and third argument
  column_in <- rlang::enquo(column_in)
  by <- rlang::enquo(by)

  # Combine rows of all resampled spectra in one data.table
  spc <- data.table::rbindlist(dplyr::pull(spc_tbl, !!column_in))

  # Add `id` group column to input spectra
  spc[, id := dplyr::pull(spc_tbl, !!by)] # spc_tbl[, by][[by]]

  # Average spectra, use `id` column as index for grouping
  data.table::setkey(spc, id)
  spc_mean <- spc[, lapply(.SD, mean), by = id]

  # Create new vector of group ID values from column `id`
  group_id_mean <- spc_mean[, id]
  # Delete sample_id column in data.table
  spc_mean_noid <- spc_mean[, id := NULL]

  # Create list of averaged spectra, one spectrum is one data.table
  # Use best performing alternative:
  # https://github.com/jennybc/row-oriented-workflows/blob/master/iterate-over-rows.md
  spc_mean_list <- stats::setNames(
    map(purrr::transpose(spc_mean_noid), data.table::as.data.table),
      group_id_mean
  )

  # Convert averaged spectra and sample_id to tibble
  spc_mean_tbl <- tibble::tibble(
    !!by := group_id_mean, spc_mean = spc_mean_list
  )
  # Join mean spectra tibble spc_tbl_mean to spc_tbl
  spc_tbl_out <- dplyr::left_join(spc_tbl, spc_mean_tbl,
    by = rlang::quo_name(by))

  return(spc_tbl_out)
}
