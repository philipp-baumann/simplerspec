#' Select every n-th spectral variable for all spectra and x-values in spectral
#' tibble (`spc_tbl`)
#'
#' @param spc_tbl Tibble data.frame containing spectra in list-column
#' @param lcol_spc List-column containing spectra, specified with column
#' name as symbols or 1L character vector.
#' @param lcol_xvalues List-column containing x-values, specified with
#' column name as symbols or 1L character vector.
#' @param every Every n-th spectral positions to keep as 1L integer vector.
#'
#' @return
#' @export
#'
#' @examples
select_spc_vars <- function(spc_tbl,
                            lcol_spc = spc_pre,
                            lcol_xvalues = xvalues_pre,
                            every = NULL) {
  lcol_spc <- rlang::enquo(lcol_spc)
  lcol_spc_nm <- rlang::quo_name(lcol_spc)
  lcol_xvalues <- rlang::enquo(lcol_xvalues)
  lcol_xvalues_nm <- rlang::quo_name(lcol_xvalues)

  stopifnot(tibble::is_tibble(spc_tbl))

  if (is.null(every)) {return(spc_tbl);
    message("Returning `spc_tbl` and keep all variables.")}

  spc_lst <- dplyr::pull(spc_tbl, !!lcol_spc)
  spc <- data.table::rbindlist(spc_lst)

  pos_sel <- seq(1L, ncol(spc), every)

  spc_sel <- spc[, ..pos_sel]

  spc_lst_out <- stats::setNames(
    map(purrr::transpose(spc_sel), data.table::as.data.table),
      names(spc_lst))

  xvalues <- dplyr::pull(spc_tbl, !!lcol_xvalues)
  xvalues_sel <- map(xvalues, ~ .x[pos_sel])

  dplyr::mutate(spc_tbl,
    !!lcol_spc_nm := spc_lst_out, !!lcol_xvalues_nm := xvalues_sel)
}
