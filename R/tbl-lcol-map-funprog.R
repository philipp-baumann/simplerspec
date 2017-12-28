# Split a tibble data.frame into a list of tibbles by a group column
#' @export
split_df2l <- function(tbl_df, group) {
  split(tbl_df, tbl_df[, group])
}

# Extract multiple tibble list columns, row bind them separately into
# single data tables and return a list of data.tables
#' @title Extract multiple tibble list-columns and return data as list of
#' data.tables
#' @description Extract multiple tibble list columns, row bind them separately
#' into single data tables and return a list of data.tables.
#' @param spc_tbl Spectral tibble (data frame) with spectral data contained
#' in list-columns
#' @param lcols Character vector containing names of list-columns to be
#' extracted into a list of data.tables
#' @return List of data.tables. Each element is a data.table derivied from a
#' list-column specified in the \code{lcols} argument.
#' @import purrr
#' @export
extract_lcols2dts <- function(spc_tbl, lcols) {
  # Below code is first part of simplerspec::bind_lcols_dts
  # todo: add warning for lcols not present in spc_tbl
  which_bind <- colnames(spc_tbl) %in% lcols
  lcols_to_bind <- colnames(spc_tbl)[which_bind]
  names(lcols_to_bind) <- lcols_to_bind
  dts <- map(lcols_to_bind,
    function(y) {
      if (is.list(spc_tbl[[y]])) {
        # todo: Test if number of columns is equal in each data.frame or matrix
        # of the list-(column); if not, return a comprehensible error
        data.table::data.table(do.call(rbind, spc_tbl[[y]]))
      } else if (is.atomic((spc_tbl[[y]]))) {
        data.table::data.table(spc_tbl[, y])
      }
    }
  )
}
