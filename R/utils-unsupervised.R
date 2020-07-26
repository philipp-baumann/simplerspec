# Helper functions to compute principal component analysis on spectral
# data and to append scores and importance -------------------------------------

pca_append_scores <- function(spc_tbl, y = "spc_pre",
                              slice = TRUE, slice_by = sample_id,
                              select_comps = 1:2,
                              scale = TRUE, center = TRUE,
                              .unnest = NULL) {
  var <- rlang::enquo(y)
  var_nm_pca <- paste0(rlang::quo_name(var), "_pca_scores")
  slice_by <- rlang::enquo(slice_by)

  if (slice) {
    spc_tbl <- spc_tbl %>% dplyr::group_by(!!slice_by) %>% dplyr::slice(1L)
  }

  sample_id <- spc_tbl %>% dplyr::pull(!!slice_by)

  # Pull the list of preprocessed spectra data.tables into one data.table
  spc <- data.table::rbindlist(dplyr::pull(spc_tbl, !!var))

  # Perform a principal component analysis
  spc_pca <- stats::prcomp(spc, scale = scale, center = center)

  # Extract PCA scores for selected principal components
  spc_pca_scores <- dplyr::as_data_frame(spc_pca$x[, select_comps])
  # Convert data frame of scores to list column
  spc_pca_scores <- stats::setNames(
    split(spc_pca_scores, seq_len(nrow(spc_pca_scores))),
    sample_id)
  ncomp <- ncol(spc_pca$rotation)

  # Add list-column with PCA scores to the spectral tibble object
  spc_tbl_pca <- tibble::add_column(spc_tbl, !!var_nm_pca := spc_pca_scores)
  if (!is.null(.unnest)) {
    # Convert variable name of new column to symbol
    var_sym_pca <- rlang::sym(var_nm_pca)
    spc_tbl_pca <- spc_tbl_pca %>% tidyr::unnest(!!var_sym_pca)
  }
  spc_tbl_pca
}

#' @noRd
pca_append_importance <- function(spc_tbl, y = "spc_pre",
                                  slice = TRUE, slice_by = "sample_id",
                                  select_comps = 1:2,
                                  scale = TRUE, center = TRUE) {
  var <- rlang::enquo(y)
  slice_by <- rlang::enquo(slice_by)

  if (slice) {
    spc_tbl <- spc_tbl %>% dplyr::group_by(!!slice_by) %>% dplyr::slice(1L)
  }

  # Pull the list of preprocessed spectra data.tables into one data.table
  spc <- data.table::rbindlist(dplyr::pull(spc_tbl, !!var))

  # Perform a principal component analysis
  spc_pca <- stats::prcomp(spc, scale = scale, center = center)
  ncomp <- ncol(spc_pca$rotation)

  # // pb 20180509: broom::tidy(spc_pca) returns error; this is a bug in
  # broom:::tidy.prcomp(); does not work if data table without row names
  # Extract variance explained
  importance_measures <- c("sd", "var_prop", "var_cum")
  # "sd" := "Standard deviation"; "var_prop" := "Proportion of Variance";
  # "var_cum" := "Cumulative Proportion"
  spc_pca_varexpl <- broom::fix_data_frame(
    t(summary(spc_pca)$importance),
    newnames = importance_measures, newcol = "PC") %>%
    tibble::add_column(ncomp = 1:ncomp, .before = 1) %>%
    tibble::as_tibble()

  # Return list of spectral tibble with pc scores attached, and
  # the PCA importance measures in a tidy data frame
  spc_pca_varexpl
}
