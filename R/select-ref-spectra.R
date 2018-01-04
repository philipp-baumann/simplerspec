# Quick fix implementation of select_ref_samples using the tibble framework ----

# Perform sampling for selection of reference samples based on spectral PCA ----
#' @title Select a set of reference spectra to be measured by reference analysis
#' methods
#' @description Select a set of calibration spectra to develop spectral models.
#' Samples in this list will be analyzed using laboratory reference methods.
#' @param spc_tbl Spectra as tibble objects that contain preprocessed spectra
#' @param ratio_ref Ratio of desired reference samples to total sample number
#' @param pc Number of principal components (numeric). If pc < 1, the number
#' of principal components kept corresponds to the number of components
#' explaining at least (pc * 100) percent of the total variance.
#' @param print logical expression whether a plot (ggplot2) of sample selection
#' for reference analysis is shown in PCA space
#' (\code{TRUE} or \code{FALSE}).
#' @usage select_ref_spc(spc_tbl, ratio_ref, pc, print = TRUE)
#' @export
select_ref_spc <- function(spc_tbl, ratio_ref = 0.15, pc = 2,
  print = TRUE) {

  # Avoid `R CMD check` NOTE: `no visible binding for global variable ...`
  PC1 <- PC2 <- type <- NULL

  if(tibble::is_tibble(spc_tbl)) {
    # Slice based on sample_id if spectral data is in tibble class
    spc_tbl <- dplyr::group_by(spc_tbl, rlang::UQ(rlang::sym("sample_id"))) %>%
      dplyr::slice(1L)
    # Bind list of data.tables in list-column spc_pre to one data table
    # containing spectral data
    spc_pre <- as.matrix(data.table::rbindlist(spc_tbl$spc_pre))
  }
  # Perform Kennard-Stone calibration sampling ---------------------------------
  sel <- prospectr::kenStone(X = spc_pre,
    k = round(ratio_ref * nrow(spc_pre)), pc = pc)
  # Select spectra tibble of reference samples based on row indices
  spc_ref <- spc_tbl[sel$model, ]
  # Select spectra tibble of prediction samples based on row indices
  spc_pred <- spc_tbl[-sel$model, ]

  # Prepare data for ggplot graphs of reference and prediction sample PC score
  # plots (PC1 and PC2) --------------------------------------------------------
  sel_df_ref <- data.frame(sel$pc[sel$model, 1:2])
  sel_df_ref$type <- as.factor(
    rep("reference analysis", nrow(sel_df_ref))
  )
  sel_df_pred <- data.frame(sel$pc[-sel$model, 1:2])
  # Create type column for visually differentiate reference and prediction
  # samples
  sel_df_pred$type <- as.factor(
    rep("model prediction", nrow(sel_df_pred)))
  # Bind rows of reference and prediction PC scores data frames
  sel_df <- rbind(sel_df_ref, sel_df_pred)
  # Compute ratio needed to make the figure square
  ratio <- with(sel_df, diff(range(PC1)) / diff(range(PC2)))
  # Create spectra PC score plots ----------------------------------------------
  p_pca <- ggplot2::ggplot(data = sel_df) +
    ggplot2::geom_point(
      ggplot2::aes(x = PC1, y = PC2, shape = type), size = 4) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::scale_shape_manual(values=c(19, 1)) +
    ggplot2::scale_colour_manual(values=c("black", "red")) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank())
  # Print reference and prediction samples in PC1 and PC2
  if (print == TRUE) {
    p_pca
  }
  # Return spectral tibbles for reference spectra (spc_ref),
  # prediction spectra (spc_pr) and ggplot object of score plots (p_pca)
  list(
    spc_ref = spc_ref,
    spc_pred = spc_pred,
    p_pca = p_pca
  )
}



