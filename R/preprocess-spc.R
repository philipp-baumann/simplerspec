#' @title Preprocess spectra from spectral data object (tibble)
#' @description Preprocesses spectra in tibble column by sample_id after
#' averaging spectra by \code{simplerspec::average_spc()}.
#' @export
preprocess_spc <- function(spc_tbl, select, column_in = "spc_mean") {

  # Convert list of spectral data.tables to one data.table
  spc_raw <- data.table::rbindlist(spc_tbl[column_in][[column_in]])

  # Perform preprocessing
  sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
    m = 0, p = 3, w = 9) # smoothing and averaging
  sg_1_w5 <- prospectr::savitzkyGolay(X = spc_raw,
    m = 1, p = 3, w = 5) # first derivative ***
  # Implement window size of 21, corresponds to ICRAF standard;
  # see e.g. Terhoeven-Urselmans et al. (2010)
  sg_1_w11 <- prospectr::savitzkyGolay(X = spc_raw,
    m = 1, p = 3, w = 11)
  sg_1_w15 <- prospectr::savitzkyGolay(X = spc_raw,
    m = 1, p = 3, w = 15)
  sg_1_w21 <- prospectr::savitzkyGolay(X = spc_raw,
    m = 1, p = 3, w = 21)
  sg_1_w35 <- prospectr::savitzkyGolay(X = spc_raw,
    m = 1, p = 3, w = 35)
  sg_1_w51 <- prospectr::savitzkyGolay(X = spc_raw,
    m = 1, p = 3, w = 51)
  # Second derivative Savitzky-Golay with a window size of 21 points
  sg_2_w21 <- prospectr::savitzkyGolay(X = spc_raw,
    m = 2, p = 3, w = 21)
  sg_2_w11 <- prospectr::savitzkyGolay(X = spc_raw,
    m = 2, p = 3, w = 11)
  # Savitzky-Golay (order 0) smoothing and derivative with a window size of
  # 21 points
  sg_0_1_w21 <- prospectr::savitzkyGolay(X = sg_0_w9,
    m = 1, p = 3, w = 21)
  # First derivative and window size of 13
  sg_1_w13 <- prospectr::savitzkyGolay(X = spc_raw,
    m = 1, p = 3, w = 13)
  # Second derivative
  sg_2 <- prospectr::savitzkyGolay(X = spc_raw,
    m = 2, p = 3, w = 5)
  # Calculate standard normal variate (SNV) after smoothing
  sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)
  sg_1_snv <- prospectr::standardNormalVariate(sg_1_w5) # added 2016-08-05
  # Standard normal variate (SNV) and first gap-segment derivative
  snv_gsd_m1_w11_s1 <- prospectr::gapDer(X = sg_0_snv, m = 1, w = 11, s = 1)
  snv_gsd_m1_w21_s5 <- prospectr::gapDer(X = sg_0_snv, m = 1, w = 21, s = 5)
  snv_gsd_m1_w31_s1 <- prospectr::gapDer(X = sg_0_snv, m = 1, w = 31, s = 5)
  snv_gsd_m1_w31_s5 <- prospectr::gapDer(X = sg_0_snv, m = 1, w = 31, s = 5)
  # Standard normal variate (SNV) and second gap-segment derivative
  snv_gsd_m2_w21_s1 <- prospectr::gapDer(X = sg_0_snv, m = 2, w = 21, s = 1)
  snv_gsd_m2_w31_s1 <- prospectr::gapDer(X = sg_0_snv, m = 2, w = 31, s = 5)
  snv_gsd_m2_w31_s5 <- prospectr::gapDer(X = sg_0_snv, m = 2, w = 31, s = 1)
  snv_gsd_m2_w51_s1 <- prospectr::gapDer(X = sg_0_snv, m = 2, w = 51, s = 1)
  snv_gsd_m2_w51_s5 <- prospectr::gapDer(X = sg_0_snv, m = 2, w = 51, s = 5)
  snv_gsd_m2_w5_s1 <- prospectr::gapDer(X = sg_0_snv, m = 2, w = 5, s = 1)
  # Gap-segement derivative; default option
  gsd <- prospectr::gapDer(X = spc_raw, m = 1, w = 11, s = 10)
  gsd_m1_w21_s21 <- prospectr::gapDer(X = spc_raw, m = 1, w = 21, s = 21)
  gsd_m1_w21_s1 <- prospectr::gapDer(X = spc_raw, m = 1, w = 21, s = 1)
  gsd_m1_w35_s21 <- prospectr::gapDer(X = spc_raw, m = 1, w = 35, s = 21)
  gsd_m1_w5_s4 <- prospectr::gapDer(X = spc_raw, m = 1, w = 5, s = 4)
  gsd_m1_w5_s21 <- prospectr::gapDer(X = spc_raw, m = 1, w = 5, s = 21)
  gsd_m1_w11_s21 <- prospectr::gapDer(X = spc_raw, m = 1, w = 11, s = 21)
  gsd_m1_w11_s5 <- prospectr::gapDer(X = spc_raw, m = 1, w = 11, s = 5)
  # 4th Gap-segment derivative
  gsd_m4_w21_s21 <- prospectr::gapDer(X = spc_raw, m = 4, w = 21, s = 21)
  # 2nd Gap-segment derivative
  gsd_m2_w21_s21 <- prospectr::gapDer(X = spc_raw, m = 2, w = 21, s = 21)
  # Continuum-removal
  cr <- prospectr::continuumRemoval(X = spc_raw,
    wav = as.numeric(colnames(spc_raw)), type = "A")
  # Select final preprocessing based on selection argument and
  # save matrix in data table
  pre <- select
  spc_pre <- data.table::as.data.table(get(pre))

  # Convert preprocessed spectra in data.table to list of data.table spectra
  spc_pre_list <- split(spc_pre, seq(nrow(spc_pre)))

  # Add list of preprocessed spectra to tibble
  spc_tbl <- tibble::add_column(spc_tbl, spc_pre = spc_pre_list)
}
