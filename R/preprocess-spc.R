#' @title Preprocess spectra
#' @description Preprocesses spectra in tibble column by sample_id after
#' averaging spectra by \code{simplerspec::average_spc()}.
#' @param spc_tbl Tibble that contains spectra to be preprocessed within
#' a list-column.
#' @param select Character vector of predefined preprocessing options to be
#' applied to the spectra list-column specified in \code{column_in}.
#' Common prefined values are stated as abbreviated preprocessing methods and
#' options such as \code{"sg_1_w21"}, where \code{"sg"} stands for
#' Savitzky-Golay and \code{1} for first derivative and \code{"w21"}
#' for a window size of 21 points.
#' @param column_in Character vector of single list-column in \code{spc_tbl} that
#' contain list of spectra (1 row matrix) to be processed by function supplied
#' in \code{select}.
#' @param custom_function A character string of a custom processing function
#' that is later parsed (produces expression in a list) and evaluated within
#'  the function \code{preprocess_spc}.
#' The character vector argument of \code{custom_function}
#' needs to contain \code{"spc_raw"}, which is the single data table of spectra
#' that results from binding a list of data.tables (spectra to preprocess)
#' from the spectra list-column specified in \code{column_in}.
#' An example for a value is
#' \code{"prospectr::savitzkyGolay(X = spc_raw, m = 0, p = 3, w = 9)"}.
#' Optional argument. Default is \code{NULL}.
#' @export
preprocess_spc <- function(spc_tbl, select, column_in = "spc_mean",
  custom_function = NULL) {

  # Convert list of spectral data.tables to one data.table
  spc_raw <- data.table::rbindlist(spc_tbl[column_in][[column_in]])

  ## Perform preprocessing =====================================================

  # Use custom function when supplying option ----------------------------------
  if (!is.null(custom_function) & select == "custom") {
    # Create full character string for parsing
    custom_fct <- paste0("custom <- ", custom_function)
    # parse converts the character string into an expression
    # eval evaluates the expression; as a result, custom object is computed
    # and saved within the current workspace
    eval(parse(text = custom_fct))
    ## x <- spc_raw
    ## custom <- eval(substitute(custom_function), envir = parent.frame())
    # -> Error in is.data.frame(X) : object 'x' not found
  }
  # -> returns error:
  # custom_function = prospectr::savitzkyGolay(X = x, m = 0, p = 3, w = 9)
  # Error in is.data.frame(X) : object 'x' not found
  # -> Maybe solution: http://stackoverflow.com/questions/30563745/non-standard-evaluation-from-another-function-in-r

  # Savitzky-Golay preprocessing
  # use different derivatives and window sizes ---------------------------------

  # Zero order Savitzky-Golay (no derivative) -> only smoothing
  if (select == "sg_0_w9") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)}
  # First derivative Savitzky-Golay
  if (select == "sg_1_w5") {
    sg_1_w5 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 5)}
  if(select == "sg_1_w9") {
    sg_1_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 9)}
  if (select == "sg_1_w11") {
    sg_1_w11 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 11)}
  if (select == "sg_1_w13") {
    sg_1_w13 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 13)}
  if (select == "sg_1_w15") {
    sg_1_w15 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 15)}
  if(select == "sg_1_w17") {
    sg_1_w17 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 17)}
  if (select == "sg_1_w19") {
    sg_1_w19 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 19)}
  # Implement window size of 21, corresponds to ICRAF standard;
  # see e.g. Terhoeven-Urselmans et al. (2010)
  if (select == "sg_1_w21") {
    sg_1_w21 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 21)}
  if (select == "sg_1_w23") {
    sg_1_w23 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 23)}
  if(select == "sg_1_w25") {
    sg_1_w25 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 25)}
  if (select == "sg_1_w27") {
    sg_1_w27 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 27)}
  if (select == "sg_1_w35") {
    sg_1_w35 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 35)}
  if (select == "sg_1_w41") {
    sg_1_w41 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 41)}
  if (select == "sg_1_w51") {
    sg_1_w51 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 51)}
  # Second derivative Savitzky-Golay
  if (select == "sg_2_w11") {
    sg_2_w11 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 2, p = 3, w = 11)}
  if (select == "sg_2_w21") {
    sg_2_w21 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 2, p = 3, w = 21)}
  # Savitzky-Golay (order 0) smoothing and derivative with a window size of
  # 21 points
  if (select == "sg_0_1_w21") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_1_w21 <- prospectr::savitzkyGolay(X = sg_0_w9,
      m = 1, p = 3, w = 21)}
  # Savitzky-Golay second derivative
  if (select == "sg_2_w5") {
    sg_2_w5 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 2, p = 3, w = 5)}
  if (select == "sg_2_w11") {
    sg_2_w11 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 2, p = 3, w = 11)}

  # Standard normal variate (SNV) ----------------------------------------------

  # Calculate standard normal variate (SNV) after Savitzky-Golay smoothing
  if (select == "sg_0_snv") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)}
  if (select == "sg_1_snv") {
    sg_1_snv <- prospectr::standardNormalVariate(sg_1_w5)}
  # Standard normal variate (SNV) and first gap-segment derivative
  if (select == "snv_gsd_m1_w11_s1") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)
    snv_gsd_m1_w11_s1 <- prospectr::gapDer(X = sg_0_snv, m = 1, w = 11, s = 1)}
  if (select == "snv_gsd_m1_w21_s5") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)
    snv_gsd_m1_w21_s5 <- prospectr::gapDer(X = sg_0_snv, m = 1, w = 21, s = 5)}
  if (select == "snv_gsd_m1_w31_s1") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)
    snv_gsd_m1_w31_s1 <- prospectr::gapDer(X = sg_0_snv, m = 1, w = 31, s = 5)}
  if (select == "snv_gsd_m1_w31_s5") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)
    snv_gsd_m1_w31_s5 <- prospectr::gapDer(X = sg_0_snv, m = 1, w = 31, s = 5)}
  # Standard normal variate (SNV) and second gap-segment derivative
  if (select == "snv_gsd_m2_w5_s1") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)
    snv_gsd_m2_w5_s1 <- prospectr::gapDer(X = sg_0_snv, m = 2, w = 5, s = 1)}
  if (select == "snv_gsd_m2_w21_s1") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)
    snv_gsd_m2_w21_s1 <- prospectr::gapDer(X = sg_0_snv, m = 2, w = 21, s = 1)}
  if (select == "snv_gsd_m2_w31_s1") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)
    snv_gsd_m2_w31_s1 <- prospectr::gapDer(X = sg_0_snv, m = 2, w = 31, s = 5)}
  if (select == "snv_gsd_m2_w31_s5") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)
    snv_gsd_m2_w31_s5 <- prospectr::gapDer(X = sg_0_snv, m = 2, w = 31, s = 1)}
  if (select == "snv_gsd_m2_w51_s1") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)
    snv_gsd_m2_w51_s1 <- prospectr::gapDer(X = sg_0_snv, m = 2, w = 51, s = 1)}
  if (select == "snv_gsd_m2_w51_s5") {
    sg_0_w9 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 0, p = 3, w = 9)
    sg_0_snv <- prospectr::standardNormalVariate(sg_0_w9)
    snv_gsd_m2_w51_s5 <- prospectr::gapDer(X = sg_0_snv, m = 2, w = 51, s = 5)}
  # 1rst Gap-segement derivative
  if (select == "gsd_m1_w5_s4") {
    gsd_m1_w5_s4 <- prospectr::gapDer(X = spc_raw, m = 1, w = 5, s = 4)}
  if (select == "gsd_m1_w11_s5") {
    gsd_m1_w11_s5 <- prospectr::gapDer(X = spc_raw, m = 1, w = 11, s = 5)}
  if (select == "gsd_m1_w11_s21") {
    gsd_m1_w11_s21 <- prospectr::gapDer(X = spc_raw, m = 1, w = 11, s = 21)}
  if (select == "gsd_m1_w21_s1") {
    gsd_m1_w21_s1 <- prospectr::gapDer(X = spc_raw, m = 1, w = 21, s = 1)}
  if (select == "gsd_m1_w21_s21") {
    gsd_m1_w21_s21 <- prospectr::gapDer(X = spc_raw, m = 1, w = 21, s = 21)}
  if (select == "gsd_m1_w35_s21") {
    gsd_m1_w35_s21 <- prospectr::gapDer(X = spc_raw, m = 1, w = 35, s = 21)}
  if (select == "gsd_m1_w5_s21") {
    gsd_m1_w5_s21 <- prospectr::gapDer(X = spc_raw, m = 1, w = 5, s = 21)}
  # 2nd Gap-segment derivative
  if (select == "gsd_m2_w21_s21") {
    gsd_m2_w21_s21 <- prospectr::gapDer(X = spc_raw, m = 2, w = 21, s = 21)}
  # 4th Gap-segment derivative
  if (select == "gsd_m4_w21_s21") {
    gsd_m4_w21_s21 <- prospectr::gapDer(X = spc_raw, m = 4, w = 21, s = 21)}

  # Savitzky-Golay combined with multiple scatter correction (MSC --------------
  # Savitzky-Golay with 3rd order polynomial, a window size of 21
  # and first derivative + MSC
  if (select == "sg_1_w21_msc") {
    sg_1_w21 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 1, p = 3, w = 21)
    # Use msc function from the pls package; use column means of X as reference
    # spectrum
    sg_1_w21_msc <- pls::msc(X = sg_1_w21, reference = NULL)
  }
  # Savitzky-Golay combined with multiple scatter correction (MSC --------------
  # Savitzky-Golay with 4th order polynomial, a window size of 21
  # and second derivative + MSC
  if (select == "sg_2_w21_msc") {
    sg_2_w21 <- prospectr::savitzkyGolay(X = spc_raw,
      m = 2, p = 4, w = 21)
    # Use msc function from the pls package; use column means of X as reference
    # spectrum
    sg_2_w21_msc <- pls::msc(X = sg_2_w21, reference = NULL)
  }

  # Continuum-removal ----------------------------------------------------------
  if (select == "cr") {
    cr <- prospectr::continuumRemoval(X = spc_raw,
      wav = as.numeric(colnames(spc_raw)), type = "A")}
  if (select == "cr_refl")
    cr_refl <- prospectr::continuumRemoval(X = spc_raw,
      wav = as.numeric(colnames(spc_raw)), type = "R")}

  # Select final preprocessing based on selection argument and
  # save matrix in data.table
  pre <- select
  spc_pre <- data.table::as.data.table(get(pre))

  # Convert preprocessed spectra in data.table to list of data.table spectra
  spc_pre_list <- split(spc_pre, seq(nrow(spc_pre)))
  # Convert x-values of preprocessed spectra in list of vectors
  # prospectr only hands over new xunits in matrix colnames of type character
  xvalues_pre_list <- lapply(spc_pre_list,
    function(x) as.numeric(colnames(x)))

  # Add list of preprocessed spectra and correspoding wavenumbers to tibble
  spc_tbl <- tibble::add_column(spc_tbl,
    spc_pre = spc_pre_list, xvalues_pre = xvalues_pre_list)
}
