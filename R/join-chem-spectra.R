## Join chemical and spectral data ==============================
#' @title Join chemical and spectral data frames
#' @description Combines spectral data (data.frame) and chemical
#' data (data.frame).
#' @param dat_chem data.frame that contains chemical values of
#' the sample
#' @param dat_spec List that contains spectral data
#' @return List: xxx
#' @param by character of column name that defines sample_ID
#' @export
join_chem_spec <- function(
  dat_chem, dat_spec, by = "sample_ID") {
  # Alternative when "no visible binding for global variable":
  data_meta <- MIR <- MIR_pre <- ori <- MIR_mean <- NULL
  # http://stackoverflow.com/questions/23475309/in-r-is-it-possible-to-suppress-note-no-visible-binding-for-global-variable
  # Replace sample_ID by ID
  if(!is.data.frame(dat_chem)) {
    stop(dat_chem, "needs to be a data.frame", call. = FALSE)
  } else {
  colnames(dat_chem)[colnames(dat_chem) == by] <- "ID"
  dat_chem$ID <- as.factor(dat_chem$ID)
  # Select only chemical data that have no outlier spectra
  dat_chem <- dat_chem[dat_spec$data_meta$ID, ]
  ID <- as.factor(dat_spec$data_meta$ID)
  # Join ref analyses
  MIRdata <- data.frame(ID = ID)
  MIRdata$MIR <- dat_spec$MIR_pre
  MIRdata$ori <- dat_spec$MIR_mean
  # Joining by ID, type = "inner"
  MIRdata_chem <- plyr::join(dat_chem, MIRdata, type = "inner")
  # before dplyr::inner_join(dat_chem, MIRdata)
  MIRdata_chem
  }
}

## Join spectra and chemical tibbles ===========================================

#' @title Join spectra data and chemical data tibbles
#' @description Combines spectral data (tibble class) and chemical
#' data (tibble class).
#' @param spc_tbl Tibble that contains spectral data
#' @param chem_tbl Tibble that contains chemical reference values of
#' the samples
#' @param by character of column name that defines sample_ID
#' @return Tibble joined by sample_id
#' @export
join_spc_chem <- function(spc_tbl, chem_tbl, by = "sample_id") {
  if(!is.data.frame(spc_tbl)) {
    stop(dat_chem, "needs to be a Tibble", call. = FALSE)
  } else {
    # Rename column sample_ID to sample_id if sample_ID exists
    if("sample_ID" %in% colnames(chem_tbl)) {
    chem_tbl <- dplyr::rename(chem_tbl, sample_id = sample_ID)
    }
    spc_tbl <- dplyr::inner_join(spc_tbl, chem_tbl)
  }
}
