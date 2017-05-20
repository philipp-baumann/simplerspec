#' @title Read tab delimited text (.txt) files from ASD field spectrometer data
#' export into simplerspec spectra tibble.
#' @description For reading ASD file spectrometer data the exported data are
#' expected to be in one .txt tab delimited file. The first row should contain
#' the name 'Wavelength' for the first column and the file names for the
#' remaining columns.
#' @return Spectra data in tibble data frame (class `tbl_df`) that contains
#' columns `sample_id`, `spc` (list of spectral matrices) and 'wavelengths'
#' (list of wavelength vectors).
#' @export
read_asd <- function(file) {

  # Read fixed with file into a tibble
  asd_tbl <- readr::read_tsv(file = "data/asd/170214_all.txt")
  # Transpose tibble and add Wavelengths as column names
  asd_tbl_t <- tibble::as_tibble(
    t(dplyr::select(asd_tbl, - Wavelength))
  )
  colnames(asd_tbl_t) <- asd_tbl[["Wavelength"]]

  # Split matrix by each row into list of matrices
  asd_m <- as.matrix(asd_tbl_t)
  asd_listofv <- split(asd_m, row(asd_m)) # List of numerical vectors
  # Convert list of vectors into list of matrices
  asd_listofm <- lapply(seq_along(asd_listofv),
    function(i) matrix(asd_listofv[[i]], nrow = 1, byrow = FALSE))
  # Assign file names as names for list of matrices
  names(asd_listofm) <- colnames(asd_tbl)[-1] # Remove "Wavelength"

  # Assign columnes for all matrices in list
  asd_listofm <- lapply(asd_listofm,
    function(x) {colnames(x) <- asd_tbl[["Wavelength"]]; x})

  # Create list of wavelengths and assign sample names
  wavelengths_list <- rep(list(asd_tbl[["Wavelength"]]), length(asd_listofm))
  names(wavelengths_list) <- names(asd_listofm)

  # Return a tibble
  spc_tbl <- tibble::data_frame(
    sample_id = names(asd_listofm),
    spc = asd_listofm,
    wavelengths = wavelengths_list
  )

}
