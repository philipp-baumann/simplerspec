#' @title Read ASD fieldspec spectrometer data export into into simplerspec
#' spectra tibble.
#' @description Read tab delimited text (.txt) files exported from ASD field
#' spectrometer into simplerspec spectra tibble.
#' ASD Fieldspec data files are expected in .txt tab delimited file format.
#' The first row should contain
#' the name 'Wavelength' for the first column and the file names for the
#' remaining columns.
#' @param file Tab delmited file from ASD software export where the first
#' column called \code{Wavelength} contais wavelengths in nanometer and the
#' remaining columns are sample spectra referred by an ID name provided in the
#' first row of these columns.
#' @return Spectra data in tibble data frame (class `tbl_df`) that contains
#' columns \code{sample_id} (derived from 2nd and following column names of
#' tab delimited ASD exported text file),
#' \code{spc} (list-column of spectral matrices)
#' and \code{wavelengths} (list-column containing wavelength vectors).
#' @export
read_asd <- function(file) {

  # Read fixed with file into a tibble
  asd_tbl <- readr::read_tsv(file = file)
  # Transpose tibble and add Wavelengths as column names
  asd_tbl_t <- tibble::as_tibble(
    t(dplyr::select(asd_tbl, - rlang::UQS(rlang::syms("Wavelength"))))
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

  # Return spectra as tibble
  tibble::data_frame(
    sample_id = names(asd_listofm),
    spc = asd_listofm,
    wavelengths = wavelengths_list
  )

}
