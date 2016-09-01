## Function 1: Read spectra in text form ========================
#' @title Read an OPUS text file and extract metadata
#' @description
#' Read single text file acquired with
#' an Bruker Vertex FTIR Instrument
#' (as exported from OPUS software) and extract sample metadata
#' provided in the filename
#' @usage
#' read_spectra(path)
#' @param path character of the directory
#' where the spectral text files are stored
#' @return
#' List that contains the following elements:
#' \itemize{
#'  \item \code{MIR}: data.frame that contains all the spectra.
#'  The columns of \code{MIR} contain absorbance values at
#'  different wavenumber in the MIR range. The wavenumbers
#'  rounded to 0.1 are given as column names. The original file
#'  names are stored as row names. One line in the data frame
#'  \code{MIR} contains one replicate scan of a sample.
#'  \item \code{data_rep}: data.frame that constists of sample
#'  metadata that was extracted from the file name of
#'  individual spectral files. The first vector \code{ID}
#'  contains the spectral file name without the repetition number
#'  supplied as \code{.<number>} in the file name.
#'  Letters 1 to 2 of the spectral
#'  file name are used for the country abbreviation, stored
#'  as in the \code{} vector \code{data_rep} . Letters
#'  4 to 5 of the file name are used for the landscape (site)
#'  abbreviation.
#' }
#' @note: This function is derived from  a re-factored and
#' simplified version of the \code{read.opus} function from the
#' \sQuote{soil.spec} package for reading OPUS VERTEX files
#' The function should also work for other OPUS files (eg alpha),
#' see \code{read.opus}. The function readOPUS() was
#' written by Antoine Stevens.
#' @export
read_spectra <- function(path){
  # Needs various utilities for spectral processing that Antoine
  # put on github
  ID <- NULL
  # Load the MIR data exported from OPUS to txt files
  # List files in the directory
  lf <- list.files(path, full.names = TRUE)
  # Read files into R with readOPUS()
  # (comes from the github file)
  MIR <- readOPUS(lf, in_format = "txt")
  # Wavenumber, from 3996.4 to 599.8 cm-1
  colnames(MIR) <- round(as.numeric(colnames(MIR)), 1)
  # Remove the txt extension
  rownames(MIR) <- sub("\\.txt", "", row.names(MIR))
  # Prepare a dataset with ID,
  # Extract country with substring
  # and repetition
  # (ID = character before the dot; rep = number after the dot)
  data_rep <- data.frame(ID = sub("(.+)\\.[[:digit:]]+$", "\\1",
    row.names(MIR)),
    rep = as.numeric(sub(".+\\.([[:digit:]])+$", "\\1",
      row.names(MIR)))
  )
  data_rep <- cbind(data_rep,
    country = substring(data_rep$ID, first = 1, last = 2),
    site = substring(data_rep$ID, first = 4, last = 5)
  )
  list_spectra <- list(
    MIR = MIR,
    data_rep = data_rep
  )
  return(list_spectra)
}
