#' @title Gather spectra from list of spectral data into a tibble object
#' @description Gather specta and spectrometer metadata from list into a tibble.
#' Spectra and wavenumbers are stored in list-columns. A tibble is an extended
#' data frame and each spectrum can contain complex data and metadata that
#' are in a rectangular data structure. List-columns is tidy
#' data structure concept that can be combined with functional programming
#' frameworks provided by e.g. the purrr package.
#' @param data list with file name elements that contain spectra and metadata
#' after reading binary OPUS files with \code{simplerspec::read_opus_univ()}
#' @usage gather_spc(data)
#' @return Spectral data and metadata in object of class tible
#' @export
gather_spc <- function(data) {

  ## Extract data from list by dplyr::map variants -----------------------------

  # Extract original spectral matrix for all scans
  # First, try to map all spectra; if some of the list elements are NULL,
  # remove all spectra and metadata for Bruker files that were
  # not successfully read (NULL in final spectra in sublist "spc")
  map_spc <- purrr::map(data, "spc")
  which_NULL <- which(sapply(map_spc, is.null))
  if (length(which_NULL > 0)) {
    message(paste0("Sample spectra originating from the following",
      " OPUS files could not gathered because extracted spectra are NULL: <",
      paste(names(which_NULL), collapse = ";"), ">. ",
      "The following list positions have therefore been omitted from
      <data> when gathering spectra from list into tibble: ",
      paste(which_NULL, collapse = "; ")), ".")
    data <- data[- which_NULL]
    map_spc <- map_spc[- which_NULL]
  }
  # Extract metadata list elements and combine into data.frame
  map_metadata_df <- purrr::map_df(data, "metadata")
  map_metadata <- purrr::map(data, "metadata")
  # Extract rownames of spectra; remove names of rownames vector
  rownames_spc <- unname(unlist(lapply(map_spc, rownames)))
  # Extract wavenumbers
  map_wavenumbers <- purrr::map(data, "wavenumbers")

  ## Create list-column tibble
  data_tibble <- tibble::as_tibble(
    map_metadata_df[c("unique_id", "file_id", "sample_id")]
  )

  ## Add spectra and wavenumbers
  tibble::add_column(.data = data_tibble,
    # raw spectra
    spc = map_spc,
    wavenumbers = map_wavenumbers,
    metadata = map_metadata
  )
}
