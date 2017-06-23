#' @title Gather spectra from list of spectral data into a tibble object
#' @description Gather specta and spectrometer metadata from list into a tibble.
#' Spectra and wavenumbers are stored in list-columns. A tibble is an extended
#' data frame and each spectrum can contain complex data and metadata that
#' are in a rectangular data structure. List-columns is tidy
#' data structure concept that can be combined with functional programming
#' frameworks provided by e.g. the purrr package.
#' @param data list with file name elements that contain spectra and metadata
#' after reading binary OPUS files with \code{simplerspec::read_opus()}
#' @usage gather_spc(data)
#' @return Spectral data and metadata in object of class tible
#' @export
gather_spc <- function(data) {

  ## Extract data from list
  # Extract metadata list elements and combine into data.frame
  map_metadata_df <- purrr::map_df(data, "metadata")
  map_metadata <- purrr::map(data, "metadata")
  # Extract original spectral matrix for all scans
  map_spc <- purrr::map(data, "spc")
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
