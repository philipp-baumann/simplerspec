#' @title Gather spectra from list of spectral data into a tibble object
#' @description
#' @param data list with file name elements that contain spectra and metadata
#' after reading binary OPUS files with \code{simplerspec::read_opus()}
#' @usage gather_spc(data)
#' @return Spectral data and metadata in object of class tible
#' @export
gather_spc <- function(data) {

  ## Extract data from list
  # Extract metadata list elements and combine into data.frame
  map_metadata_df <- data %>% purrr::map_df("metadata")
  map_metadata <- data %>% purrr::map("metadata")
  # Extract original spectral matrix for all scans
  map_spc <- data %>% purrr::map("spc")
  # Extract rownames of spectra; remove names of rownames vector
  rownames_spc <- unname(unlist(lapply(map_spc, rownames)))
  # Use rbind.fill.matrix function from plyr package to combine rows
  spc <- plyr::rbind.fill.matrix(map_spc)
  # Add rownames to resampled spectra list
  rownames(spc) <- rownames_spc
  # Extract wavenumbers
  map_wavenumbers <- data %>% purrr::map("wavenumbers")

  ## Create list-column tibble
  data_tibble <- as_tibble(
    map_metadata_df[c("unique_id", "scan_id", "sample_id")]
  )

  ## Add spectra and wavenumbers
  data_tibble %>% tibble::add_column(
    # raw spectra
    spc = map_spc,
    wavenumbers = map_wavenumbers,
    metadata = map_metadata
  )
}
