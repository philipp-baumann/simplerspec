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
#' @importFrom tidyselect one_of
#' @export
read_asd <- function(file) {

  # Read fixed with file into a tibble
  asd_tbl <- readr::read_tsv(file = file)
  # Transpose tibble and add Wavelengths as column names
  asd_tbl_t <- tibble::as_tibble(
    t(dplyr::select(asd_tbl, - one_of("Wavelength")))
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

## Simplespec spectra tibble version of ASD reader based on prospectr::readASD
## Reads binary ASD data and converts data into list-columns containing spectral
## data that can be further processed within the simplerspec spectra processing
## framework ===================================================================

#' @title Read ASD binary files and gather spectra and metadata in tibble data
#' frame.
#' @description Read multiple ASD binary files and gather spectra and metadata
#' into a simplerspec spectral tibble (data frame). The resulting spectral
#' tibble is compatible with the simplerspec spectra processing and modeling
#' framework.
#' @param fnames Character vector containing full paths of ASD binary files
#' to be read
#' @return A spectral tibble (data frame) containing the follwing columns:
#' \item{unique_id}{Character vector. Unique identifier containing file name
#' pasted with date and time.}
#' \item{file_id}{Character vector containing file names and exension}
#' \item{sample_id}{Character vector containing files names without extension}
#' \item{metadata}{List-column. List of data frames containing spectral
#' metadata}
#' \item{wavelengths}{List-column. List of wavelengths vectors (numeric).}
#' \item{spc_radiance}{List-column. List of data.tables containing
#' radiance sample spectra.}
#' \item{spc_reference}{List-column. List of data.tables containing
#' reference reflectance spectra.}
#' \item{spc}{List-column. List of data.tables containing final reflectance
#' spectra.}
#' @export
read_asd_bin <- function(fnames) {

  data <- prospectr::readASD(fnames = fnames,
    in_format = "binary", out_format = "list")
  gps <- purrr::map(data, c("header", "GPS"))
  header <- purrr::map(purrr::map(data, "header"),
    function(x) x[- which(names(x) == "GPS")])
  file_id <- purrr::map_chr(data, "name")
  sample_id <- sub("(.+)\\.[[:alpha:]]+$", "\\1", file_id) # remove ".asd"
  datetime <- purrr::map(data, "datetime")
  unique_id <- mapply(function(x, y) paste0(x, "_", y), sample_id, datetime)
  metadata <- purrr::map(header, tibble::as_tibble)
  # Add GPS to metadata
  metadata <- purrr::map2(metadata, gps, dplyr::bind_cols)
  spc_l <- purrr::transpose(
    purrr::map(data, `[`, c("radiance", "reference", "reflectance")))
  wl_l <- purrr::transpose(purrr::map(data, `[`, "wavelength"))
  spc_dt <- purrr::modify_depth(spc_l, 2,
    function(x) data.table::data.table(t(x)))

  spc_tbl <- tibble::tibble(
    unique_id = unique_id,
    file_id = file_id,
    sample_id = sample_id,
    metadata = metadata,
    wavelengths = wl_l[["wavelength"]],
    spc_radiance = spc_dt[["radiance"]],
    spc_reference = spc_dt[["reference"]],
    spc = spc_dt[["reflectance"]]
  )
}

# Helper function to remove the ".asd.xxx" (.xxx for example ".ref" or "")
# extension in id column (e.g. sample_id) strings in tibble with metadata or
# reference analysis data ------------------------------------------------------

#' @importFrom stringr str_replace
#' @importFrom dplyr pull
remove_id_extension <- function(data,
                                id_col = "sample_id",
                                id_new_nm = "sample_id",
                                extension = "\\.asd.*$") {
  id_col <- enquo(id_col)
  id_col_chr <- quo_name(id_col)
  id_col_rm <- rlang::expr(-!!rlang::sym(id_col_chr))
  id_new_nm <- quo_name(enquo(id_new_nm))

  id_new <- stringr::str_replace(string = dplyr::pull(data, !!id_col),
    pattern = extension, replacement = "")

  # Remove old id column and bind new id column to the remaining columns
  rest <- dplyr::select(data, !!id_col_rm)
  dplyr::bind_cols(!!id_new_nm := id_new, rest)
}

