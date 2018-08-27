#' @title Gather spectra of different types from nested list of spectral data
#' into a tibble object
#' @description Gather specta, spectrometer metadata, measurment id's,
#' x unit values from list into a tibble, one row per spectrum replicate
#' measurement. Spectra, x values and metadata are stored in list-columns.
#' A tibble is an extende data frame and each spectrum can contain complex data
#' and metadata that are in a rectangular data structure. List-columns is tidy
#' data structure concept that can be combined with functional programming
#' frameworks provided by e.g. the purrr package.
#' @param data list with file name (`file_id`) elements that contain spectum
#' types, corresponding x unit type values and metadata; The `data` list
#' is commonly the output after reading binary OPUS files with
#' \code{simplerspec::read_opus_univ()}.
#' @param spc_types Character vector with the spectum types to be extracted
#' and gathered into list columns. The spectrum type names need to match exactly
#' the standardised spectrum type names, which are also names of the second list
#' hierarchy elements of `data`. These values are allowed: `"spc"`: final
#' spectrum after atmospheric compensation if performed, otherwise spectrum
#' after ratioing sample and reference single channel reflectance -- referred to
#' as `AB` in Bruker OPUS software; `"spc_nocomp"`: Spectrum prior to
#' atmospheric compensation if performed; `"sc_sm"`: Single channel reflectance
#' spectrum of the sample; `"sc_rf"`: Single channel reflectance spectrum of
#' the reference (background spectrum); `"ig_sm"`: Interferogram of the sample;
#' `"ig_rf"`: Interferogram of the reference. Default is only extracting final
#' spectra, `spc_types = c("spc")`
#' @usage gather_spc(data = data, spc_types = c("spc", "spc_nocomp", "sc_sm",
#' "sc_rf", "ig_sm", "ig_rf"))
#' @return Spectra related spectral data transformed into (list-)columns
#' of object of class tibble.
#' @importFrom rlang set_names
#' @export
gather_spc <- function(data,
                       spc_types = c("spc")) {
  spc_types <- map(spc_types, rlang::sym)
  spc_types_chr <- map_chr(spc_types, rlang::quo_name)

  # Duplicate original data to refer to full `data` list in messages
  data_origin <- data

  ## Make sure first level `data` elements are named

  if (any(names(data) %in% "")) {
    which_missing <- which(names(data) %in% "")
    idx_nm_missing <- rlang::set_names(names(data)[which_missing],
      which_missing)
    message(paste0(length(which_missing), " `data` ",
      ifelse(length(which_missing) > 1, "elements", "element"),
      ifelse(length(which_missing) > 1, " have", " has"),
      " missing `file_id` names. ",
      ifelse(length(which_missing) > 1, "These", "This"), "`data` list index",
      ifelse(length(which_missing) > 1, " positions\n", " position\n"),
      ifelse(length(which_missing) > 1, "were", "was"),
      " therefore removed:\n\n",
      idx_nm_missing %>%
        purrr::imap_chr(~ paste0(.y, " : ", .x)) %>% paste(collapse = "  "),
      "\n\n"))

    data <- data[-which_missing]
  }

  ## Check whether any file names (`file_id`'s) are duplicated;
  ## remove and raise warning in case duplicated

  is_duplicated <- duplicated(names(data))
  idx_nm_duplicated <- rlang::set_names(names(data)[is_duplicated],
    which(is_duplicated))
  if (any(is_duplicated)) {
    message(paste0(sum(is_duplicated), " `file_id` named ",
      ifelse(sum(is_duplicated) > 1, "elements", "element"), " of `data` ",
      ifelse(sum(is_duplicated) > 1, "are", "is"),
      " non-unique/duplicated\nand ",
      ifelse(sum(is_duplicated) > 1, "were", "was"), " therefore removed:\n\n",
      idx_nm_duplicated %>%
        purrr::imap_chr(~ paste0(.y, " : ", .x)) %>% paste(collapse = "  "),
      "\n\n"))

    data <- data[-is_duplicated]
  }

  ## Check whether `spc_types` are allowed and present at the second level of
  ## the nested list `data`

  spc_types_allowed <- c("spc", "spc_nocomp", "sc_sm", "sc_rf", "ig_sm",
    "ig_rf")
  spc_types_present <- spc_types_chr %in% spc_types_allowed

  if (any(!spc_types_present)) {
    options(useFancyQuotes = FALSE)
    stop("The following spectrum types specified in `spc_types` are",
      " not supported:\n\n",
      paste(dQuote(spc_types_chr[!spc_types_present]), collapse = ", ")) }

  ## Collect names of spectrum types nested at second `file_id` `data` list
  ## level to check whether all specified `spc_types` are present in `data`;
  ## access sub-element names at that list level

  data_types_byfile <- purrr::modify_depth(data, .depth = 1,
    ~ purrr::imap_chr(.x, ~ .y))
  spc_types_byfile_in <- map(data_types_byfile, ~ spc_types_chr %in% .x)

  if (any(map_lgl(spc_types_byfile_in, ~ any(!.x)))) {
    which_rm <- map_lgl(spc_types_byfile_in, ~ any(!.x)) %>% which()
    rm_origin <- names(data_origin[names(which_rm)])
    rm_origin_lgl <- names(data_origin) %in% rm_origin
    which_rm_origin <- rlang::set_names(rm_origin, which(rm_origin_lgl))

    message(paste0("Spectrum types (second list level names) specified in",
      " `spc_types` were not found\nwithin all first level elements ",
      " of list `data` or are NULL (spectra data by `file_id`).\n",
      "Therefore, all data elements for the corresponding",
      " list indices and `file_id`'s\nwere removed from `data`:\n\n",
      which_rm_origin %>%
        purrr::imap_chr(~ paste0(.y, " : ", .x)) %>% paste(collapse = "  "),
      "\n\n"))

    data <- data[- which_rm]
  }

  ## Extract the spectra of different types by plucking all types specified in
  ## `spc_types_chr` argment into separate first list level elements

  spc_mapped <- map(rlang::set_names(spc_types_chr),
    ~ pluck_depth(data = data, .depth = 1, .string = .x))

  ## Match and extract the values for the x unit types

  spc_type_x_unit <- list(
    "spc" = c("wavenumbers", "wavelengths", "x_values"),
    "spc_nocomp" = c("wavenumbers", "wavelengths", "x_values"),
    "sc_sm" = c("wavenumbers_sc_sm", "wavelengths_sc_sm", "x_values_sc_sm"),
    "sc_rf" = c("wavenumbers_sc_rf", "wavelengths_sc_rf", "x_values_sc_rf")
  )

  lcols_x_values_matching <- purrr::flatten_chr(spc_type_x_unit[spc_types_chr])
  x_values_matching <- lcols_x_values_matching[lcols_x_values_matching
    %in% unique(purrr::flatten_chr(data_types_byfile))]

  # `data` has eventually been updated (removed `data` elements)
  data_types_byfile <- purrr::modify_depth(data, .depth = 1,
    ~ purrr::imap_chr(.x, ~ .y))
  x_value_types_byfile_in <- map(data_types_byfile, ~ x_values_matching %in% .x)

  if (any(map_lgl(x_value_types_byfile_in, ~ any(!.x)))) {
    which_rm <- map_lgl(x_value_types_byfile_in, ~ any(!.x)) %>% which()
    rm_origin <- names(data_origin[names(which_rm)])
    rm_origin_lgl <- names(data_origin) %in% rm_origin
    which_rm_origin <- rlang::set_names(rm_origin, which(rm_origin_lgl))

    message(paste0("These x unit types (second list level names)",
      " corresponding to \nspecified `spc_types` spectrum types were not",
      " found within\nall first level elements or are NULL",
      " (spectra data by `file_id`) of list `data`: \n\n",
      paste(dQuote(x_values_matching), collapse = ", ")), "\n\n",
      "Therefore, all data elements for the corresponding",
      " list indices\nand `file_id`'s were removed from `data`:\n\n",
      which_rm_origin %>%
        purrr::imap_chr(~ paste0(.y, " : ", .x)) %>% paste(collapse = "  "),
      "\n\n")

    data <- data[- which_rm]
    spc_mapped <- map(spc_mapped, ~ .x[- which_rm])
  }

  x_values_mapped <- map(set_names(x_values_matching),
    ~ pluck_depth(data = data, .depth = 1, .string = .x))

  ## Combine the mapped spectra and correspoding x unit types, and order
  ## list before returning spectral tibble (`spc_tbl`)

  spc_x_unit_types_order <- c("spc", "wavenumbers", "wavelengths", "x_values",
    "spc_nocomp",
    "sc_sm", "wavenumbers_sc_sm", "wavelengths_sc_sm", "x_values_sc_sm",
    "sc_rf", "wavenumbers_sc_rf", "wavelengths_sc_rf", "x_values_sc_rf",
    "ig_sm", "ig_rf")

  spc_x_values_mapped <- c(spc_mapped, x_values_mapped)
  spc_x_unit_types_matched <- spc_x_unit_types_order[spc_x_unit_types_order %in%
    names(spc_x_values_mapped)]
  spc_x_values_mapped <- spc_x_values_mapped[spc_x_unit_types_matched]

  ## Check if `metadata` elements are present for all list first level elements

  # `data` has eventually been updated (removed `data` elements)
  data_types_byfile <- purrr::modify_depth(data, .depth = 1,
    ~ purrr::imap_chr(.x, ~ .y))
  metadata_byfile_in <- map(data_types_byfile, ~ "metadata" %in% .x)

  if (any(map_lgl(metadata_byfile_in, ~ any(!.x)))) {
    which_rm <- map_lgl(metadata_byfile_in, ~ any(!.x)) %>% which()
    rm_origin <- names(data_origin[names(which_rm)])
    rm_origin_lgl <- names(data_origin) %in% rm_origin
    which_rm_origin <- rlang::set_names(rm_origin, which(rm_origin_lgl))

    message(paste0("`metadata` (second list level names) was not",
      " found within\n all first level elements (`metadata` by `file_id`)",
      " of list `data` or are NULL.\n\n",
    "Therefore, all data elements for the corresponding",
    " list indices\nand `file_id`'s were removed from `data`:\n\n",
    which_rm_origin %>%
      purrr::imap_chr(~ paste0(.y, " : ", .x)) %>% paste(collapse = "  "),
        "\n\n"))

    data <- data[- which_rm]
    spc_x_values_mapped <- map(spc_x_values_mapped, ~ .x[- which_rm])
  }

  ## Extract metadata list elements and combine into tibble

  metadata_mapped <- map(set_names(list("metadata")),
    ~ pluck_depth(data = data, .depth = 1, .string = .x))
  metadata_df <- map_df(data, "metadata")
  id_tbl <- tibble::as_tibble(
    metadata_df[c("unique_id", "file_id", "sample_id")]
  )

  ## Column bind all tibbles of metadata id's, metadata, spectra, and
  ## x unit type values; return on combined spectral tibble (`spc_tbl`)

  spc_tbl <- tibble::as_tibble(spc_x_values_mapped)
  metadata_tbl <- tibble::as_tibble(metadata_mapped)
  dplyr::bind_cols(list(id_tbl, metadata_tbl, spc_tbl))
}


# Helpers ----------------------------------------------------------------------

# Pluck out single elements from a list hierarchy level (`.depth`)
# specified in `.string` based on element names; All elements with the specified
# name will be transposed into a single list named by the element name.
# `pluck_depth` can be combined with `purrr::map` to extract elements of
# different names at all list nodes at the specified level, using e.g.
# map(rlang::set_names(c("spc", "spc_nocomp")),
#   ~ pluck_depth(data = data, .depth = 1, .string = .x))

pluck_depth <- function(data, .depth = 1, .string = NULL) {
  purrr::modify_depth(data, .depth = .depth, ~ purrr::pluck(.x, .string))
}
