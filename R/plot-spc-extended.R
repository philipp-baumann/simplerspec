################################################################################
## Helper functions to gather spectra, corresponding x-value vectors,
## metadata and measure columns (e.g. chemical reference data) from tibble
## list-columns into a single data.table or a list of data.tables conntaining
## long form data directly to be used for customized ggplot2 plotting functions
################################################################################

# bind a list-column in a tibble to a list of data.tables ----------------------

#' @title Bind list-columns within a tibble into a list of data.tables
#' @description Bind one to many list-columns in spectral tibble into a list
#' of data.tables.
#' @param spc_tbl Spectral data in a tibble data frame (classes "tibble_df",
#' "tbl" and "data.frame").
#' @param lcols Character vector of column names of list-columns to be bound
#' into a list of data.tables
#' @param spc_id Character vector denoting column name for a unique spectrum ID.
#' Default is \code{"unique_id"}.
#' @param group_id Character vector denoting column name for the spectrum group
#' ID. Default is \code{"sample_id"}. The group ID can later be used for
#' plotting spectra by group (e.g. by using different colors or panels).
#' @return A list of data.tables. Elements contain data from list-columns
#' specified in \code{lcols} argument as data.tables. All data.tables contain in
#' addition \code{spc_id} and \code{group_id} columns.
#' @export
bind_lcols_dts <- function(spc_tbl, lcols,
                           spc_id = "unique_id",
                           group_id = "sample_id") {

  # todo: add warning for lcols not present in spc_tbl
  which_bind <- colnames(spc_tbl) %in% lcols
  lcols_to_bind <- colnames(spc_tbl)[which_bind]
  names(lcols_to_bind) <- lcols_to_bind
  dts <- purrr::map(lcols_to_bind,
    function(y) {
      if (is.list(spc_tbl[, y][[y]])) {
        # todo: Test if number of columns is equal in each data.frame or matrix
        # of the list-(column); if not, return a comprehensible error
        data.table::data.table(do.call(rbind, spc_tbl[, y][[y]]))
      } else if (is.atomic((spc_tbl[, y][[y]]))) {
        data.table::data.table(spc_tbl[, y])
      }
    }
  )
  # Append IDs to data.tables in list
  spc_id <- spc_tbl[, spc_id][[spc_id]]
  lcol_types <- purrr::imap(dts, ~ rep(.y, nrow(spc_tbl)))
  group_id <- as.character(spc_tbl[, group_id][[group_id]])

  # Return list of data.tables
  purrr::imap(dts, function(dt, nm) {
    dt[, `:=` (spc_id = spc_id, group_id = group_id)]
    dt[, `:=` (lcol_type = lcol_types[[nm]])]}
  )
}


# Convert list of wide form data.tables into long form -------------------------

dts_to_long <- function(spc_tbl, lcols,
                        spc_id = "unique_id",
                        group_id = "sample_id",
                        variable_name = "variable",
                        value_name = "value") {

  dts <- bind_lcols_dts(spc_tbl = spc_tbl, lcols = lcols,
    spc_id = spc_id, group_id = group_id)
  # Convert list of data.tables into long form
  dts_long <- purrr::map(dts, function(x) {
    data.table::melt(
      x,
      id.vars = c("spc_id", "lcol_type", "group_id"),
      variable.factor = FALSE,
      variable.name = variable_name,
      value.name = value_name
    )}
  )
  # Append unique id (idx lg:= index 'long') for long form
  purrr::imap(dts_long,
    function(dt_long, nm) {
      dt_long[, `:=` (id_long = 1:nrow(dt_long))]
    }
  )
}


# Match the spectra list columns and corresponding xunit list columns ----------

match_lcols <- function(spc_tbl, lcols) {

  # Determine to which spectrum types list-columns belong
  lcols_spc_all <- c("spc",  "spc_rs", "spc_mean", "spc_nocomp", "sc_sm",
    "sc_rf", "spc_pre")
  xvalue_lookup <- list(
    "spc" = c("wavenumbers", "wavelengths"),
    "spc_rs" = c("wavenumbers_rs", "wavelengths_rs"),
    "spc_mean" = c("wavenumbers_rs", "wavelengths_rs"),
    "spc_nocomp" = c("wavenumbers", "wavelengths"),
    "sc_sm" = c("wavenumbers_sc_sm"),
    "sc_rf" = c("wavenumbers_sc_rf"),
    "spc_pre" = c("xvalues_pre")
  )
  # Create character vector of spectra type names
  spc_matched <- lcols[lcols %in% lcols_spc_all]
  spc_matched <- spc_matched[order(match(spc_matched, lcols_spc_all))]
  # Create vector of corresponding xunit types in predefined order
  xvalues <- unlist(xvalue_lookup[spc_matched])
  xvalues_matched <- colnames(spc_tbl)[colnames(spc_tbl) %in% xvalues]
  xvalues_matched <- xvalues_matched[order(match(xvalues_matched, xvalues))]
  # Return all matches as list of character vectors for spectra and x-values
  list(
    spc_matched = spc_matched,
    xvalues_matched = xvalues_matched
  )
}


# Create a list of long form data.tables containing spectra and x-values for
# a set of spectral types ------------------------------------------------------

tolist_spc_xvalues <- function(spc_tbl, lcols_spc,
                      spc_id = "unique_id",
                      group_id = "sample_id",
                      variable_name = "variable",
                      value_name = "value") {

  lcols_matched <- match_lcols(spc_tbl = spc_tbl, lcols = lcols_spc)
  # Check if length of matched spectra and xunits is equal and if not,
  # return an error
  spc_types <- lcols_matched[["spc_matched"]]
  xunit_types <- lcols_matched[["xvalues_matched"]]

  # Gather different spectra types into list of data.tables
  spc_dts <- dts_to_long(spc_tbl = spc_tbl,
    lcols = spc_types, spc_id = spc_id, group_id = group_id,
    variable_name = "spc_variable", value_name = "spc_value")
  # Gather corresponding xunit types into list of data.tables
  xvalues_dts <- dts_to_long(spc_tbl = spc_tbl,
    lcols = xunit_types, spc_id = spc_id, group_id = group_id,
    variable_name = "xvalues_variable", value_name = "xvalues_value")

  # Rename lcol_type to spc_type only for spectra data.tables
  spc_dts <- purrr::map(spc_dts,
    function(x) data.table::setnames(x, "lcol_type", "spc_type"))
  # Return data.tables in nested list
  list(
    "spc_dts" = spc_dts,
    "xvalues_dts" = xvalues_dts
  )
}


# Merge data tables of spectra, xunits, metadata and measured variables
# into a single long form data.table -------------------------------------------

#' @title Merge list-columns of spectra, x-axis values, metadata and additional
#' measured variables into a single long form data.table
#' @description Helper function that merges all spectra and related data into
#' a single long form data.table than can subsequently be used for plotting.
#' @param spc_tbl Tibble data frame containing spectra, x-axis values, metadata
#' and eventual measured variables as list-columns.
#' @param lcols_spc Character vector of spectral list-columns to be extracted.
#' Default is \code{c("spc", "spc_pre")} (raw and preprocessed spectra).
#' @param lcol_measure Character vector of length 1 denoting the column name
#' of the measure columns. This argument is optional. Default is \code{NULL},
#' which does not extract an additional measure column.
#' @param spc_id Character vector of column that contains a unique spectral
#' identifier for all spectra. Default is \code{"unique_id"}.
#' @param group_id Character vector of columns that is used assigning spectra
#' into groups. Default is \code{"sample_id"}. The \code{group_id} can be
#' used for later plotting and thereby visually separating spectral groups into
#' using different colors or panels.
#' @return A single data.table containing long form aggregated data of
#' spectra, x-axis values, metadata and an additionally measured variable.
#' @export
merge_dts <- function(spc_tbl,
                      lcols_spc = c("spc", "spc_pre"), lcol_measure = NULL,
                      spc_id = "unique_id",
                      group_id = "sample_id") {

  id_long <- NULL
  spc_xvalues <- tolist_spc_xvalues(spc_tbl = spc_tbl,
    lcols_spc = lcols_spc, spc_id = spc_id, group_id = group_id)
  # Set keys for merging list of data.tables for spectra and xunits
  purrr::imap(
    spc_xvalues,
    function(dts, nm) purrr::map(dts[[nm]],
      function(x) data.table::setkey(x = x, spc_id, id_long, group_id))
  )
  spc_xvalues <- purrr::map2(spc_xvalues[["spc_dts"]],
    spc_xvalues[["xvalues_dts"]], merge)

  # Bind metadata if present, and set keys for merging metadata to spectra
  metadata <- bind_lcols_dts(spc_tbl = spc_tbl,
    lcols = "metadata", spc_id = spc_id, group_id = group_id)
  dts <- list(
    "data" = spc_xvalues,
    "metadata" = rep(metadata, length(spc_xvalues))
  )
  if (length(metadata) == 0) dts$metadata <- NULL

  # Convert a "measure" tibble column (numeric|character) to list of data.tables
  if (!is.null(lcol_measure)) {
    measure <- bind_lcols_dts(spc_tbl = spc_tbl,
      lcols = lcol_measure, spc_id = spc_id, group_id = group_id)
    dts$measure <- rep(measure, length(spc_xvalues))
  }
  # Set keys (common columns), merge metadata with spectral data (list of
  # data tables) and combine into a single data.table that is returned
  purrr::imap(dts,
    function(dt, nm) lapply(dts[[nm]],
      function(x) data.table::setkey(x = x, spc_id, group_id))
  )
  # Merge multiple data.table by common keys
  # https://gist.github.com/reinholdsson/67008ee3e671ff23b568
  data.table::rbindlist(
    lapply(seq_along(dts[[1]]),
      function(i) Reduce(merge, lapply(dts, `[[`, i)))
  )
}


# Wrapper function around merge_dts for list of tibbles to aggregate data for
# plotting ---------------------------------------------------------------------

#' @title Wrapper function around \code{merge_dts()} for list of tibbles to
#' aggregate data for plotting.
#' @description Instead of a single spectral tibble (data frame) multiple
#' spectral tibbles can be merged into a long-form data.table for plotting
#' spectra and related data. For details, see
#' \code{\link{merge_dts}}.
#' @param spc_tbl_l List of spectral tibbles (data frames).
#' @param lcols_spc Character vector of spectral list-columns to be extracted.
#' Default is \code{c("spc", "spc_pre")} (raw and preprocessed spectra).
#' @param lcol_measure Character vector of length 1 denoting the column name
#' of the measure columns. This argument is optional. Default is \code{NULL},
#' which does not extract an additional measure column.
#' @param spc_id Character vector of column that contains a unique spectral
#' identifier for all spectra. Default is \code{"unique_id"}.
#' @param group_id Character vector of columns that is used assigning spectra
#' into groups. Default is \code{"sample_id"}. The \code{group_id} can be
#' used for later plotting and thereby visually separating spectral groups into
#' using different colors or panels.
#' @return A single data.table containing long form aggregated data of
#' spectra, x-axis values, metadata and an additionally measured variable.
#' An additional column called \code{group_id_tbl} is appended. It denotes
#' the name of the spectral tibble supplied with the list \code{spc_tbl_l}.
#' @export
merge_dts_l <- function(spc_tbl_l,
                        lcols_spc = c("spc", "spc_pre"),
                        lcol_measure = NULL,
                        spc_id = "unique_id",
                        group_id = "sample_id") {

  group_id_tbl <- NULL

  dts <- lapply(seq_along(spc_tbl_l),
    function(i) merge_dts(spc_tbl = spc_tbl_l[[i]],
      lcols_spc = lcols_spc, lcol_measure = lcol_measure,
      spc_id = spc_id, group_id = group_id))
  dts <- lapply(seq_along(dts),
    function(i) dts[[i]][, group_id_tbl := names(spc_tbl_l[i])])
  data.table::rbindlist(dts)
}


## Create plotting functions based on complete long data.table =================

# Function that reorders factor column in data.table based on ascending numeric
# order when converted to numeric type
# https://stackoverflow.com/questions/15665535/reorder-factors-numerically-in-a-data-frame
# ------------------------------------------------------------------------------
reorder_factor_num <- function(dt, column = "group_id") {
  group_id <- NULL
  if (!any(is.na(
    suppressWarnings(as.numeric(dt[, column, with = FALSE]))))
  ) {
    dt[, group_id := as.factor(group_id)]
    sorted_labels <- paste(sort(as.numeric(levels(dt$group_id))))
    dt$group_id <- factor(dt$group_id, levels = sorted_labels)
  }
  dt
}


# Custom ggplot2 labeller for spectra types ------------------------------------

relabel_spc_types <- function(lb_sc_sm = "Reflectance sample (<ScSm>)",
                              lb_sc_rf = "Reflectance reference (<ScRf>)",
                              lb_ig_sm = "Interferogram sample (<IgSm>)",
                              lb_ig_rf = "Interferogram reference (<IgRf>)",
                              lb_spc_nocomp = "Abs. before atm. comp.",
                              lb_spc = "Absorbance",
                              lb_spc_rs = "Resampled Abs.",
                              lb_spc_pre = "Preprocessed Abs.") {
  ggplot2::as_labeller(
    x = c(
      "sc_sm" = lb_sc_sm,
      "sc_rf" = lb_sc_rf,
      "ig_sm" = lb_ig_sm,
      "spc_nocomp" = lb_spc_nocomp,
      "spc" = lb_spc,
      "spc_rs" = lb_spc_rs,
      "spc_pre" = lb_spc_pre
    )
  )
}


# Main spectra explorative analysis and diagnostics plotting function ----------

#' @title ggplot2 wrapper for extended spectra plotting
#' @description \code{plot_spc_ext} is a custom plotting function developed
#' within the simplerspec framework. Returns plots based on ggplot2
#' (class "ggplot"). Different spectra types such as raw or preprocessed spectra
#' and groups can be differentiated by different colors or by using panels
#' (so called facets). Additionally, spectra can be colored based on an
#' additional measure variable, e.g. determined by chemical reference analysis.
#' @param spc_tbl Tibble data frame containing spectra, x-axis values, metadata
#' and eventual measured variables as list-columns.
#' @param spc_tbl_l List of spectral tibbles (data frames). Default is
#' \code{NULL} (argument is not used).
#' @param lcols_spc Character vector of spectral list-columns to be extracted.
#' Default is \code{"spc"} (raw spectra).
#' @param lcol_measure Character vector of length 1 denoting the column name
#' of the measure columns. This argument is optional. Default is \code{NULL},
#' which does not extract an additional measure column.
#' @param spc_id Character vector denoting column name for a unique spectrum ID.
#' Default is \code{"unique_id"}.
#' @param group_id Character vector denoting column name for the spectrum group
#' ID. Default is \code{"sample_id"}. The group ID is used for
#' plotting spectra by group (e.g. by using different colors or panels).
#' @param group_id_order Logical that specifies whether the panel names
#' derived from a numeric \code{group_id} column are reordered using ascending
#' numbers. Default is \code{TRUE}.
#' @param group_color Logical defining whether spectra are colored by the column
#' specified by \code{group_id}.
#' @param group_panel Logical defining whether spectra are arranged into panels
#' by groups specified in \code{group_id}. Default is \code{TRUE}.
#' @param group_legend Logical defining whether a legend for the \code{group_id}
#' is plotted. Default is \code{FALSE}.
#' @param ncol Integer vector of length 1. Defines number of columns when
#' plotting panels (facets). Default is \code{NULL} (argument not used).
#' @param relabel_spc Logical defining whether panels are relabeled with custom
#' names for spectra types. Default is TRUE. When \code{TRUE}, arguments
#' from \code{relabel_spc_types} can be passed to \code{plot_spc_ext}
#' (supported via the \code{...} (ellipsis) argument)
#' @param ylab Character vector or vector of type \code{"expression"} created by
#' mathematical expression created by \code{expression}. Custom annotation for
#' y-axis of spectra
#' @param alpha Integer of length 1, from 0 to 1. Defines transparency of
#' spectral lines. Default is \code{0.5} (0 is completely transparent and
#' 1 is no transparency).
#' @param line_width Numeric vector of length 1 specifying the width of the
#' spectral lines. Default is \code{0.2}.
#' @param ... Further arguments to be passed to \code{plot_spc_ext}. Currently,
#' arguments of \code{relabel_spc_types} are supported.
#' @return Object of class \code{"ggplot"} (ggplot2 graph).
#' @export
plot_spc_ext <- function(spc_tbl, spc_tbl_l = NULL,
                          lcols_spc = "spc",
                          lcol_measure = NULL,
                          spc_id = "unique_id",
                          group_id = "sample_id", group_id_order = TRUE,
                          group_color = TRUE,
                          group_panel = TRUE,
                          group_legend = FALSE,
                          ncol = NULL,
                          relabel_spc = TRUE,
                          ylab = "Spectrum value",
                          alpha = 0.5, line_width = 0.2,
                          # Further arguments to be passed to functions called
                          # within this function
                          ...) {

  # Merge spectral data, additional (measurement data) and metadata into a
  # single long-form data.table
  if (!is.null(spc_tbl_l)) {
    dt <- merge_dts_l(spc_tbl_l = spc_tbl_l,
      lcols_spc = lcols_spc, lcol_measure = lcol_measure,
      spc_id = spc_id, group_id = group_id) # see merge_dts_l wrapper function
  } else {
    dt <- merge_dts(spc_tbl = spc_tbl,
      lcols_spc = lcols_spc, lcol_measure = lcol_measure,
      spc_id = spc_id, group_id = group_id)
  }
  # Option to order originally numeric group_id factors by group
  if (group_id_order) {
    dt <- reorder_factor_num(dt = dt, column = "group_id")
  }
  brk <- pretty(dt[["xvalues_value"]], n = 10) # Pretty x-axis breaks
  p <- ggplot2::ggplot(data = dt,
    ggplot2::aes_string(x = "xvalues_value", y = "spc_value"))

  if (group_color == TRUE && is.null(lcol_measure)) {
    p <- p +
      ggplot2::geom_line(ggplot2::aes_string(colour = "group_id",
        group = "spc_id"),
        alpha = alpha, size = line_width)
    if (group_legend == FALSE) {
      p <- p + ggplot2::guides(colour = FALSE)
    }
  } else if (group_color == FALSE && is.null(lcol_measure)) {
    p <- p + ggplot2::geom_line(
      ggplot2::aes_string(group = "spc_id"),
        alpha = alpha, size = line_width)
  }

  if (!is.null(lcol_measure)) {
    p <- p + ggplot2::geom_line(
      ggplot2::aes_string(colour = lcol_measure, group = "spc_id",
        x = "xvalues_value", y = "spc_value"),
        alpha = alpha, size = line_width, inherit.aes = FALSE) +
      ggplot2::scale_colour_distiller(palette = "Spectral")
  }

  # Plot different spectral types and group_id in panels
  if (group_panel && length(lcols_spc) > 1) {
    if (relabel_spc) {
      lbl <- relabel_spc_types(...) # see this function for arguments and values
      p <- p + ggplot2::facet_grid(spc_type ~ group_id, scales = "free",
        labeller = ggplot2::labeller(spc_type = lbl))
    } else {
      p <- p + ggplot2::facet_grid(spc_type ~ group_id, scales = "free")
    }
  }

  if (group_panel && length(lcols_spc) == 1) {
    p <- p + ggplot2::facet_wrap(~ group_id, ncol = ncol, scales = "free")
  }
  # Special case when list of tibbles are supplied
  if (group_panel && !is.null(spc_tbl_l)) {
    p <- ggplot2::ggplot(data = dt,
      ggplot2::aes_string(x = "xvalues_value", y = "spc_value")) +
      ggplot2::geom_line(ggplot2::aes_string(colour = "group_id_tbl",
        group = "spc_id"), alpha = alpha, size = line_width)
    if (relabel_spc == TRUE) {
      lbl <- relabel_spc_types(...)
      p <- p + ggplot2::facet_grid(spc_type ~ group_id, scales = "free",
        labeller = ggplot2::labeller(spc_type = lbl))
    } else if (relabel_spc == FALSE) {
      p <- p + ggplot2::facet_wrap(~ group_id, ncol = ncol, scales = "free")
    }
  }

  p <- p + ggplot2::scale_x_reverse(breaks = brk) +
    ggplot2::xlab(expression(paste("Wavenumber [", cm^-1, "]"))) +
    ggplot2::ylab(ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "right",
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  p
}

