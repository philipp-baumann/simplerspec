#' @title Plot tibble spectra
#' @description Plot spectra from tibble spectra objects.
#' @param spc_tbl
#' @param spc_tbl_2
#' @param y
#' @param by
#' @param graph_id_1
#' @param graph_id_2
#' @param graph_id_1_col
#' @param graph_id_2_col
#' @param xlab
#' @param ylab
#' @param slice
#' @param alpha
#' @param legend
#' @usage plot_spc(spc_tbl, spc_tbl_2,
#' y = "spc_pre", by = "sample_id", graph_id_1 = "Calibration",
#' graph_id_2 = "Validation", slice = TRUE, alpha = 1)
#' @export
plot_spc <- function(spc_tbl, spc_tbl_2 = NULL, x = NULL,
                     x_unit = "wavenumber",
                     y = "spc", by = "unique_id",
                     graph_id_1 = "graph_1", graph_id_2 = NULL,
                     graph_id_1_col = "black", graph_id_2_col = "red",
                     xlab = expression(paste("Wavenumber [", cm^-1, "]")),
                     ylab = "Absorbance", slice = TRUE, alpha = 0.2,
                     legend = TRUE) {

  # Fix `R CMD check NOTE`: "no visible binding for global variable ‘...‘"
  graph_id <- id <- variable <- value <- NULL

  # (0) Slice spectra tibble to remove triplicate spectra (reps)
  # only sample_id level
  if(slice == TRUE) {
    spc_tbl <- dplyr::slice(spc_tbl)
  }
  # (1) Gather spectra into one data.table
  if(!is.null(spc_tbl_2)) {
    if(y == "spc") {
      # raw spectra are not yet data.tables and extraction is done alternatively
      # via do.call(rbind, list) -> a little bit slower
      dt_1 <- data.table::data.table(do.call(rbind, spc_tbl[, y][[y]]))
      dt_2 <- data.table::data.table(do.call(rbind, spc_tbl_2[, y][[y]]))
    } else {
      dt_1 <- data.table::rbindlist(spc_tbl[, y][[y]])
      dt_2 <- data.table::rbindlist(spc_tbl_2[, y][[y]])
  }

  } else {
    if(y == "spc") {
      # raw spectra are not yet data.tables and extraction is done alternatively
      # via do.call(rbind, list) -> a little bit slower
      dt_1 <- data.table::data.table(do.call(rbind, spc_tbl[, y][[y]]))
    } else {
      dt_1 <- data.table::rbindlist(spc_tbl[, y][[y]])
  }
  }

  # (2) Extract ID variable and append it to the data.table
  id_1 <- spc_tbl[, by][[by]]
  # Add a graph identity column to distiguish graphical layers for spectra
  # tibble comparisons
  graph_id_1 <- as.factor(rep(graph_id_1, nrow(spc_tbl)))
  dt_1[, graph_id:=graph_id_1]
  dt_1[, id:=id_1]

  # Only if spc_tbl_2 exists
  if(!is.null(spc_tbl_2)) {
    # (2) Extract ID variable and append it to the data.table
    id_2 <- spc_tbl_2[, by][[by]]
    # Add a graph identity column to distiguish graphical layers for spectra
    # tibble comparisons
    graph_id_2 <- as.factor(rep(graph_id_2, nrow(spc_tbl_2)))
    dt_2[, graph_id:=graph_id_2]
    dt_2[, id:=id_2]
    dt_list <- list(
      dt_1 = dt_1,
      dt_2 = dt_2
    )
    dt <- data.table::rbindlist(dt_list)
  } else {
    dt <- dt_1
  }
  # (3) Convert data.table from wide to long form
  dt_long <- data.table::melt(
    dt, measure=names(dt)[!names(dt) %in% c("id", "graph_id")]
  )
  # Convert variable column from factor to numeric
  dt_long[, variable := as.numeric(as.character(variable))]
  # (4) Plot spectra
  # Define nice breaks for x axis
  brk  <- pretty(as.numeric(names(dt)[!names(dt) %in% c("id", "graph_id")]), n = 10)
  p <- ggplot2::ggplot(dt_long, ggplot2::aes(variable, value)) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_reverse(breaks = brk) +
    # Bring graph_id_2 spectra to front
    # http://stackoverflow.com/questions/21120088/ggplot2-bring-one-line-to-the-front-but-save-the-colors
    ggplot2::geom_line(ggplot2::aes(colour = graph_id, group = id),
      alpha = alpha, size = 0.2) +
    # scale_color_manual(values = rep("black", nrow(dt)))
    ggplot2::scale_color_manual(values = c(graph_id_1_col, graph_id_2_col))

  if("wavelengths_rs" %in% names(spc_tbl) && x_unit == "wavelength") {
    p <- p +
      ggplot2::scale_x_continuous() +
      xlab("Wavelength [nm]") +
      ylab("Reflectance")
  }

  if(legend == FALSE) {
    p <- p +
    # Remove legend
    guides(colour = FALSE)
  }
  p
}
