#' @title Remove outlier spectra
#' @description Remove outlier spectra based on the
#' \code{pcout()} function of the \code{mvoutlier} package.
#' @usage remove_outliers(list_spectra, remove = TRUE)
#' @param list_spectra List that contains averaged
#' spectral information
#' in list element \code{MIR_mean} (data.frame) and metadata in
#' \code{data_meta} (data.frame).
#' @param remove logical expression (\code{TRUE} or \code{FALSE})
#' that specifies weather spectra shall be removed.
#' If \code{rm = FALSE}, there will be no outlier removal
#' @return Returns list \code{spectra_out} that contains:
#' \itemize{
#'  \item \code{MIR_mean}: Outlier removed MIR spectra as
#'  data.frame object. If \code{remove = FALSE},
#'  the function will
#'  return almost identical list identical to \code{list_spectra},
#'  except that the first \code{indices} column of the spectral
#'  data frame \code{MIR_mean} is removed
#'  (This is done for both options
#'  \code{remove = TRUE} and \code{remove = FALSE}).
#'  \item \code{data_meta}: metadata data.frame, identical
#'  as in the \code{list_spectra} input list.
#'  \item \code{plot_out}: (optional) ggplot2 graph
#'  that shows all spectra (absorbance on x-axis and wavenumber
#'  on y-axis) with outlier marked, if
#'  \code{remove = TRUE}.
#' }
#' @details This is an optional function if one wants to remove
#' outliers.
#' @export
remove_outliers <- function(list_spectra, remove = TRUE) {
  # Outlier detection
  # Use the mvoutlier package and pcout function to identify
  # multivariate outliers
  wfinal01 <- ID <-  NULL
  if (remove == TRUE) {
    # Remove the 'indices' column
    list_spectra$MIR_mean <- list_spectra$MIR_mean[, -1]
    out <- mvoutlier::pcout(list_spectra$MIR_mean, makeplot = T,
      outbound = 0.05) # parameters should be adapted
    # Plot outlying spectra
    plot_out <- plotMIR(
      list_spectra$MIR_mean[
        order(out$wfinal01, decreasing = T), ],
      col = as.factor(out$wfinal01[order(out$wfinal01,
        decreasing = T)])) +
      ggplot2::scale_colour_brewer("outlier", palette = "Set1")
    out_id <- as.character(
      list_spectra$data_meta$ID[!as.logical(out$wfinal01)]
    )
    # Remove  outliers
    MIR_mean <- list_spectra$MIR_mean[
      ! list_spectra$data_meta$ID %in% out_id, ]
    # rep ID and country name
    data_meta <- list_spectra$data_meta[
      ! list_spectra$data_meta$ID %in% out_id, ]
    spectra_out <- list(MIR_mean = MIR_mean,
      data_meta = data_meta,
      plot_out = plot_out)
  } else {
    # Remove the 'indices' column
    list_spectra$MIR_mean <- list_spectra$MIR_mean[, -1]
    spectra_out <- list(MIR_mean = list_spectra$MIR_mean,
      data_meta = list_spectra$data_meta)
  }
  spectra_out
}

## plotMIR function of Antoine Stevens; don't export this
## function to the NAMESPACE
plotMIR <- function(spc, group = NULL, col = NULL,
  linetype = NULL, wr = NULL, brk = NULL,
  ylab = "Absorbance", xlab = "Wavenumber /cm-1",
  by = NULL, by.wrap = T, ...){
  # Function to plot spectra, based on the ggplot2 package
  # spc = spectral matrix, with colnames = wavelengths
  # group = grouping variable, usually the id's of the sample
  # wr = wavelength range to plot
  # brk = breaks of the x-axis
  # by = factor variable for which the mean and sd of
  # each level will be computed and plotted (optional)
  # Requires packages ggplot2; data.table; reshape2
  # Workaround to pass R CMD check:
  # http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
  # Setting the variables to NULL first
  variable <- value <- colour  <- NULL
  spc <- as.data.frame(spc)
  if (!is.null(wr))
    spc <- spc[, as.numeric(colnames(spc)) >= min(wr) &
        as.numeric(colnames(spc)) <= max(wr)]
  if (is.null(brk))
    brk  <- pretty(as.numeric(colnames(spc)), n = 10)
  if (!is.null(by)) {
    spc$by <- by
    spc <- data.table::data.table(spc, check.names = F)
    mean.spc <- reshape2::melt(
      spc[, lapply(data.table::.SD, mean), by = by],
      id.vars = "by"
    )
    sd.spc <- reshape2::melt(
      spc[, lapply(data.table::.SD, sd), by = by],
      id.vars = "by"
    )
    mean.spc$min <- mean.spc$value - sd.spc$value
    mean.spc$max <- mean.spc$value + sd.spc$value
    mean.spc$variable <-  as.numeric(
      as.character(mean.spc$variable)
    )
    if (by.wrap) {
      p <- ggplot2::ggplot(data = mean.spc) +
        ggplot2::geom_ribbon(
          ggplot2::aes(x = variable, ymin = min, ymax = max),
          fill = "grey", col = "black", size = 0.15)  +
        ggplot2::theme_bw()
      p <-  p +  ggplot2::geom_line(
        ggplot2::aes(x = variable, y = value),
        size = 0.25) +
        ggplot2::facet_wrap(~ by) +
        ggplot2::labs(x = xlab, y = ylab) +
        ggplot2::scale_x_reverse(breaks = brk)
    } else {
      p <- ggplot2::ggplot(data = mean.spc,
        ggplot2::aes(x = variable, y = value, group = by, col = by)) +
        ggplot2::geom_line(size = 0.25)  +
        ggplot2::labs(x = xlab, y = ylab) +
        ggplot2::scale_x_reverse(breaks = brk) +
        ggplot2::theme_bw()
    }
    return(p)
  } else {
    if (is.null(group))
      group  <- as.character(1:nrow(spc))
    spc$group <- group
    spc$colour <- col
    spc$linetype <- linetype
    id.var  <- colnames(spc)[
      grep("group|colour|linetype",colnames(spc))]
    tmp <- reshape2::melt(spc, id.var = id.var)
    tmp$variable <- as.numeric(as.character(tmp$variable))
    p <- ggplot2::ggplot(tmp,
      ggplot2::aes(variable, value, group = group)) +
      ggplot2::labs(x = xlab, y = ylab) +
      ggplot2::theme_bw() +
      ggplot2::scale_x_reverse(breaks = brk)
    if (is.null(col) & is.null(linetype))
      p <- p + ggplot2::geom_line(
        ggplot2::aes(colour = group))
    else if (!is.null(col) & is.null(linetype))
      p <- p + ggplot2::geom_line(
        ggplot2::aes(colour = colour))
    else if (is.null(col) & !is.null(linetype))
      p <- p + ggplot2::geom_line(
        ggplot2::aes(colour = group,
        linetype = linetype))
    else  p <- p + ggplot2::geom_line(
      ggplot2::aes(colour = colour,
      linetype = linetype))
    return(p)
  }
}

