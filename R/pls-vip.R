### VIP.R: Implementation of VIP (variable importance in projection)(*) for the
### `pls' package.
### $Id: VIP.R,v 1.2 2007/07/30 09:17:36 bhm Exp $

### Copyright: 2006,2007 Bjoern-Helge Mevik
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License version 2 as
### published by the Free Software Foundation.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.

### A copy of the GPL text is available here:
### http://www.gnu.org/licenses/gpl-2.0.txt

### Contact info:
### Boejrn-Helge Mevik
### bhx6@mevik.net
### Roedtvetvien 20
### N-0955 Oslo
### Norway

### (*) As described in Chong, Il-Gyo & Jun, Chi-Hyuck, 2005, Performance of
### some variable selection methods when multicollinearity is present,
### Chemometrics and Intelligent Laboratory Systems 78, 103--112.

## VIP returns all VIP values for all variables and all number of components,
## as a ncomp x nvars matrix.
VIP <- function(object) {
    # pb: added to avoid `R CMD check` note
    method <- Yloadings <- scores <- loading.weights <- NULL
    if (object$method != "oscorespls")
        stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
    if (nrow(object$Yloadings) > 1)
        stop("Only implemented for single-response models")

    SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
    Wnorm2 <- colSums(object$loading.weights^2)
    SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
    sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}


## VIPjh returns the VIP of variable j with h components
VIPjh <- function(object, j, h) {
    # pb: added to avoid `R CMD check` note
    method <- Yloadings <- scores <- loading.weights <- NULL
    if (object$method != "oscorespls")
        stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
    if (nrow(object$Yloadings) > 1)
        stop("Only implemented for single-response models")

    b <- c(object$Yloadings)[1:h]
    T <- object$scores[,1:h, drop = FALSE]
    SS <- b^2 * colSums(T^2)
    W <- object$loading.weights[,1:h, drop = FALSE]
    Wnorm2 <- colSums(W^2)
    sqrt(nrow(W) * sum(SS * W[j,]^2 / Wnorm2) / sum(SS))
}

#' @title Extract VIPs (variable importance in the projection) for a PLS
#' regression model output returned from model fitting with
#' \code{simplerspec::fit_pls()}
#' @description VIPs are extracted based on the \code{finalModel} sublist
#' in the \code{caret::train} output contained in the \code{model} element
#' of the \code{simplerspec::fit_pls()} model output list. The VIPs for
#' derived number of PLS components in the \code{finalModel} are computed.
#' @param mout Model output list returned from \code{simplerspec::fit_pls()}.
#' @usage extract_pls_vip(mout)
#' @return A tibble data frame with columns \code{wavenumber} and correponding
#' VIP values in the column \code{vip} for the finally chosen PLS regression
#' model at the final number of PLS components.
#' @export
extract_pls_vip <- function(mout) {
  # Compute VIP for all wavenumbers and select only VIPs with ncomp in final
  # model
  final_model <- mout$model$finalModel
  vip <- VIP(object = final_model)[final_model$ncomp, ]
  # Collect wavenumbers from preprocessed spectra
  wn <- as.numeric(colnames(mout$data$calibration$spc_pre[[1]]))
  # Create a data frame with wavenumbers and VIP scores
  tibble::data_frame_(lazyeval::lazy_dots(wavenumber = wn, vip = vip))
}


create_vip_rects <- function(df_vip) {

  # Avoid `R CMD check NOTE`: `no visible binding for global variable`
  VIP <- wavenumber <- tail <- NULL

  # Highlight region data
  # https://stackoverflow.com/questions/32543176/highlight-areas-within-certain-x-range-in-ggplot2
  v <- rep(0, nrow(df_vip))
  v[df_vip$vip > 1] <- 1
  # Get the start and end points for highlighted regions
  inds <- diff(c(0, v))
  start <- df_vip$wavenumber[inds == 1]
  end <- df_vip$wavenumber[inds == -1]
  if (length(start) > length(end)) {
    end <- c(end, tail(df_vip$wavenumber, 1))
  }
  # Create data frame for rectangle layer (geom_rects)
  data.frame(start = start, end = end, group = seq_along(start))
}

#' @title Plot stacked ggplot2 graphs with the Variable Importance for the
#' Projection (VIP) scores, mean replicate spectra (absorbance) per sample_id,
#'and the preprocessed spectra.
#' @description Plot stacked ggplot2 graphs of VIP for the final
#' PLS regression model output of the calibration (training) data set for the
#' final number of components, raw (replicate mean) spectra, and preprocessed
#' spectra. Regions with VIP > 1 are highlighted across the stacked graphs
#' in beige colour rectangles. VIP calculation is implemented as described in
#' Chong, I.-G., and Jun, C.-H. (2005). Performance of some variable selection
#' methods when multicollinearity is present. Chemometrics and Intelligent
#' Laboratory Systems, 78(1--2), 103--112. https://doi.org/10.1016/j.chemolab.2004.12.011
#' @param mout Model output list that is returned from
#' \code{simplerspec::fit_pls()}. This object contains a nested list with
#' the \code{caret::train()} object (class \code{train}), based on which
#' VIPs at finally selected number of PLS components are computed.
#' @param y1 Character vector of list-column name in
#' \code{mout$data$calibration}, where spectra for bottom graph are extracted.
#' Default is \code{"spc_mean"}, which plots the mean calibration spectra after
#' resampling.
#' @param y2 Character string of list-column name in
#' \code{mout$data$calibration}, where spectra for bottom graph are extracted.
#' Default is \code{"spc_pre"}, which plots the preprocessed calibration
#' spectra after resampling.
#' @param by Character string that is used to assign spectra to the same group
#' and therefore ensures that all spectra are plotted with the same colour.
#' Default is \code{"sample_id"}
#' @param xlab Character string of X axis title for shared x axis of stacked
#' graphs. Default is \code{expression(paste("Wavenumber [", cm^-1, "]"))}
#' @param ylab1 Y axis title of bottom spectrum. Default is \code{"Absorbance"}.
#' @param ylab2 Y axis title of bottom spectrum. Default is
#' \code{"Preprocessed Abs."}.
#' @param alpha Double between 0 and 1 that defines transparency of spectra
#' lines in returned graph (ggplot plot object).
#' @usage plot_pls_vip(mout, y1 = "spc_mean", y2 = "spc_pre",
#'   by = "sample_id",
#'   xlab = expression(paste("Wavenumber [", cm^-1, "]")),
#'   ylab1 = "Absorbance", ylab2 = "Preprocessed Abs.",
#'   alpha = 0.2)
#' @export
plot_pls_vip <- function(mout, y1 = "spc_mean", y2 = "spc_pre",
                         by = "sample_id",
                         xlab = expression(paste("Wavenumber [", cm^-1, "]")),
                         ylab1 = "Absorbance", ylab2 = "Preprocessed Abs.",
                         alpha = 0.2) {

  # Avoid `R CMD check` NOTE: no visible binding for global variable `...`
  variable <- wavenumber <- value <- group <- vip <- NULL

  # Extract spectra tibble for calibration
  spc_tbl <- mout$data$calibration
  # Gather spectra in one data.table

  ifelse(y1 == "spc",
    {dt1 <- data.table::data.table(do.call(rbind, spc_tbl[, y1][[y1]]))},
    {dt1 <- data.table::rbindlist(spc_tbl[, y1][[y1]])}
  )

  ifelse(y2 == "spc",
    {dt2 <- data.table::data.table(do.call(rbind, spc_tbl[, y2][[y2]]))},
    {dt2 <- data.table::rbindlist(spc_tbl[, y2][[y2]])}
  )

  # Extract ID variable and append it to the data.table
  id <- spc_tbl[, by][[by]]
  dt1[, id := id]
  dt2[, id := id]

  # Convert data table to long form for ggplot2 plotting
  dt1_long <- data.table::melt(
    dt1, measure = names(dt1)[!names(dt1) %in% c("id")]
  )
  dt1_long[, variable := as.numeric(as.character(variable))]
  data.table::setnames(dt1_long, old = "variable", new = "wavenumber")

  dt2_long <- data.table::melt(
    dt2, measure = names(dt2)[!names(dt2) %in% c("id")]
  )
  dt2_long[, variable := as.numeric(as.character(variable))]
  data.table::setnames(dt2_long, old = "variable", new = "wavenumber")

  # Plot spectra and VIP scores
  brk <- pretty(as.numeric(names(dt1)[!names(dt1) %in% c("id")]),
    n = 10)
  maxmin <- function(x) {c(max(x), min(x))}
  x_lim <- maxmin(as.numeric(names(dt1)[!names(dt1) %in% c("id")]))

  # Extract VIP (Variable Importance in Projection) scores for the given
  # model
  df_vip <- extract_pls_vip(mout)

  # Determine highlighted regions above VIP = 1
  rects <- create_vip_rects(df_vip)

  # Plot for resampled and mean replicate spectra
  p_spc <- ggplot2::ggplot(dt1_long, ggplot2::aes(x = wavenumber, y = value)) +
    ggplot2::geom_rect(data = rects, inherit.aes = FALSE,
      ggplot2::aes(xmin = start, xmax = end, ymin = min(dt1_long$value),
        ymax = max(dt1_long$value), group = group), color = "transparent",
        fill = "orange", alpha = 0.3) +
    ggplot2::geom_line(data = dt1_long, inherit.aes = FALSE,
      ggplot2::aes(x = wavenumber, y = value, group = id),
        alpha = alpha, size = 0.2) +
    ggplot2::labs(x = xlab, y = ylab1) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.margin = ggplot2::unit(c(0, 5, 1, 1),
      units = "points")) +
    ggplot2::scale_x_reverse(limits = x_lim, breaks = brk)

    p_spc_pre <- ggplot2::ggplot(dt2_long, ggplot2::aes(wavenumber, value)) +
    ggplot2::geom_rect(data = rects, inherit.aes = FALSE,
      ggplot2::aes(xmin = start, xmax = end, ymin = min(dt2_long$value),
        ymax = max(dt2_long$value), group = group), color = "transparent",
        fill = "orange", alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(group = id), alpha = alpha, size = 0.2) +
    ggplot2::labs(x = xlab, y = ylab2) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.margin = ggplot2::unit(c(0, 5, 1, 1),
      units = "points")) +
    ggplot2::scale_x_reverse(limits = x_lim, breaks = brk) +
    ggplot2::theme(plot.margin = ggplot2::unit(c(1, 5, -30, 6),
      units = "points"),
      axis.title.y = ggplot2::element_text(vjust = 0.25),
      axis.text.x = ggplot2::element_blank())

  # Plot for VIP
  # Extract PLS model response variable and number of components
  response <- as.character(unique(mout$stats$response))
  ncomp <- mout$model$finalModel$ncomp
  p_vip <- ggplot2::ggplot(data = df_vip,
      ggplot2::aes(x = wavenumber, y = vip)) +
    ggplot2::geom_rect(data = rects, inherit.aes = FALSE,
      ggplot2::aes(xmin = start, xmax = end, ymin = min(df_vip$vip),
      ymax = max(df_vip$vip), group = group), color = "transparent",
      fill = "orange", alpha = 0.3) +
    ggplot2::geom_hline(yintercept = 1, colour = "red") +
    ggplot2::geom_line() +
    ggplot2::ylab(
      bquote(paste(VIP[PLSR], " (", .(response),
        ", ", .(ncomp), " comps)"))) +
    ggplot2::scale_x_reverse(limits = x_lim, breaks = brk) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.margin = ggplot2::unit(c(1, 5, -30, 6),
      units = "points"),
      axis.title.y = ggplot2::element_text(vjust = 0.25),
    axis.text.x = ggplot2::element_blank())

  # Arrange plots in two panels without any margins in between
  # Hints from
  # https://stackoverflow.com/questions/42567045/alignment-of-two-plots-using-grid-arrange
  cowplot::plot_grid(
    p_spc, p_spc_pre, p_vip, rel_heights = c(0.4, 0.3, 0.3),
    ncol = 1, align = "v")
}

