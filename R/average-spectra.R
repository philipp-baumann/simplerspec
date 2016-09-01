# Helper function written by Antoine Stevens that
# is used for averaging replication scans of one sample
#' @import data.table
#' @import Rcpp
by_spc <- function(spc, indices, fun = mean){
  # Fast summary of spectral data
  # spc = spectral matrix
  # indices  = factor variable used to summarize data
  # fun = summary function
  # Avoid NOTE "no visible binding for global variable '.SD'"
  # when checking package by devtools::check()
  # .SD <- NULL
  spc <- data.table::data.table(indices, spc, check.names = F)
  if(is.null(ncol(indices))){
    x <- 1
  } else {
    x <- ncol(indices)
  }
  as.data.frame(spc[, lapply(.SD, fun),
    by = eval(names(spc)[1:x])])
}
#' @title Calculate mean of spectra
#' @description Calculate the mean of each spectral repetitions
#' (absorbance average per wavenumber)
#' @param in_spectra List that contains spectral data in the
#' element \code{MIR} (data.frame) and sample metadata in the
#' list element \code{data_rep} (data.frame).
#' The data.frame \code{data_meta}
#' contains the sample ID stored in the \code{ID}
#' vector (originally from spectral file names),
#' country abbreviation stored in \code{contry} (2 letters),
#' and the vector \code{site} (2 letters) that is the country
#' abbreviation.
#' @return \code{out_spectra}: List that contains:
#' \itemize{
#'  \item \code{data_meta}: metadata of sample (data.frame)
#'  that is
#'  taken from the element \code{rep} of the input list argument
#'  \code{in_spectra}
#'  \item \code{MIR_mean}: average spectra from replicates of
#'   sample ID
#'  (data.frame)
#'  \item \code{MIR_sd}: standard deviation of spectra calculated
#'  from replicates of sample ID (data.frame)
#'  \item \code{cvar} coefficient of variance over all
#'  wavenumbers of spectra
#'  calculated from replicates of sample ID (vector)
#' }
#' @export
average_spectra <- function(in_spectra) {
  # Compute mean per sample with by_spc,
  # provided by Antoine
  # spc = spectral data
  # indices = character vector(s) or factor(s) to group the rows
  # by fun = summary function
  # Also compute the standard deviation (SD)
  # of the three measurements
  # Identify samples in which the spectrum has SD higher than 1.5
  # and need to be re-scanned - ?
  MIR <- NULL
  data_rep <- NULL
  ID <- NULL
  MIR_mean <- by_spc(spc = in_spectra$MIR[, ],
    indices = in_spectra$data_rep$ID[], fun = mean)
  MIR_sd <- by_spc(spc = in_spectra$MIR[, ],
    indices = in_spectra$data_rep$ID[], fun = sd)
  # MIR_sd[order(rowMeans(MIR_sd[, -1])), 1]
  # rowMeans(MIR_sd[, -1])[order(rowMeans(MIR_sd[, -1]))] /
  # rowMeans(MIR_mean[, -1])[order(rowMeans(MIR_sd[, -1]))]
  # compute the coefficient of variation
  cvar <- rowMeans(MIR_sd[, -1])/rowMeans(MIR_mean[, -1])
  # add metadata for each sample; take strings of rownames
  # from spectra; first two characters of the sample_ID code
  # is country, character pos. 4 and 5 is site abbreviation
  data_meta <- data.frame(ID = MIR_mean[, 1])
  data_meta <- cbind(data_meta,
    country = substring(data_meta$ID, first = 1, last = 2),
    site = substring(data_meta$ID, first = 4, last = 5)
  )
  out_spectra <- list(data_meta = data_meta,
    MIR_mean = MIR_mean,
    MIR_sd = MIR_sd,
    cvar = cvar)
  return(out_spectra)
}
