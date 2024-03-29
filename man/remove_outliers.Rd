% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove-outl-spectra.R
\name{remove_outliers}
\alias{remove_outliers}
\title{Remove outlier spectra}
\usage{
remove_outliers(list_spectra, remove = TRUE)
}
\arguments{
\item{list_spectra}{List that contains averaged
spectral information
in list element \code{MIR_mean} (data.frame) and metadata in
\code{data_meta} (data.frame).}

\item{remove}{logical expression (\code{TRUE} or \code{FALSE})
that specifies weather spectra shall be removed.
If \code{rm = FALSE}, there will be no outlier removal}
}
\value{
Returns list \code{spectra_out} that contains:
\itemize{
\item \code{MIR_mean}: Outlier removed MIR spectra as
data.frame object. If \code{remove = FALSE},
the function will
return almost identical list identical to \code{list_spectra},
except that the first \code{indices} column of the spectral
data frame \code{MIR_mean} is removed
(This is done for both options
\code{remove = TRUE} and \code{remove = FALSE}).
\item \code{data_meta}: metadata data.frame, identical
as in the \code{list_spectra} input list.
\item \code{plot_out}: (optional) ggplot2 graph
that shows all spectra (absorbance on x-axis and wavenumber
on y-axis) with outlier marked, if
\code{remove = TRUE}.
}
}
\description{
Remove outlier spectra based on the
\code{pcout()} function of the \code{mvoutlier} package.
}
\details{
This is an optional function if one wants to remove
outliers.
}
