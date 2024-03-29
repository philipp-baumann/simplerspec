% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read-asd.R
\name{read_asd}
\alias{read_asd}
\title{Read ASD fieldspec spectrometer data export into into simplerspec
spectra tibble.}
\usage{
read_asd(file)
}
\arguments{
\item{file}{Tab delmited file from ASD software export where the first
column called \code{Wavelength} contais wavelengths in nanometer and the
remaining columns are sample spectra referred by an ID name provided in the
first row of these columns.}
}
\value{
Spectra data in tibble data frame (class \code{tbl_df}) that contains
columns \code{sample_id} (derived from 2nd and following column names of
tab delimited ASD exported text file),
\code{spc} (list-column of spectral matrices)
and \code{wavelengths} (list-column containing wavelength vectors).
}
\description{
Read tab delimited text (.txt) files exported from ASD field
spectrometer into simplerspec spectra tibble.
ASD Fieldspec data files are expected in .txt tab delimited file format.
The first row should contain
the name 'Wavelength' for the first column and the file names for the
remaining columns.
}
