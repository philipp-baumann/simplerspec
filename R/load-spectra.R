## Soil spectroscopy related functions that were compiled by
## Antoine Stevens ==============================================

#' @title Read an OPUS text file
#' @description
#' Read single text file acquired with
#' an Bruker Vertex FTIR Instrument
#' (as exported from OPUS software)
#' @param file.name Character vector with path to files
#' @usage read_opus_text(file.name)
#' @export
read_opus_text <- function(file.name){
  if (file.exists(file.name)) {
    out <- read.csv(file.name, header=F,
      col.names = c("wavenumber", "absorbance")
    )
    return(out)
  } else {
    warning(paste("File", file.name, "does not exist"))
  }
}

#' @title Read an OPUS binary file
#' @description
#' Read single binary file acquired with an
#' Bruker Vertex FTIR Instrument
#' @param file.name Character vector with path to files
#' @usage read_opus_bin(file.name)
#' @export
read_opus_bin <- function(file.name){
  size <- fileRaw <- NULL
  if (file.exists(file.name)) {
    try(
      pa <- hexView::readRaw(file.name, offset = 0,
        nbytes = file.info(file.name)$size, human = "char",
        size = 1, endian = "little"), silent = TRUE)
    if (!class(.Last.value)[1] == "try-error") {

      pr <- pa$fileRaw
      # Get source of instrument
      ins <- grepRaw("INS", pr, all = TRUE)
      ins <- hexView::readRaw(
        file.name, offset = ins[length(ins)] + 7,
        nbytes = 3, human = "char", size = 1, endian = "little"
      )
      ins <- hexView::blockString(ins)
      # Get source of infrared to know if NIR or MIR
      src <- grepRaw("SRC", pr, all = TRUE)
      src <- hexView::readRaw(
        file.name, offset = src[length(src)] + 4,
        nbytes = 3, human = "char", size = 1, endian = "little"
      )
      src <- hexView::blockString(src)
      instr.range <- tolower(paste(ins, src, sep = "-"))
      # Get Beam Splitter
      bms <- grepRaw("BMS", pr, all = TRUE)
      bms <- hexView::readRaw(
        file.name, offset = bms[length(bms)] + 4,
        nbytes = 4, human = "char", size = 1, endian = "little"
        )
      bms <- hexView::blockString(bms)

      z <- grepRaw("ZFF", pr, all = TRUE)[1] + 5
      re <- grepRaw("RES", pr, all = TRUE)[1] + 5
      snm <- grepRaw("SNM", pr, all = TRUE)[1] + 7
      lwn <- grepRaw("LWN", pr, all = TRUE)[1] + 7
      mxy <- grepRaw("MXY", pr, all = TRUE)[1] + 7
      mny <- grepRaw("MNY", pr, all = TRUE)[3] + 7
      end <- grepRaw("END", pr, all = TRUE) + 11
      dat <- grepRaw( "DAT", pr, all = TRUE)[1] + 7
      tim <- grepRaw("TIM", pr, all = TRUE) + 11
      # Select matching positions based on file structure
      # Calculate how many number of points blocks are present in file
      npt_length <- length(grepRaw("NPT", pr, all = TRUE))
      # If spectra are from ETH Alpha
      if (npt_length < 9) {
        fx_sample <- grepRaw("FXV", pr, all = TRUE)[1] + 7
        lx_sample <- grepRaw("LXV", pr, all = TRUE)[1] + 7
        npt_sample <- grepRaw("NPT", pr, all = TRUE)[1] + 7
      # If spectra are from ICRAF Alpha
      } else {
        fx_sample <- grepRaw("FXV", pr, all = TRUE)[3] + 7
        lx_sample <- grepRaw("LXV", pr, all = TRUE)[3] + 7
        npt_sample <- grepRaw("NPT", pr, all = TRUE)[3] + 7
      }

      # byts <- diff(offs)
      ZFF <- hexView::readRaw(file.name, offset = z, nbytes = 4,
        human = "int", size = 2)[[5]][1]
      RES <- hexView::readRaw(file.name, offset = re, nbytes = 4,
        human = "int", size = 2)[[5]][1]
      snm.lab.material <- hexView::blockString(
        hexView::readRaw(file.name, offset = snm, nbytes = 22,
          human = "char", size = 1, endian = "little")
      )
      if (!nzchar(snm.lab.material)) {
        SSN <- ""
        Material <- ""
        warning("Product name not found inside OPUS file...")
      }
      else {
        if (!length(grep(snm.lab.material, pattern = ";")) == 0) {
          snm.lab.material <- as.vector(
            strsplit(snm.lab.material, ";")
            )[[1]]
          SSN <- paste0(snm.lab.material[2], snm.lab.material[1])
          Material <- snm.lab.material[3]
        }   else {
          if (!length(grep(snm.lab.material, pattern = "_")) == 0) {
            # Don't remove "_" from unique id SSN (@baumann)
            # SSN <- sub("_", "", snm.lab.material)
            SSN <- snm.lab.material
            Material <- ""
          } else {
            if (!length(snm.lab.material) == 0) {
              SSN <- snm.lab.material
              Material <- ""
            }
          }
        }
      }
      # Set three SSN first three characters to lower
      # Don't convert to lowercase
      # SSN <- paste0(tolower(substr(SSN, 1, 3)),
      #  substr(SSN, 4, 20))
      Scandate <- hexView::blockString(
        hexView::readRaw(file.name, offset = dat,
        nbytes = 10, human = "char", size = 1,
        endian = "little")
      )
      Scantime <- hexView::blockString(
        hexView::readRaw(file.name,
        offset = tim[2] - 4, nbytes = 8, human = "char",
        size = 1, endian = "little")
      )
      Scandate <- paste(Scandate, Scantime)
      LWN <- hexView::readRaw(
        file.name, offset = lwn, nbytes = 8,
        human = "real", size = 8)[[5]][1]
      # Combine the above parameters
      spectrum.meta <- c(SSN, Material, Scandate, ZFF, RES, LWN)
      # Get number of data points for each spectra data block
      NPT_sample <- hexView::readRaw(
        file.name, offset = npt_sample, nbytes = 12,
        human = "int", size = 4)[[5]][1]
      # Add extra step for spectrum "CI_tb_05_soil_cal.2" measured at ETH
      if(NPT_sample > 10000) {
        npt_sample <- grepRaw("NPT", pr, all = TRUE)[2] + 7
        NPT_sample <- hexView::readRaw(
          file.name, offset = npt_sample, nbytes = 12,
          human = "int", size = 4)[[5]][1]
      }
      # fxv:	Frequency of first point
      fxv_sample <- hexView::readRaw(
        file.name, offset = fx_sample, nbytes = 16,
        human = "real", size = 8)[[5]][1]
      # lxv:	Frequency of last point
      lxv_sample <- hexView::readRaw(
        file.name, offset = lx_sample, nbytes = 16,
        human = "real", size = 8)[[5]][1]
      # Read all through all the data blocks inside the OPUS file
      nbytes_spc_sample <- NPT_sample * 4

      # Offsets (start) of sample spectra
      if (length(end) == 17) {
      offs_spc_sample <- end[16]
      }
      else if (length(end) == 16) {
        offs_spc_sample <- end[15]
      } else {
        offs_spc_sample <- end[7]
      }

      ## Calculate wavenumbers
      wavenumbers <- rev(seq(lxv_sample, fxv_sample,
        (fxv_sample - lxv_sample)/(NPT_sample - 1)))

      # Read spectra
      spectra <- hexView::readRaw(file.name, width = NULL,
        offset = offs_spc_sample, nbytes = nbytes_spc_sample,
        human = "real",
        size = 4, endian = "little")[[5]]

      # File name
      file_name <- sub(".+/(.+)", "\\1", file.name)

      # Create date_time object
      date_time <- as.POSIXct(spectrum.meta[3],
        format = "%d/%m/%Y %H:%M:%S ")

      # Create unique_id using file_name and time
      ymd_id <- format(date_time, "%Y%m%d")
      unique_id <- paste0(file_name, "_", ymd_id)

      # Add sample_id: remove extension .0, .1 etc. from OPUS files
      sample_id <- sub("(.+)\\.[[:digit:]]+$", "\\1", file_name)

      # Extract repetition number (rep_no) from file name
      rep_no <- sub(".+\\.([[:digit:]])+$", "\\1", file.name)

      # Convert spectra to matrix and add dimnames (wavenumbers for columns
      # and unique_id for rows)
      spc_m <- matrix(spectra, ncol = length(spectra), byrow = FALSE)
      rownames(spc_m) <- unique_id
      colnames(spc_m) <- round(wavenumbers, 1)

      out <- list(
        metadata = data.frame(
          unique_id = unique_id,
          scan_id = file_name, # changed file_name to scan_id in output list
          sample_id = sample_id,
          rep_no = rep_no,
          date_time = date_time,
          sample_info = spectrum.meta[1],
          instrument_name = instr.range,
          resolution = spectrum.meta[5],
          bms = bms,
          lwn = spectrum.meta[6]
          ),
        spc = spc_m,
        wavenumbers =  wavenumbers
      )

      # names(out)[-c(1:9)] <- as.character(round(wavenumbers, 1))
      return(out)
    }
  } else {
    warning(paste("File", file.name, "does not exist"))
  }
}

#' @title Read OPUS binary and ASCII files
#' @description
#' Read single or multiple binary and ASCII files acquired with
#' an Bruker Vertex FTIR Instrument
#' @usage
#' read_opus(fnames, in_format, out_format)
#' @param fnames character \code{vector} of the name(s)
#' (with absolute path) of the file(s) to read
#' @param in_format format of the input file: \code{'binary'} or
#' \code{'txt'}
#' @param out_format format of the output:
#' \code{'matrix'} (default) or \code{'list'} (see below)
#' @return
#' if \code{out_format} = \code{'matrix'}, absorbance values
#' of the input file(s) in a single \code{matrix}.
#'
#' if \code{out_format} = \code{'list'}, a \code{list} of the
#' input file(s) data consisting of a \code{list} with components:
#' \itemize{
#'  \item{\code{Name}}{ name of the file imported}
#'  \item{\code{datetime}}{ date and time of acquisition in
#'  \code{POSIXct} format (available only when
#'  \code{in_format} = 'binary')}
#'  \item{\code{metadata}}{ \code{list} with information
#'  on instrument configuration (available only when
#'  \code{in_format} = 'binary')}
#'  \item{\code{absorbance}}{  a numeric \code{vector}
#'  of absorbance values}
#'  \item{\code{wavenumbers}}{ numeric \code{vector}
#' of the band positions}
#' }
#' @author Antoine Stevens and Andrew Sila (soil.spec package)
#' @note
#' This is essentially a re-factored and simplified version of
#' the \code{read.opus} function from the
#' \sQuote{soil.spec} package for reading OPUS VERTEX files
#' The function should also work for other OPUS files (eg alpha),
#' see \code{read.opus}.
#' @export
read_opus <- function(fnames, in_format = c("binary", "txt"),
  out_format = c("matrix", "list")) {
  # hexView and plyr are required

  wavenumbers <- NULL
  absorbance <- NULL

  in_format <- match.arg(in_format)
  out_format <- match.arg(out_format)

  spc <- vector("list", length(fnames))
  i <- 1
  for (file.name in fnames) {
    if (in_format == "binary") {
      spc[[i]] <- read_opus_bin(file.name)
    } else {
      spc[[i]] <- read_opus_text(file.name)
    }
    i <- i + 1
  }
  names(spc) <- sub(".+/(.+)(\\.txt)?$", "\\1", fnames)
  if (out_format == "matrix") {
    test <- sapply(spc, function(x) class(x) != "character")
    # warning(
    # paste0(paste(names(spc)[!test], collapse = ","),
    # " do not exist")
    # )
    spc <- spc[test]
    if(in_format == "binary"){
      spc <- do.call(plyr::rbind.fill, lapply(spc, function(x){
        x <- t(data.frame(
          wav = x$wavenumbers, absorbance = x$absorbance))
        colnames(x) <- x[1,]
        data.frame(x[2, , drop = F], check.names = F)}))

    } else {
      spc <- do.call(plyr::rbind.fill, lapply(spc, function(x) {
        x <- t(x)
        colnames(x) <- x[1,]
        data.frame(x[2, , drop = F], check.names = F)}))
    }
    rownames(spc) <- sub(".+/(.+)(\\.txt)?$", "\\1", fnames)

  }
  return(spc)

}
