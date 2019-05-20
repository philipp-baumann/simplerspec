#' @title Read a Bruker OPUS spectrum binary file
#' @description
#' Read single binary file acquired with an
#' Bruker Vertex FTIR Instrument
#' @param file_path Character vector with path to file
#' @param extract Character vector of spectra types to extract from OPUS binary
#' file. Default is \code{c("spc")}, which will extract the final spectra, e.g.
#' expressed in absorbance (named \code{AB} in Bruker OPUS programs). Possible
#' additional values for the character vector supplied to extract are
#' \code{"ScSm"} (single channel spectrum of the sample measurement), \
#' code{"ScRf"} (single channel spectrum of the reference measurment),
#' \code{"IgSm"} (interferogram of the sample measurment) and \code{"IgRf"}
#' (interferogram of the reference measurement).
#' @param print_progress Logical (default \code{TRUE}) whether a message is
#' printed when an OPUS binary file is parsed into an R list entry.
#' @param atm_comp_minus4offset Logical whether spectra after atmospheric
#' compensation are read with an offset of \code{-4} bites from Bruker OPUS
#' files. Default is \code{FALSE}.
#' @usage read_opus_bin_univ(file_path, extract = c("spc"),
#' print_progress = TRUE, atm_comp_minus4offset = FALSE)
# Importing functions `%do%` and foreach::`%dopar%` does not work, see
# http://stackoverflow.com/questions/30216613/how-to-use-dopar-when-only-import-foreach-in-description-of-a-package
# Got the following error:
# "Error : object '`%do%`' is not exported by 'namespace:foreach'"
#' @importFrom foreach %dopar% %do%
#' @export
read_opus_bin_univ <- function(file_path, extract = c("spc"),
                               print_progress = TRUE,
                               atm_comp_minus4offset = FALSE) {

  # Avoid `R CMD check` NOTE: no visible binding for global variable ...
  x <- y <- i <- npt <- NULL
  if (!file.exists(file_path)) {
    stop(paste0("File does not exist"))
  }

  try({

    # file_path <- "data/soilspec_background/yamsys_bg_gold/BF_lo_15_soil_cal.0"
    # Read entire content of file as bytes
    pa <- hexView::readRaw(file_path, offset = 0,
      nbytes = file.info(file_path)$size, human = "char",
      size = 1, endian = "little")
    # Get raw vector
    pr <- pa$fileRaw

    # Read byte positions for selected 3 letter strings that flag important
    # spectral information -----------------------------------------------------

    # Get positions of "END" strings
    end <- grepRaw("END", pr, all = TRUE) + 11
    # Get all positions of "NPT" (number of points) string
    npt_all <- grepRaw("NPT", pr, all = TRUE) + 3
    # Get frequency of first (FXV) and last point (LXV) positions
    fxv_all <- grepRaw("FXV", pr, all = TRUE) + 7
    lxv_all <- grepRaw("LXV", pr, all = TRUE) + 7

    # For some files, the number of positions where "FXV" and "LXV" occur
    # are not equal, e.g. for the file in
    # data/soilspec_esal_bin/BF_mo_01_soil_cal.0 ; As a consequence, the
    # fist and last point numbers (e.g. wavenumber or points for interferograms)
    # are not correctly read. This results in an error when trying to calculate
    # the wavenumbers; The below code is a quick and dirty fix to remove
    # FXV values that don't have LXV values and vice versa
    # (difference between "LXV" and "FXV" for a spectral data block
    # should be 16) ------------------------------------------------------------
    if (length(fxv_all) > length(lxv_all)) {
      diff_lxv_fxv <- lapply(lxv_all, function(x) x - fxv_all)
      # Return list of logical vectors indicating whether difference of fxv
      # and lxv is 16 (distance of 16 bytes)
      lxv_fxv_min <- lapply(diff_lxv_fxv, function(x) x == 16)
      fxv_list <- rep(list(fxv_all), length(fxv_all))
      fxv_all <- foreach::foreach(
        x = 1:length(fxv_list), y = 1:length(lxv_fxv_min),
        .combine = 'c') %do% {
          fxv_list[[x]][lxv_fxv_min[[y]]]
        }
    }

    if (length(lxv_all) > length(fxv_all)) {
      diff_fxv_lxv <- lapply(fxv_all, function(x) x - lxv_all)
      # Return list of logical vectors indicating whether difference of fxv
      # and lxv is 16 (distance of 16 bytes)
      fxv_lxv_min <- lapply(diff_fxv_lxv, function(x) x == -16)
      lxv_list <- rep(list(lxv_all), length(lxv_all))
      lxv_all <- foreach::foreach(
        x = 1:length(lxv_list), y = 1:length(fxv_lxv_min),
        .combine = 'c') %do% {
          lxv_list[[x]][fxv_lxv_min[[y]]]
        }
    }

    # Reduce size of npt_all -----------------------------------------------------
    # Error in reading file: "data/spc_NABO_error/30014 KB014 1-14A O-1_CN0104.4"
    if (length(npt_all) != length(fxv_all)) {
      diff_npt_fxv <- lapply(npt_all, function(x) fxv_all - x)
      is16 <- lapply(diff_npt_fxv, function(x) x == 16)
      which16 <- sapply(is16, function(x) any(x == TRUE))
      npt_all <- npt_all[which16]
    }

    # --------------------------------------------------------------------------

    ## Read basic spectral information =========================================

    # Read all number of points (NPT) at once
    NPT <- foreach::foreach(npt = npt_all, .combine = 'c') %do% {
      hexView::readRaw(
        file_path, offset = npt, nbytes = 12, human = "int", size = 4)[[5]][2]
    }

    # Specific error for file: <"data/soilspec_eth_bin/CI_tb_05_soil_cal.2">
    # "Invalid number of bytes" when trying to read spectra
    # -> Reason: NPT at position 1 is 995236000 !!!
    # Omit this entry in NPT and corresponding byte position in npt_all
    # Quick fix ----------------------------------------------------------------
    npt_all <- npt_all[NPT < 40000]
    NPT <- NPT[NPT < 40000]
    # --------------------------------------------------------------------------

    # Figure out how many spectral blocks exist and select final spectra
    # positions; end_spc is vector of offsets where spectra start
    end_spc <- end[diff(end) > 4 * min(NPT)]

    ## Find final spectra information block positions
    ## that belong to spectra data =============================================

    # Save positions that contain possible spectra data block
    # standard parameters
    spc_param_list <- list(
      'npt' = npt_all,
      'fxv' = fxv_all,
      'lxv' = lxv_all
    )

    ## Return list of final parameters corresponding to data blocks that contain
    ## spectra, elements are npt (number of points),
    ## fxv (frequency of first point) and lxv (frequency of last point);
    ## returned values represent byte positions in the file where spectra
    ## parameters are stored. --------------------------------------------------
    return_spc_param <- function(end_spc, spc_param_list) {

      # Difference between any NPT position vector elements end_spc element
      # (end_spc[i] is a scalar, constant value at iteration i)
      diff_l <- lapply(end_spc, function(x) npt_all - x)
      # Test of any vector in list contains -164 (returns list of vectors
      # TRUE or FALSE)
      isminus164 <- lapply(diff_l, function(x) x == -164)

      # Find minimum positive difference within each list
      sel_min <- lapply(diff_l,
        function(x) {if (any(x > 0)) {x == min(x[x > 0])} else {x == -164}}
      )
      # Set FALSE repeated vector in sel_min element where TRUE positions are
      # duplicated
      which_elem_dupl <- which(duplicated(sapply(sel_min, which)))
      if (length(which_elem_dupl) > 1) {
        sel_min[which_elem_dupl] <- NULL
        # Reduce end_spc with duplicated elements
        end_spc <- end_spc[- which_elem_dupl]
      }

      # Select minimum difference NPT position for each END position
      npt_min <- Map(function(x, y) x[y],
        rep(list(npt_all), length(end_spc)), sel_min)
      npt_min <- Filter(length, npt_min)

      # Select spectra parameters that immediately follow END positions before
      # corresponding spectra
      param_min <- foreach::foreach(i = 1:length(spc_param_list),
        .final = function(i) setNames(i, names(spc_param_list))) %do% {
          Map(function(x, y) x[y],
            rep(list(spc_param_list[[i]]), length(end_spc)), sel_min)
      }

      # Test if any difference in list is -164
      if (any(unlist(isminus164) == TRUE)) {
        # Find all list element that contain TRUE in logical vector
        minus164 <- lapply(isminus164, function(x) Find(isTRUE, x))
        # Return element position of last TRUE in list
        where <- function(f, x) {
          vapply(x, f, logical(1))
        }
        last_minus164 <- Position(isTRUE, where(isTRUE, minus164),
          right = TRUE)
        # Replace positions in parameter list are at positions of last
        # -164 difference between end_spc element and NPT position
        param_min <- foreach::foreach(i = 1:length(spc_param_list),
          .final = function(i) setNames(i, names(spc_param_list))) %do% {
          param_min[[i]][[last_minus164]] <-
            spc_param_list[[i]][isminus164[[last_minus164]]]
          param_min[[i]]
        }
      }
      # Return list of final parameters corresponding to data blocks that
      # contain spectra
      param_spc <- lapply(param_min, unlist)
      param_spc$end_spc <- end_spc
      param_spc
    }
    # Save spectra parameter list
    param_spc <- return_spc_param(end_spc, spc_param_list)

    # Create individual vectors containing spectra parameters
    npt_spc <- param_spc[["npt"]]
    fxv_spc <- param_spc[["fxv"]]
    lxv_spc <- param_spc[["lxv"]]
    end_spc <- param_spc[["end_spc"]]

    # Read number of points corresponding to spectra in file -------------------

    NPT_spc <- foreach::foreach(i = 1:length(npt_spc), .combine = 'c') %do% {
      hexView::readRaw(
        file_path, offset = npt_spc[i],
        nbytes = 12, human = "int", size = 4)[[5]][2]
    }

    # Delete NPT with negative signs
    NPT_spc <- NPT_spc[NPT_spc > 0]

    ## Read all spectra ========================================================

    spc <- Map(function(end, NPT) hexView::readRaw(file_path, width = NULL,
      offset = end - 4, nbytes = NPT * 4,
      human = "real", size = 4, endian = "little")[[5]], end_spc, NPT_spc)

    # Read FXV and LXV and calculate wavenumbers  ------------------------------

    FXV_spc <- foreach::foreach(i = 1:length(fxv_spc), .combine = 'c') %do% {
      hexView::readRaw(file_path,
        offset = fxv_spc[i], nbytes = 16, human = "real", size = 8)[[5]][1]
    }
    LXV_spc <- foreach::foreach(i = 1:length(lxv_spc), .combine = 'c') %do% {
      hexView::readRaw(file_path,
        offset = lxv_spc[i], nbytes = 16, human = "real", size = 8)[[5]][1]
    }
    # Calculate wavenumbers
    wavenumbers <- foreach::foreach(i = 1:length(FXV_spc)) %do% {
      rev(seq(LXV_spc[i], FXV_spc[i],
        (FXV_spc[i] - LXV_spc[i]) / (NPT_spc[i] - 1)))
    }

    ## Assigning list of intially read spectra depending on block type =========

    # Assign an index name to the spectra and parameters for reading
    names(end_spc) <- paste0("idx", 1:length(end_spc))
    names(spc) <- paste0("idx", 1:length(spc))
    names(NPT_spc) <- paste0("idx", 1:length(NPT_spc))
    names(FXV_spc) <- paste0("idx", 1:length(FXV_spc))
    names(wavenumbers) <- paste0("idx", 1:length(wavenumbers))

    # Check if elements in FXV_spc (frequency of first point) are equal to 0;
    # these are interferogram spectra ------------------------------------------
    which_Ig <- FXV_spc[which(FXV_spc == 0)]
    Ig_assigned <- if (length(which_Ig) == 0) {
        NULL
      } else if (length(which_Ig) == 1) {
        list(
          spc_idx = names(which_Ig),
          spc_code = "IgSm"
        )
      } else if (length(which_Ig) == 3) {
        list(
          spc_idx = names(which_Ig)[c(1, 3)],
          spc_code = c("IgSm", "IgRf")
        )
      } else {
      list(
        spc_idx = names(which_Ig),
        spc_code = c("IgSm", "IgRf")
      )
    }

    na_assigned <- list(
      spc_idx = NULL,
      spc_code = NULL
    )
    if (length(which_Ig) == 3) {
      # Assign duplicated interferogram spectrum to 'not available' assigned
      na_assigned <- list(
        spc_idx = names(which_Ig)[2],
        spc_code = NA
      )
    }

    # Remove NA assigned spectra in spc list -------------------------------------
    if (!is.null(na_assigned$spc_idx)) {
      spc[na_assigned$spc_idx] <- NULL
    # Remove wavenumbers with NA assigned spectra in spc list
      wavenumbers[na_assigned$spc_idx] <- NULL
    }

    # Assign single channel spectra if present in file -------------------------
    # Return idx (index names) of all remaining spectra that are not
    # interferograms
    notIg <- names(spc)[!names(spc) %in%
      c(Ig_assigned$spc_idx, na_assigned$spc_idx)]
    # Calculate peak ratio for absorbance at around 2392 cm^(-1)
    # and 2358 cm^(-1)
    peak_ratio <- lapply(
      lapply(names(wavenumbers[notIg]),
        function(i) spc[[i]][wavenumbers[notIg][[i]] < 2392 &
          wavenumbers[notIg][[i]] > 2358]),
      function(j) j[[1]] / j[[length(j)]]
    )
    names(peak_ratio) <- names(spc[notIg])
    # Single channel (Sc) assignment list
    which_Sc <- names(which(peak_ratio > 2))
    # Check for single channel, exclude spectral blocks already assigned to
    # interferograms
    Sc_assigned <- if (length(which_Sc) == 0) {
        NULL
      } else if (length(which_Sc) == 1) {
        list(
          spc_idx = which_Sc,
          spc_code = "ScSm"
      )
      } else {
        list(
          spc_idx = which_Sc,
          spc_code = c("ScSm", "ScRf")
      )
      }
    # Assign corrected and uncorrected (if present) ----------------------------
    # AB spectra list
    which_AB <- names(spc)[!names(spc) %in%
      c(Ig_assigned[["spc_idx"]], na_assigned[["spc_idx"]],
        Sc_assigned[["spc_idx"]])]
    AB_assigned <- if (length(which_AB) == 1) {
      list(
        spc_idx = which_AB,
        spc_code = "spc"
      )
    } else {
      list(
        spc_idx = which_AB,
        spc_code = c("spc_nocomp", "spc")
      )
    }

    # Read result spectrum with new offset (no `-4`) when atmospheric
    # compensation was done by the OPUS software; replace the spectrum position
    # with index name idx that corresponds to final spectrum after atmospheric
    # compensation; OPUS files from particular spectrometers/OPUS software
    # versions do still need the same offset end_spc[[spc_idx]] - 4 as the other
    # spectra types; new argument atm_comp_minus4offset (default FALSE) is a
    # quick fix to read files with different offsets after atmospheric
    # compensation -------------------------------------------------------------
    if (length(which_AB) == 2 && !atm_comp_minus4offset) {
      spc[[which_AB[length(which_AB)]]] <-
        hexView::readRaw(file_path, width = NULL,
        offset = end_spc[which_AB[length(which_AB)]],
        nbytes = NPT_spc[which_AB[length(which_AB)]] * 4,
        human = "real", size = 4, endian = "little")[[5]]
    }

    # Assign spectra type for final spectra in element names of spc list -------
    # Combine spectral assignments lists
    list_assigned <- list(
      'Ig' = Ig_assigned,
      'Sc' = Sc_assigned,
      'AB' = AB_assigned
    )
    # Transpose spectra assignment list, first remove NULL elements in list
    list_assigned_t <- purrr::transpose(
      Filter(Negate(function(x) is.null(unlist(x))), list_assigned)
    )
    # Save spectra index (spc_idx) and spectra code (spc_code)
    # in character vector
    spc_idx <- unlist(list_assigned_t[["spc_idx"]])
    spc_code <- unlist(list_assigned_t[["spc_code"]])
    # Order spc_idx from 1 to n spectra (n = length of end_spc)
    order_spc <- as.numeric(
      sub(".*idx", "", unlist(list_assigned_t[["spc_idx"]])))
    spc_type <- spc_code[order(order_spc)]
    # Set spectrum type as element names of spectra list (spc)
    names(spc) <- spc_type
    # Set spectrum type in wavenumbers list
    names(wavenumbers) <- spc_type

    # Read with new offset when first value of
    # ScSm  single channel sample spectrumspectrum is 0 and replace previous ---
    if (any(names(spc) %in% "ScSm" & spc[["ScSm"]][1] == 0)) {
      spc[["ScSm"]] <-
        hexView::readRaw(file_path, width = NULL,
          offset = end_spc[Sc_assigned$spc_idx[Sc_assigned$spc_code == "ScSm"]],
          nbytes = NPT_spc[Sc_assigned$spc_idx[Sc_assigned$spc_code == "ScSm"]]
            * 4,
          human = "real", size = 4, endian = "little")[[5]]
    }

    ## Get additional parameters from OPUS binary file =========================

    # Instrument parameters ----------------------------------------------------
    ins <- grepRaw("INS", pr, all = TRUE) # Instrument type
    INS <- hexView::blockString(
      hexView::readRaw(
        file_path, offset = ins[length(ins)] + 7,
        nbytes = 10, human = "char", size = 1, endian = "little"))
    lwn <- grepRaw("LWN", pr, all = TRUE)[1] + 7 # Laser wavenumber
    LWN <- hexView::readRaw(file_path, offset = lwn,
      nbytes = 8, human = "real", size=8)[[5]][1]
    tsc <- grepRaw("TSC", pr, all = TRUE) + 7 # Scanner temperature
    TSC_all <- lapply(tsc, function(tsc)
      hexView::readRaw(file_path, offset = tsc,
        nbytes = 16, human = "real", size = 8)[[5]][[1]] # can include sample
      # and background temperature
    )
    # Read relative humidity of the interferometer during measurement
    hum_rel <- grepRaw("HUM", pr, all = TRUE) + 7
    HUM_rel <- lapply(hum_rel, function(hum_rel)
      hexView::readRaw(
        file_path, offset = hum_rel, nbytes = 16,
        human = "int", size = 8)[[5]][[1]] # can include sample and background
        # humidity
    )
    # Read absolute humidity of the interferometer during measurement
    hum_abs <- grepRaw("HUA", pr, all = TRUE) + 7
    HUM_abs <- lapply(hum_abs, function(hum_abs)
      hexView::readRaw(
        file_path, offset = hum_abs, nbytes = 16,
        human = "real", size = 8)[[5]][[1]] # can include sample and background
        # humidity
    )

    # Optics parameters --------------------------------------------------------
    src <- grepRaw("SRC", pr, all = TRUE) # Source: MIR or NIR
    SRC <- hexView::blockString(
      hexView::readRaw(
        file_path, offset = src[length(src)] + 4,
        nbytes = 3, human = "char", size = 1, endian = "little"))
    instr_range <- tolower(paste(INS, SRC, sep = "-")) # instrument range
    bms <- grepRaw("BMS", pr, all = TRUE) # Beamsplitter
    BMS <- hexView::blockString(
      hexView::readRaw(file_path, offset = bms[length(bms)] + 4,
        nbytes = 3, human = "char", size = 1, endian = "little"))

    # Fourier transform parameters ---------------------------------------------
    zff <- grepRaw("ZFF", pr, all = TRUE)[1] + 5 # Zero filling factor (numeric)
    ZFF <- hexView::readRaw(file_path, offset = zff,
      nbytes = 4, human = "int", size=2)[[5]][1]

    # (Additional) Standard parameters -----------------------------------------
    csf_all <- grepRaw("CSF", pr, all = TRUE) + 7 # y-scaling factor
    # Read only CSF byte positions that correspond to final spectra
    CSF <- lapply(csf_all[npt_all %in% npt_spc],
      function(csf) hexView::readRaw(
        file_path, offset = csf, nbytes = 8, human = "real", size = 8)[[5]][1])
    mxy_all <- grepRaw("MXY", pr, all = TRUE) + 7 # Y-maximum
    MXY <- unlist(lapply(mxy_all[npt_all %in% npt_spc],
      function(mxy) hexView::readRaw(
        file_path, offset = mxy, nbytes = 8, human = "real", size = 8)[[5]][1]))
    mny <- grepRaw("MNY", pr, all = TRUE) + 7 # Y-minimum
    dxu_all <- grepRaw("DXU", pr, all = TRUE) + 7 # X units
    DXU <- lapply(dxu_all, function(dxu)
      hexView::blockString(
        hexView::readRaw(file_path, offset = dxu,
        nbytes = 3, human = "char", size = 1, endian = "little")
      )
    )
    # Y units -> there is no DYU present in file
    dyu_all <- grepRaw("DYU", pr, all = TRUE) + 7
    dat <- grepRaw("DAT", pr, all = TRUE) + 7 # Date
    tim <- grepRaw("TIM", pr, all = TRUE) + 7 # Time
    time <- unlist(lapply(tim, function(tim)
      hexView::blockString(
        hexView::readRaw(file_path, offset = tim,
          nbytes = 22, human = "char",
        size = 1, endian = "little")))
    )
    # Only select "DAT" string positions that are immediately before time
    dat_sel <- foreach::foreach(i = 1:length(tim), .combine = 'c') %do% {
      diff_sel <- dat - tim[i]
      dat[which(diff_sel <= 32 & diff_sel >= -20)]
    }
    date <- lapply(dat_sel, function(dat) hexView::blockString(
      hexView::readRaw(file_path, offset = dat,
      nbytes = 10, human = "char", size = 1,
      endian = "little"))
    )

    date_time <- unique(paste(date, time))
    # Convert date_time from character to class POSIXct (calendar date and time)
    date_time <- as.POSIXct(date_time, format = "%d/%m/%Y %H:%M:%S")
      # , tz = "GMT+1") # tz is argument for time zone

    # Scale all spectra with y-scaling factor if any of spectra types present
    # in file are not 1 --------------------------------------------------------
    # Set names of CSF elements equal to spectra list element names
    names(CSF) <- names(spc)
    if (any(unlist(CSF) != 1)) {
      # Return all elements in CSF that have scaling value not equal to 1
      CSF_toscale <- Filter(function(x) x != 1, CSF)
      # Apply scaling for spectra with CSF value not equal to 1;
      # Map() returns list
      spc_scaled <- Map(function(CSF, spc) CSF * spc,
        unlist(CSF_toscale), spc[names(CSF_toscale)])
      # Replace all spc list elements that have CSF not equal 1 with
      # scaled values
      spc <- replace(x = spc, list = names(CSF_toscale), values = spc_scaled)
    }

    # Data aquisition parameters -----------------------------------------------

    plf <- grepRaw("PLF", pr, all = TRUE) + 4 # Result spectrum
    PLF_all <- lapply(plf, function(plf) hexView::blockString(
      hexView::readRaw(file_path, offset = plf,
        nbytes = 2, human = "char", size = 1,
        endian = "little"))
    )
    # Select only result spectra abbreviations that are more than 0 characters
    # long
    PLF <- unlist(PLF_all[lapply(PLF_all, nchar) > 0])
    res <- grepRaw("RES", pr, all = TRUE)[1] + 5 # Resolution (wavenumber)
    RES <- hexView::readRaw(
      file_path, offset = res, nbytes = 4, human = "int", size = 2)[[5]][1]

    ## Create sample metadata objects ==========================================
    # File name
    file_name_nopath <- sub(".+/(.+)", "\\1", file_path)
    # Create sample id from file name;
    # remove extension .0, .1 etc. from OPUS files
    sample_id <- sub("(.+)\\.[[:digit:]]+$", "\\1", file_name_nopath)
    # Extract sample repetition number (rep_no) from file name
    rep_no <- sub(".+\\.([[:digit:]])+$", "\\1", file_path)
    snm <- grepRaw("SNM", pr, all = TRUE)[1] + 7
    SNM <- hexView::blockString(
      hexView::readRaw(file_path, offset = snm,
        nbytes = 30, human = "char", size = 1, endian = "little")
    )
    # Create unique_id using file_name and time
    # ymd_id <- format(max(date_time), "%Y%m%d")
    ymdhms_id <- max(date_time)
    unique_id <- paste0(file_name_nopath, "_", ymdhms_id)

    ## Convert all spectra in list spc into a matrix of 1 row ==================
    spc_m <- lapply(spc,
      function(x) matrix(x, ncol = length(x), byrow = FALSE))
    # Add dimnames (wavenumbers for columns and unique_id for rows
    spc_m <- foreach::foreach(i = 1:length(spc_m),
      .final = function(i) setNames(i, names(spc_m))) %do% {
      colnames(spc_m[[i]]) <- round(wavenumbers[[i]], 1)
      rownames(spc_m[[i]]) <- unique_id
      data.table::as.data.table(spc_m[[i]])
      }

    # Save all relevant data parameters (metadata)
    # in tibble data frame (class "data.frame" and "tbl_diff" ==================
    metadata <- tibble::data_frame(
      unique_id = unique_id,
      file_id = file_name_nopath, # pb (20170514): changed `scan_id` to `file_id`
      sample_id = sample_id,
      rep_no = as.numeric(rep_no),
      date_time_sm = max(date_time),
      date_time_rf = min(date_time),
      sample_name = SNM,
      instr_name_range = instr_range,
      resolution_wn = RES,
      # Result spectrum; e.g. "AB" = Absorbance
      result_spc = ifelse(length(unique(PLF)) == 1, unique(PLF), unique(PLF)[2]),
      beamspl = BMS,
      laser_wn = LWN,
      # `spc_in_file`: character vector of spectra found in OPUS file
      spc_in_file = paste(unlist(list_assigned_t[["spc_code"]]),
        collapse = ";", sep = ";"),
      zero_filling = ZFF, # Zero filling factor for fourier transformation
      # Temperature of scanner during sample measurement
      temp_scanner_sm = TSC_all[[length(TSC_all)]], # select last element
      # Temperature of scanner during reference measurement;
      # if there is only one element in TSC_all, temperature during reference
      # mesurement is not saved
      temp_scanner_rf = ifelse(length(TSC_all) == 1, NA, TSC_all[[1]]),
      # Relative humidity
      hum_rel_sm = HUM_rel[[length(HUM_rel)]], # sample measurement
      hum_rel_rf = ifelse(length(HUM_rel) == 1, NA, HUM_rel[[1]]), # reference
      # measurement
      # Absolute humidity; sample measurement (sm); reference measurment (rf);
      # note: for Vertex 70 instrument HUA is not present, in this case,
      # HUM_abs is a list without elements
      hum_abs_sm = ifelse(length(HUM_abs) != 0, HUM_abs[[length(HUM_abs)]], NA),
      hum_abs_rf = ifelse(length(HUM_abs) == 1 | length(HUM_abs) == 0, NA,
        HUM_abs[[1]]) # reference measurement
    )

    ## Allocate and return data from spectra in output list (out) ==============
    out <- list(
      'metadata' = metadata,
      'spc' = spc_m[["spc"]],
      'spc_nocomp' = if ("spc_nocomp" %in% extract &&
        "spc_nocomp" %in% names(spc_m)) {
          spc_m[["spc_nocomp"]]} else {NULL},
      'sc_sm' = if ("ScSm" %in% extract && "ScSm" %in% names(spc_m)) {
        spc_m[["ScSm"]]} else {NULL},
      'sc_rf' = if ("ScRf" %in% extract && "ScRf" %in% names(spc_m)) {
	      spc_m[["ScRf"]]} else {NULL},
      'ig_sm' = if ("IgSm" %in% extract && "IgSm" %in% names(spc_m)) {
	      spc_m[["IgSm"]]} else {NULL},
      'ig_rf' = if ("IgRf" %in% extract && "IgRf" %in% names(spc_m)) {
	      spc_m[["IgRf"]]} else {NULL},
      # Wavenumbers of final AB spectra
      wavenumbers = wavenumbers[["spc"]],
      wavenumbers_sc_sm = if ("ScSm" %in% extract) {
        wavenumbers[["ScSm"]]} else {NULL},
      wavenumbers_sc_rf = if ("ScRf" %in% extract) {
        wavenumbers[["ScRf"]]} else {NULL}
    )

    # Print message that file was read if option is set
    if (print_progress == TRUE) {
      message(
        paste0("Extracted spectra data from file: <", file_name_nopath, ">")
      )
    }
  # Return spectra data and metadata contained as elements in list out
  out
  }) # closes try() function

}

#' @title Read a list of Bruker OPUS spectrum binary files.
#' @description
#' Read multiple spectral files measured with a Bruker FTIR Instrument. Files
#' containing spectra are in OPUS binary format.
#' \code{read_opus_univ} is a wrapper for \code{read_opus_bin_univ()})
#' @param fnames List of character vectors containing full path names of spectra
#' @param extract Character vector of spectra types to extract from file.
#' Possible values are: "spc" (AB block in Bruker Opus software), "spc_nocomp"
#' (Spectra before final atmospheric compensation; only present if background
#' correction has been set in Opus), "ScSm" (Single channel spectrum of the
#' sample), "ScRf" (Single channel spectrum of the sample), "IgSm" (Interferogram
#' of the sample), "IgRf" (Interferogram of the reference). Default is
#'  \code{extract = c("spc")}.
#' @param parallel Logical (\code{TRUE} or \code{FALSE} indicating whether
#' files are read in parallel (multiple processors or multiple cores)).
#' Default is \code{parallel = FALSE}. If \code{TRUE} a parallel backend needs
#' to be registered, e.g. by using the \code{doParallel} package.
#' @param atm_comp_minus4offset Logical whether spectra after atmospheric
#' compensation are read with an offset of \code{-4} bites from Bruker OPUS
#' files. Default is \code{FALSE}.
#' @usage read_opus_univ(fnames, extract = c("spc"), parallel = FALSE,
#' atm_comp_minus4offset = FALSE)
#' @return out List spectra and metadata (parameters) extracted from
#' Bruker OPUS spectrometer files. List names are the names of the OPUS
#' files whose spectral data were extracted.
#' @export
read_opus_univ <- function(fnames, extract = c("spc"), parallel = FALSE,
                           atm_comp_minus4offset = FALSE) {

  # Avoid `R CMD check` NOTE: ``no visible binding for variable ...
  i <- NULL

  if (parallel == TRUE) {
    foreach::foreach(i = 1:length(fnames),
      .export = "read_opus_bin_univ",
      .errorhandling = "pass", # error object generated by task evaluation will
      # be included with the rest of the results
      # export the foreach package to the individual workers
      .packages = c("foreach", "simplerspec"),
      .final = function(i) setNames(i, sub(".+/(.+)", "\\1", fnames))) %dopar% {
        try(
          read_opus_bin_univ(file_path = fnames[[i]], extract = extract,
            atm_comp_minus4offset = atm_comp_minus4offset)
        )
    }
  } else if (parallel == FALSE) {
    spc_list <- lapply(fnames,
      function(x) try(read_opus_bin_univ(file_path = x, extract = extract,
        atm_comp_minus4offset = atm_comp_minus4offset)))
    names(spc_list) <- sub(".+/(.+)", "\\1", fnames)
    spc_list
  }
}
