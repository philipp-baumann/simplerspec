### VIP.R: Implementation of VIP (variable importance in projection)(*) for the
### `pls' package.
### $Id: VIP.R,v 1.2 2007/07/30 09:17:36 bhm Exp $

### Copyright ? 2006,2007 Björn-Helge Mevik
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
### Bj?rn-Helge Mevik
### bhx6@mevik.net
### R?dtvetvien 20
### N-0955 Oslo
### Norway

### (*) As described in Chong, Il-Gyo & Jun, Chi-Hyuck, 2005, Performance of
### some variable selection methods when multicollinearity is present,
### Chemometrics and Intelligent Laboratory Systems 78, 103--112.

## VIP returns all VIP values for all variables and all number of components,
## as a ncomp x nvars matrix.
#' @export
VIP <- function(object) {
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

#' @export
extract_pls_vip <- function(mout) {
  # Compute VIP for all wavenumbers and select only VIPs with ncomp in final
  # model
  final_model <- mout$pls_model$finalModel
  vip <- VIP(object = final_model)[final_model$ncomp, ]
  # Collect wavenumbers from preprocessed spectra
  wn <- as.numeric(colnames(mout$data$calibration$spc_pre[[1]]))
  # Create a data frame with wavenumbers and VIP scores
  tibble::data_frame_(lazyeval::lazy_dots(wavenumber = wn, vip = vip))
}


