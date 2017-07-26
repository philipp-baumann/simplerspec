#' @title Predict soil properties of new spectra based on a list of calibration models
#' @description Append predictions for a set of responses specified by a list
#' of calibration models and a tibble containing preprocessed spectra as
#' list-columns
#' @param model_list List of model output generated from calibration step
#' (\code{pls_ken_stone()}
#' @param spc_tbl Tibble of spectra after preprocessing
#' (\code{preprocess_spc()})
#' @param slice Logical expression wheather only one row per sample_id returned.
#' @usage predict_from_spc(model_list, spc_tbl, slice = TRUE)
#' @return tibble with new columns \code{model}, and predicted values with
#' column names of model list.
#' @export
predict_from_spc <- function(model_list, spc_tbl, slice = TRUE) {

  # Helper function that extracts all pls_model elements (outputs from caret)
  # for a list of models
  models <- lapply(model_list, function(x) x[["model"]])

  # Group by spectra tibble by sample_id and keep one row per sample_id
  if (slice == TRUE) {
    spc_tbl <- spc_tbl %>% dplyr::group_by(sample_id) %>%
      slice(1L)
  }

  # Collect preprocessed spectra in one data.table
  spc <- data.table::rbindlist(spc_tbl$spc_pre)
  pred_caret <- caret::extractPrediction(
    models,
    unkX = spc
  )

  # Number of caret model objects used to predict
  n <- length(unique(pred_caret$object))
  # Add sample_id from metadata of spectra to predicted values
  sample_id <- spc_tbl$sample_id
  # id column to long form data frame
  id <- rep(sample_id, n)
  pred_id <- cbind(pred_caret, sample_id = id)
  # Get data into wide form
  pred_wide <- tidyr::spread(
    data = pred_id, key = "object", value = "pred"
  )

  # Join predictions with tibble
  dplyr::inner_join(spc_tbl, pred_wide)
}
