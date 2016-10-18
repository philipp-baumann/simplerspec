#' @title Predict soil properties of new spectra based on calibration models
#' @description
#' Function that uses pre-processed spectra, additional metadata of new
#' samples, and caret model output for the different soil property models
#' to create predicted values.
#' @param model_list List that contains caret output objects
#' of the different calibration models to predict (one model per soil property)
#' @param spectra_list List that contains spectra and additional data
#' after pre-processing (\code{do_pretreatment()}including metadata
#' (\code{sample_ID})
#' @usage predict_from_spectra(model_list, spectra_list)
#' @export
predict_from_spectra <- function(model_list, spectra_list) {

  # Use extractPrediction function (caret) and supply model_list that contains
  # caret calibration outputs; use pre-processed spectra dataset (list
  # resulting from do_pretreatment())
  predictions_caret <- caret::extractPrediction(
    models_prediction,
    unkX = model_list$MIR0
  )

  # Convert data.frame into long form; one sample should be represented by
  # one single row and the predicted values of soil properties should be
  # in the different columns
  # Use the tidyr::spread() function (from tidyr packge)
  # to gather columns into rows
  # Add sample_ID column to uniquely identify observations

  # Number of caret model objects used to predict
  n <- length(unique(predictions_caret$object))
  # Add sample_ID from metadata of spectra to predicted values
  sample_ID <- spectra_list$data_meta$ID
  # Repeat meta_data for each of the additional model rows and add
  # ID column to long form data frame
  id <- rep(sample_ID, n)
  predictions_metadata <- cbind(predictions_caret, sample_ID = id)
  # Get data into wide form
  predictions_wide <- tidyr::spread(
    data = predictions_metadata, key = "object", value = "pred"
  )

  # Join predictions with tibble
  dplyr::inner_join(spc_tbl, pred_wide)
}
