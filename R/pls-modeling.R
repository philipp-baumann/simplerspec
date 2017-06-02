# Perform calibration sampling based on spectral PCA ------------
#' @title Split spectra into calibration and validation sets
#' @description Perform calibration sampling based on
#' the Kennard-Stones algorithm.
#' @param spec_chem data.frame that contains chemical
#' and IR spectroscopy data
#' @param ratio_val Ratio of number of validation and all samples.
#' @param pc Number of principal components (numeric)
#' @param print logical expression weather calibration
#' @param validation Logical expression weather
#' calibration sampling is performed
#' (\code{TRUE} or \code{FALSE}).
#' @usage ken_stone(spec_chem, ratio_val, pc, print = TRUE,
#' validation = TRUE)
#' @export
ken_stone_q <- function(spec_chem, ratio_val, split_method, pc = 2,
  print = TRUE,
  validation = TRUE, evaluation_method = "test_set",
  invert = FALSE, env = parent.frame()) {
  MIR <- model <- type <- PC1 <- PC2 <- NULL

  # Evaluate the invert argument in the parent function (pls_ken_stone)
  invert <- eval(invert, envir = parent.frame())
  # Evaluate the validation argument in the parent function (pls_ken_stone)
  # 20170602: revise argument name and values of validation;
  # Suggestion: replace validation = TRUE or FALSE with
  # new argument evaluation_method = “test_set” or “resampling”
  if (!missing(validation)) {
    evaluation_method <- validation
  }
  evaluation_method <- eval(evaluation_method, envir = parent.frame())

  # Slice based on sample_id if spectral data is in tibble class
  if(tibble::is_tibble(spec_chem)) {
    spec_chem <- spec_chem %>% group_by(sample_id) %>%
      slice(1L)
  }

  if(evaluation_method == "test_set") {
    # pc = 0.99 before !!!
    pc_number <- eval(pc, envir = parent.frame())

    if(invert == FALSE) {
    ## Select calibration set by Kennard-Stones algorithm
    # Check if tibble; if yes slice tibble and bind list of data.tables in
    # one data table for spectral data
      if(tibble::is_tibble(spec_chem)) {
        spc_pre <- as.matrix(data.table::rbindlist(spec_chem$spc_pre))
        # k = number of samples to select
        # pc = if provided, the number of principal components
        # (see ?kenStone)
        sel <- prospectr::kenStone(X = spc_pre,
          k = round((1 - ratio_val) * nrow(spec_chem)), pc = substitute(pc_number))
      } else {
      sel <- prospectr::kenStone(X = spec_chem$MIR,
        k = round((1 - ratio_val) * nrow(spec_chem)), pc = substitute(pc_number))
      }
      # Split MIR data into calibration and validation set using
      # the results of Kennard-Stone Calibration Sampling
      # Selct by row index of calibration samples
      val_set <- spec_chem[- sel$model, ]
      cal_set <- spec_chem[sel$model, ]

      # Optionally split up calibation (train) and validation (test) sets
      # randomly; use function from modelr package
      # !!! Important note: The option to split up the calibration and
      # sets randomly is still experimental and a modification for the
      # PC space projection is not yet implemented for graphical output.
      # p_pc ggplot2 output needs to be updated for split_method = "random"
      if(split_method == "random") {
        # Split data sets into test and traing using modelr package
        df_split <- modelr::crossv_mc(spec_chem, n = 1, test = ratio_val)
        # Select train of df_split and convert back into tibble,
        # assign to calibration set
        cal_set <-  df_split[1, ] %>% .$train %>% .[[1]] %>% as_tibble()
        # Select test of df_split and convert back into tibble,
        # assign to validation set
        val_set <- df_split[1, ] %>% .$test %>% .[[1]] %>% as_tibble()
      }
      sel_df_cal <- data.frame(sel$pc[sel$model, 1:2])
      sel_df_val <- data.frame(sel$pc[- sel$model, 1:2])
    } else {

        if(tibble::is_tibble(spec_chem)) {
          spc_pre <- as.matrix(data.table::rbindlist(spec_chem$spc_pre))
          sel <- prospectr::kenStone(X = spc_pre,
            k = round(ratio_val * nrow(spec_chem)), pc = substitute(pc_number))
        } else {
        ## Select validation set by Kennard-Stones algorithm
        sel <- prospectr::kenStone(X = spec_chem$MIR,
          k = round(ratio_val * nrow(spec_chem)), pc = substitute(pc_number))
        }
    sel_df_cal <- data.frame(sel$pc[- sel$model, 1:2])
    sel_df_val <- data.frame(sel$pc[sel$model, 1:2])
    # Split MIR data into calibration and validation set using
    # the results of Kennard-Stone Calibration Sampling
    # Selct by row index of calibration samples
    val_set <- spec_chem[sel$model, ]
    cal_set <- spec_chem[- sel$model, ]
    }

    # Add additional columns to calibration set and validation sets for plotting
    sel_df_cal$type <- as.factor(
      rep("calibration", nrow(sel_df_cal))
    )
    sel_df_val$type <- as.factor(
      rep("validation", nrow(sel_df_val)))
    sel_df <- rbind(sel_df_cal, sel_df_val)
    # Compute ratio needed to make the figure square
    ratio <- with(sel_df, diff(range(PC1))/diff(range(PC2)))
    # Save graph showing the selected calibration and validation samples
    # for the first two principal components (pc)
    p_pc <- ggplot2::ggplot(data = sel_df) +
      ggplot2::geom_point(
        ggplot2::aes(x = PC1, y = PC2, shape = type), size = 4) +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::scale_shape_manual(values = c(1, 19)) +
      ggplot2::scale_colour_manual(values = c("black", "red")) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.title = ggplot2::element_blank())

    # Print outputs to list
    list_out <- list(
      calibration = cal_set,
      validation = val_set,
      p_pc = p_pc
    )
    list_out
    # Check number of observations (rows) for calibration set
    # nrow(cal_set)
  } else {
    cal_set <- spec_chem
    list(
      calibration = cal_set,
      validation = NULL,
      p_pc = NULL
    )
  }
}

#' @title Perform model tuning
#' @description Uses function from caret to to model tuning
#' for PLS regression.
#' @param x list from calibration sampling
#' @param variable response variable for PLS regression, supplied
#' as character expression
#' @param validation Logical expression weather an independent
#' validation is performed.
#' @param env Environment where function is evaluated
#' @export
tune_model_q <- function(x, variable,
  env = parent.frame(), validation = TRUE, evaluation_method = "test_set") {
  calibration <- NULL
  # List of calibration and validation samples
  # set up a cross-validation scheme
  # create 10 folds that we will keep for the different
  # modeling approaches to allow comparison
  # randomly break the data into 10 partitions
  # note that k is the total number of samples for leave-one-out
  # use substitute function to make non-standard evaluation
  # of variable argument (looks at a function as argument,
  # sees code used to compute value;
  # see chapter 13.1 Capturing expressions
  # in Advanced R (Hadley Wickham)
  # !! p. 270
  r <- eval(variable, x$calibration, env)
  idx <- caret::createFolds(y = r, k = 10, returnTrain = T) # update ***
  idx
  # inject the index in the trainControl object
  tr_control <- caret::trainControl(method = "cv", index = idx,
    savePredictions = T)
  if (!missing(validation)) {
    evaluation_method <- validation
  }
  if (evaluation_method == "test_set") {
  tr_control
  } else {
  tr_control
  }
}

## Adapt model tuning to leave-one-out cross-validation ========================

#' @title Perform model tuning
#' @description Uses function from caret to to model tuning
#' for PLS regression.
#' @param x list from calibration sampling
#' @param variable response variable for PLS regression, supplied
#' as character expression
#' @param validation Logical expression weather an independent
#' validation is performed.
#' @param env Environment where function is evaluated
#' @export
tune_model_loocv_q <- function(x, variable,
  env = parent.frame(),
  validation = TRUE, evaluation_method = "test_set") {
  calibration <- NULL
  # List of calibration and validation samples
  # set up a cross-validation scheme
  # create 10 folds that we will keep for the different
  # modeling approaches to allow comparison
  # randomly break the data into 10 partitions
  # note that k is the total number of samples for leave-one-out
  # use substitute function to make non-standard evaluation
  # of variable argument (looks at a function as argument,
  # sees code used to compute value;
  # see chapter 13.1 Capturing expressions
  # in Advanced R (Hadley Wickham)
  # !! p. 270
  r <- eval(variable, x$calibration, env)
  # idx <- caret::createFolds(y = r, k = 10, returnTrain = T) # update ***
  # idx
  # inject the index in the trainControl object
  tr_control <- caret::trainControl(method = "LOOCV", # index = idx,
    savePredictions = T)
  if (!missing(validation)) {
    evaluation_method <- validation
  }
  if (evaluation_method == "test_set") {
    tr_control
  } else {
    tr_control
  }
}

## Adapt model tuning to repeated k-fold cross-validation ======================

#' @title Perform model tuning
#' @description Uses function from caret to to model tuning
#' for PLS regression.
#' @param x list from calibration sampling
#' @param variable response variable for PLS regression, supplied
#' as character expression
#' @param validation Logical expression weather an independent
#' validation is performed.
#' @param env Environment where function is evaluated
#' @export
tune_model_rcv_q <- function(x, variable,
  env = parent.frame(),
  validation = TRUE, evaluation_method = "test_set") {
  calibration <- NULL
  # List of calibration and validation samples
  # set up a cross-validation scheme
  # create 10 folds that we will keep for the different
  # modeling approaches to allow comparison
  # randomly break the data into 10 partitions
  # note that k is the total number of samples for leave-one-out
  # use substitute function to make non-standard evaluation
  # of variable argument (looks at a function as argument,
  # sees code used to compute value;
  # see chapter 13.1 Capturing expressions
  # in Advanced R (Hadley Wickham)
  # !! p. 270
  r <- eval(variable, x$calibration, env)
  idx <- caret::createMultiFolds(y = r, k = 10, times = 5) # update ***
  # inject the index in the trainControl object
  tr_control <- caret::trainControl(method = "repeatedcv", index = idx,
    savePredictions = T)
  if (!missing(validation)) {
    evaluation_method <- validation
  }
  if (evaluation_method == "test_set") {
    tr_control
  } else {
    tr_control
  }
}

## Fitting models without parameter tuning =====================================
# 5.9; https://topepo.github.io/caret/model-training-and-tuning.html

#' @title Perform model fitting without parameter tuning
#' @description Uses function from caret to set model tuning to none
#' for PLS regression.
#' @param x list from calibration sampling
#' @param variable response variable for PLS regression, supplied
#' as character expression
#' @param validation Logical expression weather an independent
#' validation is performed.
#' @param env Environment where function is evaluated
#' @export
tune_model_none_q <- function(x, variable,
  env = parent.frame(),
  validation = TRUE, evaluation_method = "test_set") {
  calibration <- NULL
  # List of calibration and validation samples
  # set up a cross-validation scheme
  # create 10 folds that we will keep for the different
  # modeling approaches to allow comparison
  # randomly break the data into 10 partitions
  # note that k is the total number of samples for leave-one-out
  # use substitute function to make non-standard evaluation
  # of variable argument (looks at a function as argument,
  # sees code used to compute value;
  # see chapter 13.1 Capturing expressions
  # in Advanced R (Hadley Wickham)
  # !! p. 270
  r <- eval(variable, x$calibration, env)
  idx <- caret::createFolds(y = r, k = 10, returnTrain = T) # update ***
  # inject the index in the trainControl object
  tr_control <- caret::trainControl(method = "none", index = idx,
    savePredictions = T)
  if (!missing(validation)) {
    evaluation_method <- validation
  }
  if (evaluation_method == "test_set") {
    tr_control
  } else {
    tr_control
  }
}

#' @title Perform model tuning
#' @description Uses function from caret to to model tuning
#' for PLS regression.
#' @param x list from calibration sampling
#' @param variable response variable for PLS regression, supplied
#' as character expression
#' @param validation Logical expression weather an independent
#' validation is performed.
#' @param env Environment where function is evaluated
#' @export
tune_model <- function(x, variable,
  env = parent.frame(),
  validation = TRUE, evaluation_method = "test_set") {
  tune_model_q(x, substitute(variable), env)
}

# Fit a PLS regression model using the caret package ------------

#' @title Fit a PLS regression model
#' (quoted version of the function)
#' @description Uses the caret package to perform PLS modeling.
#' Spectra are centered and scaled prior to modeling.
#' @param x List that contains calibration
#' set, validation set, and model tuning options
#' @param validation Logical expression weather independent
#' validation is performed
#' @param variable Response variable to be modeled
#' @param tr_control Object that defines controlling parameters
#' of the desired internal validation framework
#' @param env Environment where function is evaluated
#' @export
fit_pls_q <- function(x,
  validation = TRUE, evaluation_method = "test_set",
  variable, tr_control, env = parent.frame(),
  pls_ncomp_max = 20, ncomp_fixed = 5,
  center, scale, tuning_method = "resampling") {
  # Fit a partial least square regression (pls) model
  # center and scale MIR (you can try without)
  calibration <- MIR <- NULL
  v <- eval(variable, x$calibration, env)
  # ? Is it really necessary to evaluate this in the parent frame?
  pls_ncomp_max <- eval(pls_ncomp_max, envir = parent.frame())
  # Evaluate fixed number of PLS regression components
  # from ncomp_fixed object in parent frame (pls_ken_stone function)
  ncomp_fixed <- eval(ncomp_fixed, envir = parent.frame())
  # Test whether the spectral object has the class "tibble"
  if (tibble::is_tibble(x$calibration)) {
    spc_pre <- data.table::rbindlist(x$calibration$spc_pre)
    if(scale == TRUE & center == TRUE) {
      if(tuning_method == "resampling") {
        # Fit model with parameter tuning
        pls_model <- caret::train(x = spc_pre, y = v,
          method = "pls",
          tuneLength = pls_ncomp_max,
          trControl = tr_control,
          preProcess = c("center", "scale"))
      } else if (tuning_method == "none") {
        # Fit model without parameter tuning
        pls_model <- caret::train(x = spc_pre, y = v,
          method = "pls",
          trControl = tr_control,
          preProcess = c("center", "scale"),
          tuneGrid = data.frame(ncomp = ncomp_fixed))
      }
    } else {
      # No centering and scaling!
      pls_model <- caret::train(x = spc_pre, y = v,
        method = "pls",
        tuneLength = pls_ncomp_max,
        trControl = tr_control)
    }
  } else {
    # depreciated list interface
    pls_model <- caret::train(x = x$calibration$MIR, y = v,
      method = "pls",
      tuneLength = pls_ncomp_max,
      trControl = tr_control,
      preProcess = c("center", "scale")
    )
  }
  # Collect fitted object into a list
  # fitList_cal <- list(pls = fit_pls)
  # fitList_cal
  pls_model
}


#' @title Fit a PLS regression model
#' @description Uses the caret package to perform PLS modeling.
#' Spectra are centered and scaled prior to modeling.
#' @param x List that contains calibration
#' set, validation set, and model tuning options
#' @param validation Logical expression weather independent
#' validation is performed
#' @param variable Response variable to be modeled
#' @param env Environment where function is evaluated
#' @export
fit_pls <- function(x, validation = TRUE,
  variable, env = parent.frame()) {
  q(x = x, validation = TRUE, evaluation_method = "test_set",
    variable = substitute(variable), env
  )
}

# Fit a random forest model using the caret package -------------

#' @title Fit a random forest model
#' (quoted version of the function)
#' @description Uses the caret package to perform random forest
#' modeling.
#' Spectra are centered and scaled prior to modeling.
#' @param x List that contains calibration
#' set, validation set, and model tuning options
#' @param validation Logical expression weather independent
#' validation is performed
#' @param variable Response variable to be modeled
#' @param tr_control Object that defines controlling parameters
#' of the desired internal validation framework
#' @param env Environment where function is evaluated
#' @export
fit_rf_q <- function(x,
  validation = TRUE, evaluation_method = "test_set",
  variable, tr_control, ntree_max = 500, env = parent.frame()) {
  # Fit a partial least square regression (pls) model
  # center and scale MIR (you can try without)
  calibration <- MIR <- NULL
  v <- eval(variable, x$calibration, env)
  ntree_max <- eval(ntree_max, envir = parent.frame())
  if (tibble::is_tibble(x$calibration)) {
    spc_pre <- data.table::rbindlist(x$calibration$spc_pre)
    rf_model <- caret::train( x = spc_pre, y = v,
      method = "rf",
      ntree = ntree_max,
      trControl = tr_control,
      preProcess = c("center", "scale")
    )

  } else {
    rf_model <- caret::train(x = x$calibration$MIR, y = v,
      method = "rf",
      ntree = ntree_max,
      trControl = tr_control,
      preProcess = c("center", "scale")
    )
  }
  rf_model
}


# Evaluate PLS performance (validation and cross-validation) ----

#' @title Evaluate PLS performance
#' @description Calculate model performance indices based
#' on observed and predicted values of validation and calibration
#' set, and internal cross-validation
#' @param x List that contains calibration and validation data
#' frame with combined spectral and chemical data
#' @param pls_model List with PLS regression model output from
#' the caret package
#' @param variable Response variable (e.g. chemical property) to be
#' modelled (needs to be non-quoted expression). \code{variable}
#' needs to be a column name in the \code{validation} data.frame
#' (element of \code{x})
#' @param validation Logical expression if independent validation
#' is performed (split data set into calibration set and
#' validation set)
#' @param print Print observed vs. predicted for calibration
#' and validation. Default is \code{TRUE}.
#' @param env Specifiy the environment in which the function is
#' called. Default argument of \code{env} is
#' \code{parent.frame()}
#' @export
evaluate_pls_q <- function(x, pls_model, variable,
  validation, evaluation_method,
  tuning_method, print = TRUE, env = parent.frame()) {
  # Set global variables to NULL to avoid R CMD check notes
  MIR <- object <- model <- dataType <- obs <- pred <- NULL
  ncomp <- finalModel <- rmsd <- r2 <- r2 <- rpd <- n <- NULL
  rmse <- calibration <- NULL
  # Collect fitted object into a list
  list_models <- list(pls = pls_model)
  # Evaluate validation argument in parent.frame !!!
  if (!missing(validation)) {
    evaluation_method <- validation
  }
  evaluation_method <- eval(evaluation_method, envir = parent.frame())
  # Evaluate tuning_method argument in parent.frame
  tuning_method <- eval(tuning_method, envir = parent.frame())
  # Extract best tuning parameters and associated cv predictions
  if(evaluation_method == "test_set") {
    # !!! experimental: ignore predobs_cal
    # predobs_cal <- plyr::ldply(list_models,
    #  function(x) plyr::match_df(x$pred, x$bestTune),
    #  .id = "model"
    #)
    # Calculate training (calibration) and test (validation) data
    # predictions based on pls model with calibration data
    v <- eval(variable, x$validation, env)
    if (tibble::is_tibble(x$validation)) {
      spc_pre <- data.table::rbindlist(x$validation$spc_pre)
      predobs_val <- caret::extractPrediction(list_models,
        testX = spc_pre, testY = v) # update ***
      # Append sample_id column to predobs_val data.frame
      # extract sample_id from validation set
      predobs_val$sample_id <- c(
        x$calibration$sample_id, x$validation$sample_id)
    } else {
      # depreciated
      predobs_val <- caret::extractPrediction(list_models,
        testX = x$validation$MIR, testY = v) # update ***
    }
    # Create new data frame column <object>
    predobs_val$object <- predobs_val$model

    # Replace levels "Training" and "Test" in dataType column
    # by "Calibration" and "Validation" (rename levels of factor)
    predobs_val$dataType <- plyr::revalue(predobs_val$dataType,
      c("Test" = "Validation", "Training" = "Calibration")
    )
    # Change the order of rows in the data frame
    # Calibration as first level (show Calibration in ggplot graph
    # on left panel)
    predobs_val$dataType <- factor(predobs_val$dataType,
      levels = c("Calibration", "Validation"))
    # Calculate model performance indexes by model and dataType
    # uses package plyr and function summary.df of SPECmisc.R
    stats <- plyr::ddply(predobs_val, c("model", "dataType"),
      function(x) summary_df(x, "obs", "pred")
    )
  # Check whether method = "none" argument is selected in train();
  # this is the case when ncomp_fixed argument in pls_ken_stone() is
  # evaluated
  # Checking for the existence of a <pred> element in the train function output
  # list can be dangerous and doesn't work in all cases when using
  # e.g. is.element('pred', x) or is.null(x$pred);
  # Problems can occur e.g. if a list element contains NULL element;
  # see
  # http://stackoverflow.com/questions/7719741/how-to-test-if-list-element-exists
  } else if (evaluation_method == "resampling" & tuning_method == "resampling") {
    # Good discussion on which cross-validation results are returned from caret
    # Extract best tuning parameters and associated cv predictions
    # http://stats.stackexchange.com/questions/219154/how-does-cross-validation-in-train-caret-precisely-work
    # Alternative solution for one model: conformal::GetCVPreds(model) function
    # see https://github.com/cran/conformal/blob/master/R/misc.R
    predobs_cv <- plyr::ldply(list_models,
      function(x) plyr::match_df(x$pred, x$bestTune),
      .id = "model"
    )
    # Extract auto-prediction
    predobs <- caret::extractPrediction(list_models)
    # !!! new ---
    # Replace levels "Training" dataType column
    # by "Calibration" (rename levels of factor)
    predobs$dataType <- plyr::revalue(predobs$dataType,
      c("Training" = "Calibration")
    )
    # Append sample_id column to predobs data.frame
    # extract sample_id from calibration set
    predobs$sample_id <- x$calibration$sample_id
    # Create rowIndex for calibration tibble
    x$calibration$rowIndex <- 1:nrow(x$calibration)
    # Generate sample_id column for rowIndex of pred list element of
    # train object; select only rowIndex and sample_id of calibration tibble
    cal_index <- x$calibration %>% dplyr::select(rowIndex, sample_id)
    # Make a full join of cal_index and predobs_cv
    predobs_cv <- dplyr::full_join(predobs_cv, cal_index) %>%
      # average observed and predicted values by sample_id
      dplyr::group_by(sample_id) %>%
      # Average observed and predicted values
      dplyr::mutate(obs = mean(obs), pred = mean(pred), pred_sd = sd(pred)) %>%
      # slice data set
      dplyr::slice(1L)
      # optinally: state number of observations per group
      # The function n() in dlpyr gives you the number of observations
    # predobs_cv$sample_id <- x$calibration[predobs_cv$rowIndex, ]
    # ---
    predobs_cv$object <- predobs_cv$model
    predobs_cv$dataType <- "Cross-validation"
    predobs_cv <- dplyr::select(
      # !!! sample_id newly added
      predobs_cv, obs, pred, pred_sd, model, dataType, object
    )
    # Add column pred_sd to predobs data frame (assign values to NA) so that
    # column pred_sd is retained in predobs_cv after dplyr::bind_rows
    predobs$pred_sd <- NA
    # Desn't work because some columns are turned into numeric;
    # resulting data frame has only two rows
    # predobs_val <- rbind(predobs, predobs_cv) # !!! check columns
    predobs_val <- dplyr::bind_rows(predobs, predobs_cv)
    # Calculate model performance indexes by model and dataType
    # uses package plyr and function summary.df of SPECmisc.R
    stats <- plyr::ddply(predobs_val, c("model", "dataType"),
      function(x) summary_df(x, "obs", "pred")
    )
  }
  # Experimental: return predobs_cv
  # !!! Experimental: return predicted values
  # return(x$pred)
  # return(predobs)
  # Add number of components to stats; from finalModel list item
  # from train() function output (function from caret package)
  stats$ncomp <- rep(pls_model$finalModel$ncomp, nrow(stats))
  # !!! Experimental: return stats
  # return(stats)
  # Add range of observed values for validation and calibraton
  # get range from predicted vs. observed data frame
  # stored in object predobs
  obs_cal <- subset(predobs_val, dataType == "Calibration")$obs
  # Get name of predicted variable; see p. 261 of book
  # "Advanced R" (Hadley Wickham)
  variable_name <- deparse(variable)

  if(evaluation_method == "test_set") {
    # Assign validation set to separate data frame
    obs_val <- subset(predobs_val, dataType == "Validation")$obs
    # before: deparse(substitute(variable))
    df_range <- data.frame(
      variable = rep(variable_name, 2),
      dataType = c("Calibration", "Validation"),
      min_obs = c(range(obs_cal)[1], range(obs_val)[1]),
      median_obs = c(median(obs_cal), median(obs_val)),
      max_obs = c(range(obs_cal)[2], range(obs_val)[2]),
      mean_obs = c(mean(obs_cal), mean(obs_val)),
      CV = c(sd(obs_cal) / mean(obs_cal) * 100,
        sd(obs_val) / mean(obs_val) * 100)
    )
  } else if (evaluation_method == "resampling" & tuning_method == "resampling") {
    # Assign cross-validation set to separate data frame
    obs_val <- subset(predobs_val, dataType == "Cross-validation")$obs
    df_range <- data.frame(
      variable = rep(variable_name, 2),
      dataType = c("Calibration", "Cross-validation"),
      min_obs = c(range(obs_cal)[1], range(obs_val)[1]),
      median_obs = c(median(obs_cal), median(obs_val)),
      max_obs = c(range(obs_cal)[2], range(obs_val)[2]),
      mean_obs = c(mean(obs_cal), mean(obs_val)),
      CV = c(sd(obs_cal) / mean(obs_cal) * 100,
        sd(obs_val) / mean(obs_val) * 100)
    )
  }

  # Join stats with range data frame (df_range)
  stats <- plyr::join(stats, df_range, type = "inner")
  annotation <- plyr::mutate(stats,
    rmse = as.character(as.expression(paste0("RMSE == ",
      round(rmsd, 2)))),
    r2 = as.character(as.expression(paste0("italic(R)^2 == ",
      round(r2, 2)))),
    rpd = as.character(as.expression(paste("RPD == ",
      round(rpd, 2)))),
    n = as.character(as.expression(paste0("italic(n) == ", n))),
    ncomp = as.character(as.expression(paste0("ncomp == ",
      ncomp)))
  )

  # Plot predicted vs. observed values and model indexes
  # update label, xlim, and ylim ***
  # Add label number of samples to facet_grid using a
  # labeling function
  # ! Update labeller API:
  # https://github.com/hadley/ggplot2/commit/ef33dc7
  # http://sahirbhatnagar.com/facet_wrap_labels

  # Prepare lookup character vector
  make_label <- function(x, validation = TRUE, evaluation_method = "test_set") {
    dataType <- n <- NULL
    if (!missing(validation)) {
      evaluation_method <- validation
    }
    if (evaluation_method == "test_set") {
      c(`Calibration` = paste0("Calibration", "~(",
        x[x$dataType == "Calibration", ]$n, ")"
      ),
        `Validation` = paste0("Validation", "~(",
          x[x$dataType == "Validation", ]$n, ")"
        )
      )
    } else {
      c(`Calibration` = paste0("Calibration", "~(",
        x[x$dataType == "Calibration", ]$n, ")"
      ),
        `Cross-validation` = paste0("Cross-validation", "~(",
          x[x$dataType == "Cross-validation", ]$n, ")"
        )
      )
    }
  }
  if (evaluation_method == "test_set") {
    label_validation <- make_label(x = annotation,
      evaluation_method = "test_set"
    )
  } else {
    label_validation <- make_label(x = annotation,
      evaluation_method = "resampling"
    )
  }

  # Rename labels on the fly with a lookup character vector
  to_string <- ggplot2::as_labeller(
    x = label_validation, ggplot2::label_parsed
  )

  # -------------------------------------------------------------


  # http://docs.ggplot2.org/0.9.3.1/label_parsed.html
  # some other info: https://coderclub.b.uib.no/tag/plotmath/
  # !!! now depreciated in ggplot2 >= 2.0.0
  # dataType_labeller <- function(variable, value){
  #   new <- paste0(dataType_names[value], "~(", annotation$n, ")")
  #   plyr::llply(as.character(new), function(x) parse(text = x))
  # }
  p_pred_obs <- ggplot2::ggplot(data = predobs_val) +
    ggplot2::geom_point(ggplot2::aes(x = obs, y = pred),
      shape = 1, size = 4) +
    ggplot2::geom_text(data = annotation,
      ggplot2::aes(x = -Inf, y = Inf, label = r2), size = 7,
      hjust = -0.1, vjust = 1.5, parse = TRUE) +
    ggplot2::geom_text(data = annotation,
      ggplot2::aes(x = -Inf, y = Inf, label = rmse), size = 7,
      hjust = -0.075, vjust = 4.25, parse = TRUE) +
    ggplot2::geom_text(data = annotation,
      ggplot2::aes(x = -Inf, y = Inf, label = rpd), size = 7,
      hjust = -0.1, vjust = 6.5, parse = TRUE) +
    ggplot2::facet_grid(~ dataType,
      labeller = ggplot2::as_labeller(to_string)) +
    # ggplot2::facet_grid(~ dataType,
    #   labeller = dataType_labeller) +
    ggplot2::theme_bw() +
    ggplot2::geom_abline(col = "red") +
    ggplot2::labs(x = "Observed", y = "Predicted") +
    ggplot2::xlim(c(min(predobs_val$obs) -
        0.05 * diff(range(predobs_val$obs)),
      max(predobs_val$obs) +
        0.05 * diff(range(predobs_val$obs)))) +
    ggplot2::ylim(c(min(predobs_val$obs) -
        0.05 * diff(range(predobs_val$obs)),
      max(predobs_val$obs) +
        0.05 * diff(range(predobs_val$obs)))) # +
    # theme.user

  ## ggplot graph for model comparison
  ## (arranged later in panels)
  x_label <- paste0("Observed ",
    as.character(variable_name))
  y_label <- paste0("Predicted ",
    as.character(variable_name))

  # Create model evaluation plot -----------------------------------------------
  if(pls_model$method == "pls") {
  p_model <- ggplot2::ggplot(data = predobs_val) +
    ggplot2::geom_point(ggplot2::aes(x = obs, y = pred),
      shape = 1, size = 2, alpha = 1/2) +
    ggplot2::geom_text(data = annotation,
      ggplot2::aes(x = Inf, y = -Inf, label = ncomp), size = 3,
      hjust = 1.15, vjust = -4.5, parse = TRUE) + # !!! additional label
    ggplot2::geom_text(data = annotation,
      ggplot2::aes(x = Inf, y = -Inf, label = r2), size = 3,
      hjust = 1.15, vjust = -3, parse = TRUE) +
    ggplot2::geom_text(data = annotation,
      ggplot2::aes(x = Inf, y = -Inf, label = rmse), size = 3,
      hjust = 1.12, vjust = -2.5, parse = TRUE) +
    ggplot2::geom_text(data = annotation,
      ggplot2::aes(x = Inf, y = -Inf, label = rpd), size = 3,
      hjust = 1.15, vjust = -1.25, parse = TRUE) +
    ggplot2::facet_grid(~ dataType,
      labeller = ggplot2::as_labeller(to_string)) +
    # ggplot2::facet_grid(~ dataType,
    #   labeller = dataType_labeller) +
    ggplot2::theme_bw() +
    ggplot2::geom_abline(col = "red") +
    ggplot2::labs(x = x_label, y = y_label) +
    ggplot2::xlim(c(min(predobs_val$obs) -
        0.05 * diff(range(predobs_val$obs)),
      max(predobs_val$obs) +
        0.05 * diff(range(predobs_val$obs)))) +
    ggplot2::ylim(c(min(predobs_val$obs) -
        0.05 * diff(range(predobs_val$obs)),
      max(predobs_val$obs) +
        0.05 * diff(range(predobs_val$obs)))) +
    ggplot2::coord_fixed()
  } else {

  p_model <- ggplot2::ggplot(data = predobs_val) +
    ggplot2::geom_point(ggplot2::aes(x = obs, y = pred),
      shape = 1, size = 2, alpha = 1/2) + # without ncomp label
    ggplot2::geom_text(data = annotation,
      ggplot2::aes(x = Inf, y = -Inf, label = r2), size = 3,
      hjust = 1.15, vjust = -3, parse = TRUE) +
    ggplot2::geom_text(data = annotation,
      ggplot2::aes(x = Inf, y = -Inf, label = rmse), size = 3,
      hjust = 1.12, vjust = -2.5, parse = TRUE) +
    ggplot2::geom_text(data = annotation,
      ggplot2::aes(x = Inf, y = -Inf, label = rpd), size = 3,
      hjust = 1.15, vjust = -1.25, parse = TRUE) +
    ggplot2::facet_grid(~ dataType,
      labeller = ggplot2::as_labeller(to_string)) +
    # ggplot2::facet_grid(~ dataType,
    #   labeller = dataType_labeller) +
    ggplot2::theme_bw() +
    ggplot2::geom_abline(col = "red") +
    ggplot2::labs(x = x_label, y = y_label) +
    ggplot2::xlim(c(min(predobs_val$obs) -
        0.05 * diff(range(predobs_val$obs)),
      max(predobs_val$obs) +
        0.05 * diff(range(predobs_val$obs)))) +
    ggplot2::ylim(c(min(predobs_val$obs) -
        0.05 * diff(range(predobs_val$obs)),
      max(predobs_val$obs) +
        0.05 * diff(range(predobs_val$obs)))) +
    ggplot2::coord_fixed()
  }

  if(print == TRUE) {
    print(p_model)
  }

  list(stats = stats, p_model = p_model, predobs_val = predobs_val)
}


## PLS regression modeling in one function ======================

#' @title Calibration sampling, model tuning, and PLS regression
#' @description Perform calibration sampling and use selected
#' calibration set for model tuning
#' @param spec_chem data.frame that contains IR spectroscopy
#' and chemical data
#' @param k Number of validation samples
#' @param pc Number of Principal Components used for Calibration
#' sampling (Kennard-Stones algorithm)
#' @param ratio_val Ratio of number of validation and all samples.
#' @param print Logical expression weather graphs shall be printed
#' @param validation Logical expression weather independent
#' validation is performed
#' @param variable Response variable (without quotes)
#' @param env Environment where function is evaluated
#' @export
# Note: check non standard evaluation, argument passing...
pls_ken_stone <- function(spec_chem, split_method = "ken_stone",
  ratio_val, pc = 2,
  print = TRUE, variable,
  validation = TRUE, evaluation_method = "test_set",
  invert = TRUE,
  center = TRUE, scale = TRUE,
  env = parent.frame(), pls_ncomp_max = 20, ncomp_fixed = 5,
  cv = "kfold_cv", resampling_method = "kfold_cv",
  tuning_method = "resampling") {
  calibration <- 0

  # Perform calibration sampling
  if (!missing(validation)) {
    warning("argument validation is deprecated; please use evaluation_method instead.",
      call. = FALSE)
    evaluation_method <- validation
  }
  list_sampled <- ken_stone_q(
    spec_chem, split_method, ratio_val = ratio_val, pc = substitute(pc),
    evaluation_method = substitute(evaluation_method),
    invert = substitute(invert)
  )
  # Check on method for cross-validation to be used in caret model tuning ------
  if(cv == "LOOCV") {
    # leave-one-out cross-validation
    tr_control <- tune_model_loocv_q(list_sampled,
      substitute(variable), env)
  } else if (cv == "repeatedcv") {
    # repeated k-fold cross-validation
    tr_control <- tune_model_rcv_q(list_sampled,
      substitute(variable), env)
  } else if (cv == "none") {
    # no resampling; calls caret::train(..., method = "none");
    # fixed number of PLS components; tuning_method argument has also
    # to be set to "none"
    tr_control <- tune_model_none_q(list_sampled,
      substitute(variable), env)
  } else {
    # k-fold cross validation
    tr_control <- tune_model_q(list_sampled,
      substitute(variable), env)
  }
  # Fit a pls calibration model; pls object is output from caret::train()
  # and has class train
  if(tuning_method == "resampling") {
    pls <- fit_pls_q(x = list_sampled,
      validation = TRUE, evaluation_method = "test_set",
      variable = substitute(variable), tr_control = tr_control,
      center = center, scale = scale,
      pls_ncomp_max = substitute(pls_ncomp_max), env
      )
  } else if (tuning_method == "none") {
    pls <- fit_pls_q(x = list_sampled,
      validation = TRUE, evaluation_method = "test_set",
      variable = substitute(variable), tr_control = tr_control,
      center = center, scale = scale, tuning_method = "none",
      ncomp_fixed = substitute(ncomp_fixed), env
      )
  }
  # Evaluate model accuracy (predicted vs. observed)
  stats <- evaluate_pls_q(x = list_sampled, pls_model = pls,
    variable = substitute(variable),
    validation = substitute(validation),
    evaluation_method = substitute(evaluation_method),
    tuning_method = substitute(tuning_method),
    env = parent.frame()
  )
  list(data = list_sampled, p_pc = list_sampled$p_pc,
    pls_model = pls, stats = stats$stats, p_model = stats$p_model,
    predobs_val = stats$predobs_val)
}

## Random forest modeling in one function =======================

#' @title Calibration sampling, model tuning, and random forest modeling
#' @description Perform calibration sampling and use selected
#' calibration set for model tuning
#' @param spec_chem data.frame that contains IR spectroscopy
#' and chemical data
#' @param k Number of validation samples
#' @param pc Number of Principal Components used for Calibration
#' sampling (Kennard-Stones algorithm)
#' @param ratio_val Ratio of number of validation and all samples.
#' @param print Logical expression weather graphs shall be printed
#' @param validation Logical expression weather independent
#' validation is performed
#' @param variable Response variable (without quotes)
#' @param env Environment where function is evaluated
#' @export
# Note: check non standard evaluation, argument passing...
rf_ken_stone <- function(spec_chem, split_method = "ken_stone", ratio_val,
    pc = 2, print = TRUE,
    validation = TRUE, evaluation_method = "test_set",
    variable,
    tuning_method = "resampling", ntree_max = 500,
    env = parent.frame()) {
  calibration <- 0
  # Calibration sampling
  list_sampled <- ken_stone_q(
    spec_chem, split_method, ratio_val = ratio_val, pc = substitute(pc),
    validation = TRUE
  )
  tr_control <- tune_model_q(list_sampled,
    substitute(variable), env
  )
  rf <- fit_rf_q(x = list_sampled,
    validation = TRUE, evaluation_method = "test_set",
    variable = substitute(variable), tr_control = tr_control, env,
    ntree_max = substitute(ntree_max)
  )
  stats <- evaluate_pls_q(x = list_sampled, pls_model = rf,
    variable = substitute(variable), validation = substitute(validation),
    tuning_method = substitute(tuning_method),
    env = parent.frame()
  )
  list(data = list_sampled, p_pc = list_sampled$p_pc,
    rf_model = rf, stats = stats$stats, p_model = stats$p_model)
}
