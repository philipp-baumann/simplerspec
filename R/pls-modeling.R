## Perform calibration sampling based on spectral PCA
## or random split -------------------------------------------------------------
#' @importFrom magrittr %>%
split_data_q <- function(
  spec_chem,
  split_method,
  evaluation_method = "test_set",
  ratio_val,
  ken_sto_pc = 2,
  print = TRUE,
  invert = FALSE, env = parent.frame()) {
  MIR <- model <- type <- PC1 <- PC2 <- NULL

  # Evaluate the invert argument in the parent function (fit_pls)
  invert <- eval(invert, envir = parent.frame())
  # Evaluate the validation argument in the parent function (fit_pls)
  evaluation_method <- eval(evaluation_method, envir = parent.frame())

  # Slice based on sample_id if spectral data is in tibble class
  if (tibble::is_tibble(spec_chem)) {
    spec_chem <- spec_chem %>%
      dplyr::group_by(!!rlang::sym("sample_id")) %>%
      dplyr::slice(1L)
  }

  if (evaluation_method == "test_set") {
    # pc = 0.99 before !!!
    ken_sto_pc <- eval(ken_sto_pc, envir = parent.frame())

    if (invert == FALSE) {
    ## Select calibration set by Kennard-Stones algorithm
    # Check if tibble; if yes slice tibble and bind list of data.tables in
    # one data table for spectral data
      if (tibble::is_tibble(spec_chem)) {
        spc_pre <- as.matrix(data.table::rbindlist(spec_chem$spc_pre))
        # k = number of samples to select
        # ken_sto_pc = if provided, the number of principal components
        # (see ?kenStone)
        sel <- prospectr::kenStone(X = spc_pre,
          k = round((1 - ratio_val) * nrow(spec_chem)),
          pc = substitute(ken_sto_pc))
      } else {
      sel <- prospectr::kenStone(X = spec_chem$MIR,
        k = round((1 - ratio_val) * nrow(spec_chem)),
        pc = substitute(ken_sto_pc))
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
      if (split_method == "random") {
        # Split data sets into test and traing using modelr package
        df_split <- modelr::crossv_mc(spec_chem, n = 1, test = ratio_val)
        # Select train of df_split and convert back into tibble,
        # assign to calibration set
        cal_set <- tibble::as_tibble(df_split[1, ][["train"]][[1]])
        # Select test of df_split and convert back into tibble,
        # assign to validation set
        val_set <- tibble::as_tibble(df_split[1, ][["test"]][[1]])
      }
      sel_df_cal <- data.frame(sel$pc[sel$model, 1:2])
      sel_df_val <- data.frame(sel$pc[- sel$model, 1:2])
    } else {

        if (tibble::is_tibble(spec_chem)) {
          spc_pre <- as.matrix(data.table::rbindlist(spec_chem$spc_pre))
          sel <- prospectr::kenStone(X = spc_pre,
            k = round(ratio_val * nrow(spec_chem)),
            pc = substitute(ken_sto_pc))
        } else {
        ## Select validation set by Kennard-Stones algorithm
        sel <- prospectr::kenStone(X = spec_chem$MIR,
          k = round(ratio_val * nrow(spec_chem)),
          pc = substitute(ken_sto_pc))
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
    ratio <- with(sel_df, diff(range(PC1)) / diff(range(PC2)))
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

# trainControl generating helper function
control_train_q <- function(x, response, resampling_seed,
                            env = parent.frame()) {
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
  response <- eval(response, x$calibration, env)
  # Set seed for creating resampling indices
  set.seed(eval(resampling_seed, env))
  idx <- caret::createFolds(y = response, k = 10, returnTrain = TRUE)
  # inject the index in the trainControl object
  caret::trainControl(method = "cv", index = idx,
    savePredictions = TRUE, selectionFunction = "oneSE")
}

## Adapt model tuning to leave-one-out cross-validation ========================

## trainControl generating helper function
control_train_loocv_q <- function(x, response, env = parent.frame()) {
  calibration <- NULL
  # r: response
  response <- eval(response, x$calibration, env)
  # Set up leave-one-out cross-validation
  caret::trainControl(method = "LOOCV", savePredictions = TRUE,
    selectionFunction = "oneSE")
}

## Adapt model tuning to repeated k-fold cross-validation ======================

## trainControl generating helper function
control_train_rcv_q <- function(x, response, resampling_seed,
                                env = parent.frame()) {
  calibration <- NULL
  # r: response
  response <- eval(response, x$calibration, env)
  # Set seed for creating resampling indices
  set.seed(eval(resampling_seed, env))
  # Set up 5 times repeated 10-fold cross-validation
  idx <- caret::createMultiFolds(y = response, k = 10, times = 5) # update ***
  # Inject the index in the trainControl object
  caret::trainControl(method = "repeatedcv", index = idx,
    savePredictions = TRUE, selectionFunction = "oneSE")
}

## Fitting models without parameter tuning =====================================
# 5.9; https://topepo.github.io/caret/model-training-and-tuning.html

## trainControl generating helper function
control_train_none_q <- function(x, response, resampling_seed,
                                 env = parent.frame()) {
  calibration <- NULL
  # r: response
  response <- eval(response, x$calibration, env)
  # Set seed for creating resampling indices
  set.seed(eval(resampling_seed, env))
  # Set trainControl argument to "none" so that caret::train will only fit
  # one model to the entire training set;
  # use a fixed number of PLS components instead
  idx <- caret::createFolds(y = response, k = 10, returnTrain = TRUE) # update ***
  # inject the index in the trainControl object
  caret::trainControl(method = "none", index = idx, savePredictions = TRUE,
    selectionFunction = "oneSE")
}

## Standard evlauation version of trainContol helper function
control_train <- function(x, response, env = parent.frame()) {
  control_train_q(x, substitute(response), env)
}

# Fit a PLS regression model using the caret package ---------------------------

## Train a PLS regression model
train_pls_q <- function(x,
  evaluation_method = "test_resampling",
  response, tr_control, env = parent.frame(),
  pls_ncomp_max = 20, ncomp_fixed = 5,
  center, scale, tuning_method = "resampling") {
  # Fit a partial least square regression (pls) model
  # center and scale MIR (you can try without)
  calibration <- MIR <- NULL
  r <- eval(response, x$calibration, env)
  # ? Is it really necessary to evaluate this in the parent frame?
  pls_ncomp_max <- eval(pls_ncomp_max, envir = parent.frame())
  # Evaluate fixed number of PLS regression components
  # from ncomp_fixed object in parent frame (fit_pls function)
  ncomp_fixed <- eval(ncomp_fixed, envir = parent.frame())

  # Test whether the spectral object has the class "tibble"
  if (!tibble::is_tibble(x$calibration)) {
    stop("spec_chem needs to be of class tibble")
  }

  spc_pre <- data.table::rbindlist(x$calibration$spc_pre)
  if (scale == TRUE && center == TRUE) {
    if (tuning_method == "resampling") {
      # Fit model with parameter tuning
      pls_model <- caret::train(x = spc_pre, y = r,
        method = "pls",
        tuneLength = pls_ncomp_max,
        trControl = tr_control,
        preProcess = c("center", "scale"))
    } else if (tuning_method == "none") {
      # Fit model without parameter tuning
      pls_model <- caret::train(x = spc_pre, y = r,
        method = "pls",
        trControl = tr_control,
        preProcess = c("center", "scale"),
        tuneGrid = data.frame(ncomp = ncomp_fixed))
    }
  } else {
    # No centering and scaling!
    pls_model <- caret::train(x = spc_pre, y = r,
      method = "pls",
      tuneLength = pls_ncomp_max,
      trControl = tr_control)
  }
}


## Standard evaluation version for training a PLS regression model
train_pls <- function(x, response, evaluation_method = "resampling",
  env = parent.frame()) {
  train_pls_q(x = x, evaluation_method = substitute(evaluation_method),
    response = substitute(response), env
  )
}

# Fit a random forest model using the caret package ----------------------------

## Train a random forest model
train_rf_q <- function(x,
  validation = TRUE, evaluation_method = "resampling",
  response, tr_control, ntree_max = 500, env = parent.frame()) {
  # Fit a partial least square regression (pls) model
  # center and scale MIR (you can try without)
  calibration <- MIR <- NULL
  response <- eval(response, x$calibration, env)
  ntree_max <- eval(ntree_max, envir = parent.frame())
  if (tibble::is_tibble(x$calibration)) {
    spc_pre <- data.table::rbindlist(x$calibration$spc_pre)
    rf_model <- caret::train( x = spc_pre, y = response,
      method = "rf",
      ntree = ntree_max,
      trControl = tr_control,
      preProcess = c("center", "scale")
    )

  } else {
    rf_model <- caret::train(x = x$calibration$MIR, y = response,
      method = "rf",
      ntree = ntree_max,
      trControl = tr_control,
      preProcess = c("center", "scale")
    )
  }
  rf_model
}


# Evaluate model performance (validation and cross-validation) -----------------

## Helper function to transform repeated k-fold cross-validation hold-out
## predictions
transform_cvpredictions <- function(cal_index, predobs_cv) {

  predobs_cv <- dplyr::full_join(cal_index, predobs_cv, by = "rowIndex") %>%
    dplyr::group_by(!!rlang::sym("sample_id")) %>%
    # Average observed and predicted values
    dplyr::mutate("obs" = mean(!!rlang::sym("obs")),
      "pred_sd" = sd(!!rlang::sym("pred"))) %>%
    # Add 95% confidence interval for mean hold-out predictions from
    # repeated k-fold cross-validation
    dplyr::mutate_at(.vars = dplyr::vars(!!rlang::sym("pred")),
      .funs = dplyr::funs("pred_sem_ci" = sem_ci)) %>%
    # Add mean hold-out predictions from repeated k-fold cross-validation
    dplyr::mutate("pred" = mean(!!rlang::sym("pred"))) %>%
    # Slice data set to only have one row per sample_id
    dplyr::slice(1L)
}

## Evaluate PLS performance
evaluate_model_q <- function(x, model, response,
  evaluation_method, tuning_method, resampling_method,
  print = TRUE, env = parent.frame()) {
  # Set global variables to NULL to avoid R CMD check notes
  MIR <- object <- dataType <- obs <- pred_sem_ci <- pred <- NULL
  ncomp <- finalModel <- rmse <- r2 <- r2 <- rpd <- n <- NULL
  rmse <- calibration <- NULL
  # Collect fitted object into a list
  list_models <- list("final_model" = model)
  # Evaluate validation argument in parent.frame !!!
  evaluation_method <- eval(evaluation_method, envir = parent.frame())
  # Evaluate tuning_method argument in parent.frame
  tuning_method <- eval(tuning_method, envir = parent.frame())
  # Evaluate resampling_method argument in parent.frame
  resampling_method <- eval(resampling_method, envir = parent.frame())
  # Extract best tuning parameters and associated cv predictions
  if (evaluation_method == "test_set") {
    # Calculate training (calibration) and test (validation) data
    # predictions based on pls model with calibration data
    r <- eval(response, x$validation, env)

    if (!tibble::is_tibble(x$validation)) {
      stop("Spectra and reference data need to be provided as tibble
        (class `tbl_df`, `tbl`, `data.frame`")
    }
    spc_pre <- data.table::rbindlist(x$validation$spc_pre)
    predobs <- caret::extractPrediction(list_models,
      testX = spc_pre, testY = r) # update ***
    # Append sample_id column to predobs data.frame
    # extract sample_id from validation set
    predobs$sample_id <- c(
      x$calibration$sample_id, x$validation$sample_id)

    # Create new data frame column <object>
    predobs$object <- predobs$model

    # Replace levels "Training" and "Test" in dataType column
    # by "Calibration" and "Validation" (rename levels of factor)
    predobs$dataType <- plyr::revalue(predobs$dataType,
      c("Test" = "Validation", "Training" = "Calibration")
    )
    # Change the order of rows in the data frame
    # Calibration as first level (show Calibration in ggplot graph
    # on left panel)
    predobs$dataType <- factor(predobs$dataType,
      levels = c("Calibration", "Validation"))
    # Calculate model performance indexes by model and dataType
    # uses package plyr and function summary.df of SPECmisc.R
    stats <- plyr::ddply(predobs, c("model", "dataType"),
      function(x) summary_df(x, "obs", "pred")
    )
  # Check whether method = "none" argument is selected in train();
  # this is the case when ncomp_fixed argument in fit_pls() is
  # evaluated
  # Checking for the existence of a <pred> element in the train function output
  # list can be dangerous and doesn't work in all cases when using
  # e.g. is.element('pred', x) or is.null(x$pred);
  # Problems can occur e.g. if a list element contains NULL element;
  # see
  # http://stackoverflow.com/questions/7719741/how-to-test-if-list-element-exists
  } else if (evaluation_method == "resampling" && tuning_method == "resampling") {
    # Good discussion on which cross-validation results are returned from caret
    # Extract best tuning parameters and associated cv predictions
    # http://stats.stackexchange.com/questions/219154/how-does-cross-validation-in-train-caret-precisely-work
    # Alternative solution for one model: conformal::GetCVPreds(model) function
    # see https://github.com/cran/conformal/blob/master/R/misc.R
    predobs_cv <- plyr::ldply(list_models,
      function(x) dplyr::anti_join(x$pred, x$bestTune, by = "ncomp"),
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
    vars_indexing <- c("rowIndex", "sample_id")
    cal_index <- dplyr::select(x$calibration,
      !!!rlang::syms(vars_indexing))

    # Transform cross-validation hold-out predictions --------------------------
    predobs_cv <- transform_cvpredictions(cal_index = cal_index,
      predobs_cv = predobs_cv)

    predobs_cv$object <- predobs_cv$model
    predobs_cv$model <- factor(predobs_cv$model)
    predobs_cv$dataType <- factor("Cross-validation")
    vars_keep <- c("obs", "pred", "pred_sd", "pred_sem_ci",
      "model", "dataType", "object")
    predobs_cv <- dplyr::select(predobs_cv,
      # !!! sample_id newly added
      !!!rlang::syms(vars_keep)
    )
    # Add column pred_sd to predobs data frame (assign values to 0) so that
    # column pred_sd is retained in predobs_cv after dplyr::bind_rows
    predobs$pred_sd <- NA
    # Desn't work because some columns are turned into numeric;
    # resulting data frame has only two rows
    # pb_2018-11-09: Model evaluation graph should only show cross-validation
    #   results when arguments `evaluation_method` == "resampling" &&
    #   `tuning_method` == "resampling"
    predobs <- predobs_cv
    # predobs <- dplyr::bind_rows(predobs, predobs_cv)
    # Calculate model performance indexes by model and dataType
    # uses package plyr and function summary.df of SPECmisc.R
    stats <- suppressWarnings(plyr::ddply(predobs, c("model", "dataType"),
      function(x) summary_df(x, "obs", "pred"))
    )
  }
  # Add number of components to stats; from finalModel list item
  # from train() function output (function from caret package)
  stats$ncomp <- rep(model$finalModel$ncomp, nrow(stats))
  # !!! Experimental: return stats
  # return(stats)
  # Add range of observed values for validation and calibraton
  # get range from predicted vs. observed data frame
  # stored in object predobs
  obs_cal <- subset(predobs, dataType == "Calibration")$obs
  # Get name of predicted variable; see p. 261 of book
  # "Advanced R" (Hadley Wickham)
  response_name <- deparse(response)

  if (evaluation_method == "test_set") {
    # Assign validation set to separate data frame
    obs_val <- subset(predobs, dataType == "Validation")$obs
    # before: deparse(substitute(variable))
    df_range <- data.frame(
      response = rep(response_name, 2),
      dataType = c("Calibration", "Validation"),
      min_obs = c(range(obs_cal)[1], range(obs_val)[1]),
      median_obs = c(median(obs_cal), median(obs_val)),
      max_obs = c(range(obs_cal)[2], range(obs_val)[2]),
      mean_obs = c(mean(obs_cal), mean(obs_val)),
      CV = c(sd(obs_cal) / mean(obs_cal) * 100,
        sd(obs_val) / mean(obs_val) * 100)
    )
  } else if (evaluation_method == "resampling" && tuning_method == "resampling") {
    # Assign cross-validation set to separate data frame
    obs_val <- subset(predobs, dataType == "Cross-validation")$obs
    df_range <- data.frame(
      response = rep(response_name, 2),
      dataType = factor(c("Calibration", "Cross-validation")),
      min_obs = c(range(obs_cal)[1], range(obs_val)[1]),
      median_obs = c(median(obs_cal), median(obs_val)),
      max_obs = c(range(obs_cal)[2], range(obs_val)[2]),
      mean_obs = c(mean(obs_cal), mean(obs_val)),
      CV = c(sd(obs_cal) / mean(obs_cal) * 100,
        sd(obs_val) / mean(obs_val) * 100)
    )
  }

  # Join stats with range data frame (df_range)
  stats <- suppressWarnings(dplyr::inner_join(stats, df_range, by = "dataType"))
  annotation <- plyr::mutate(stats,
    rmse = as.character(as.expression(paste0("RMSE == ",
      round(rmse, 2)))),
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
  make_label <- function(x, evaluation_method = "test_set",
                         resampling_method = "kfold_cv") {
    dataType <- n <- NULL

    if (evaluation_method == "test_set") {
      c(`Calibration` = paste0("Calibration", "~(",
        x[x$dataType == "Calibration", ]$n, ")"
      ),
        `Validation` = paste0("Validation", "~(",
          x[x$dataType == "Validation", ]$n, ")"
        )
      )
    } else if (evaluation_method == "resampling" &&
        resampling_method == "rep_kfold_cv") {
      c(`Calibration` = paste0("Calibration", "~(",
        x[x$dataType == "Calibration", ]$n, ")"
      ),
        `Cross-validation` = paste0("5%*%repeated~10*-fold~CV", "~(",
          x[x$dataType == "Cross-validation", ]$n, ")"
        )
      )
    } else {
      c(`Calibration` = paste0("Calibration", "~(",
        x[x$dataType == "Calibration", ]$n, ")"
      ),
        `Cross-validation` = paste0("10*-fold~CV", "~(",
          x[x$dataType == "Cross-validation", ]$n, ")"
        )
      )
    }
  }
  if (evaluation_method == "test_set") {
    label_validation <- make_label(x = annotation,
      evaluation_method = "test_set"
    )
  } else if (evaluation_method == "resampling" &&
    resampling_method == "rep_kfold_cv") {
    label_validation <- make_label(x = annotation,
      evaluation_method = "resampling", resampling_method = "rep_kfold_cv"
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

  # Create model evaluation plot -----------------------------------------------

  ## ggplot graph for model comparison
  ## (arranged later in panels)
  x_label <- paste0("Observed ",
    as.character(response_name))
  y_label <- paste0("Predicted ",
    as.character(response_name))

  ## Create x and y minimum and maximum for plotting range; use either
  ## observed or predicted data, depending on what minimum and maximum values
  ## are
  xy_min <- if (min(predobs$obs) < min(predobs$pred))
    {predobs$obs} else {predobs$pred}
  xy_max <- if (max(predobs$obs) > max(predobs$pred))
    {predobs$obs} else {predobs$pred}
  xy_range <- ifelse(diff(range(xy_min) > diff(range(xy_max))),
    diff(range(xy_min)), diff(range(xy_max)))

  if (model$method == "pls") {
  p_model <- ggplot2::ggplot(data = predobs) +
    ggplot2::geom_point(ggplot2::aes(x = obs, y = pred),
      shape = 1, size = 2, alpha = 1/2, data = predobs) +
    ggplot2::geom_text(data = annotation,
      ggplot2::aes(x = Inf, y = -Inf, label = ncomp), size = 5,
      hjust = 1.15, vjust = -4.5, parse = TRUE) + # !!! additional label
    ggplot2::geom_text(data = annotation,
      ggplot2::aes(x = Inf, y = -Inf, label = r2), size = 5,
      hjust = 1.15, vjust = -3, parse = TRUE) +
    ggplot2::geom_text(data = annotation,
      ggplot2::aes(x = Inf, y = -Inf, label = rmse), size = 5,
      hjust = 1.12, vjust = -2.5, parse = TRUE) +
    ggplot2::geom_text(data = annotation,
      ggplot2::aes(x = Inf, y = -Inf, label = rpd), size = 5,
      hjust = 1.15, vjust = -1.25, parse = TRUE) +
    ggplot2::facet_grid(~ dataType,
      labeller = ggplot2::as_labeller(to_string)) +
    ggplot2::geom_abline(col = "red") +
    ggplot2::labs(x = x_label, y = y_label) +
    ggplot2::xlim(c(min(xy_min) - 0.05 * xy_range,
      max(xy_max) + 0.05 * xy_range)) +
    ggplot2::ylim(c(min(xy_min) - 0.05 * xy_range,
      max(xy_max) + 0.05 * xy_range)) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"))

    if (evaluation_method == "resampling") {
    p_model <- p_model +
      ggplot2::geom_errorbar(
        ggplot2::aes(x = obs, ymin = pred - pred_sem_ci,
        ymax = pred + pred_sem_ci),
        width = 0, data = predobs, inherit.aes = FALSE)
    }

  } else {

  p_model <- ggplot2::ggplot(data = predobs) +
    ggplot2::geom_point(ggplot2::aes(x = obs, y = pred),
      shape = 1, size = 2, alpha = 1/2) + # without ncomp label
    ggplot2::geom_text(data = annotation,
      ggplot2::aes(x = Inf, y = -Inf, label = r2), size = 5,
      hjust = 1.15, vjust = -3, parse = TRUE) +
    ggplot2::geom_text(data = annotation,
      ggplot2::aes(x = Inf, y = -Inf, label = rmse), size = 5,
      hjust = 1.12, vjust = -2.5, parse = TRUE) +
    ggplot2::geom_text(data = annotation,
      ggplot2::aes(x = Inf, y = -Inf, label = rpd), size = 5,
      hjust = 1.15, vjust = -1.25, parse = TRUE) +
    ggplot2::facet_grid(~ dataType,
      labeller = ggplot2::as_labeller(to_string)) +
    ggplot2::geom_abline(col = "red") +
    ggplot2::labs(x = x_label, y = y_label) +
    ggplot2::xlim(c(min(xy_min) - 0.05 * xy_range,
      max(xy_max) + 0.05 * xy_range)) +
    ggplot2::ylim(c(min(xy_min) -
      0.05 * xy_range, max(xy_max) + 0.05 * xy_range)) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw()
  }

  if (print == TRUE) {
    print(p_model)
  }

  list(stats = stats, p_model = p_model, predobs = predobs)
}


## PLS regression modeling in one function ======================

#' @name fit_pls
#' @title Calibration sampling, model tuning, and PLS regression
#' @description Perform calibration sampling and use selected
#' calibration set for model tuning
#' @param spec_chem Tibble that contains spectra, metadata and chemical
#' reference as list-columns. The tibble to be supplied to \code{spec_chem} can
#' be generated by the `join_chem_spc() function`
#' @param response Response variable as symbol or name
#' (without quotes, no character string). The provided response symbol needs to be
#' a column name in the \code{spec_chem} tibble.
#' @param variable Depreciated and replaced by `response`
#' @param center Logical whether to perform mean centering of each spectrum column
#' (e.g. wavenumber or wavelength) after common spectrum preprocessing. Default is
#' \code{center = TRUE}
#' @param scale Logical whether to perform standard deviation scaling
#' of each spectrum column (e.g. wavenumber or wavelength) after common
#' spectrum preprocessing. Default is \code{scale = TRUE}
#' @param evaluation_method Character string stating evaluation method.
#' Either \code{"test_set"} (default) or \code{"resampling"}. \code{"test_set"}
#' will split the data into a calibration (training) and validation (test) set,
#' and evaluate the final model by predicting on the validation set.
#' If \code{"resampling"}, the finally selected model will be evaluated based
#' on the cross-validation hold-out predictions.
#' @param validation Depreciated and replaced by \code{evaluation_method}.
#' Default is \code{TRUE}.
#' @param split_method Method how to to split the data into a independent test
#' set. Default is \code{"ken_sto"}, which will select samples for calibration
#' based on Kennard-Stone sampling algorithm of preprocessed spectra. The
#' proportion of validation to the total number of samples can be specified
#' in the argument \code{ratio_val}.
#' \code{split_method = "random"} will create a single random split.
#' @param ratio_val Ratio of validation (test) samples to
#' total number of samples (calibration (training) and validation (test)).
#' @param ken_sto_pc Number of component used
#' for calculating mahalanobsis distance on PCA scores for computing
#' Kennard-Stone algorithm.
#' Default is \code{ken_sto_pc = 2}, which will use the first two PCA
#' components.
#' @param pc Depreciated; renamed argument is `ken_sto_pc`.
#' @param invert Logical
#' @param tuning_method Character specifying tuning method. Tuning method
#' affects how caret selects a final tuning value set from a list of candidate
#' values. Possible values are \code{"resampling"}, which will use a
#' specified resampling method such as repeated k-fold cross-validation (see
#' argument \code{resampling_method}) and the generated performance profile
#' based on the hold-out predictions to decide on the final tuning values
#' that lead to optimal model performance. The value \code{"none"} will force
#' caret to compute a final model for a predefined canditate PLS tuning
#' parameter number of PLS components. In this case, the value
#' supplied by \code{ncomp_fixed}` is used to set model complexity at
#' a fixed number of components.
#' @param resampling_method Character specifying resampling method. Currently,
#' \code{"kfold_cv"} (default, performs 10-fold cross-validation),
#' \code{"rep_kfold_cv"} (performs 5-times repeated 10-fold cross-validation),
#' \code{"loocv"} (performs leave-one-out cross-validation), and \code{"none"}
#' (if \code{resampling_method = "none"}) are supported.
#' @param resampling_seed Random seed (integer) that will be used for generating
#' resampling indices, which will be supplied to \code{caret::trainControl}.
#' This makes sure that modeling results are constant when re-fitting.
#' Default is \code{resampling_seed = 123}.
#' @param cv Depreciated. Use \code{resampling_method} instead.
#' @param pls_ncomp_max Maximum number of PLS components that are evaluated
#' by caret::train. Caret will aggregate a performance profile using resampling
#' for an integer sequence from 1 to \code{pls_ncomp_max}
#' @param ncomp_fixed Integer of fixed number of PLS components. Will only be
#' used when \code{tuning_method = "none"} and  \code{resampling_method = "none"}
#' are used.
#' @param print Logical expression whether model evaluation graphs shall be
#' printed
#' @param env Environment where function is evaluated. Default is
#' \code{parent.frame}.
#' @export
# Note: check non standard evaluation, argument passing...
fit_pls <- function(
  spec_chem,
  response, variable = NULL, # variable depreciated, will not work
  center = TRUE, scale = TRUE, # center and scale all predictors (wavenumbers)
  evaluation_method = "test_set", validation = TRUE, # validation depreciated
  split_method = "ken_stone",
  ratio_val = 1/3, # is only used if evaluation_method = "test_set"
  ken_sto_pc = 2, pc, # only if split_method = "ken_stone"; number of component
  # used for calculating mahalanobsis distance on PCA scores. pc is depreciated.
  invert = TRUE, # only if split_method = "ken_stone"
  tuning_method = "resampling",
  resampling_method = "kfold_cv", cv = NULL, # cv depreciated
  resampling_seed = 123, # Seed for creating resampling indices
  pls_ncomp_max = 20, # Maximal number of PLS components used by model tuning
  ncomp_fixed = 5, # only fit and evaluate one model, if tuning_method = "none"
  print = TRUE, # print model summary and evaluation graphs
  env = parent.frame())

{
  calibration <- 0

  # Warning messages and reassignment for depreciated arguments ----------------
  # Depreciate argument variable, use more specific term for the response
  # to be predicted by spectral modeling
  if (!is.null(variable)) {
    stop("argument variable has been replaced by response for simplerspec_0.1.0")
  }
  # 20170602: revise argument name and values of validation;
  # Replace validation = TRUE or FALSE with
  # new argument evaluation_method = "test_set" or "resampling"
  if (!missing(validation)) {
    warning("argument validation is depreciated; please use evaluation_method instead.",
      call. = FALSE)
    evaluation_method <- validation
  }
  # Depreciate argument pc, use more consistent and verbose argument ken_sto_pc
  if (!missing(pc)) {
    warning("argument pc is depreciated; please use ken_sto_pc instead.",
      call. = FALSE)
    ken_sto_pc <- pc
  }
  # Depreciate argument cv, use more consistent and verbose argument
  # resampling_method
  if (!missing(cv)) {
    warning("argument cv is depreciated; please use resampling_method instead.",
      call. = FALSE)
    resampling_method <- cv
  }
  # Change values for resampling_method argument
  if (resampling_method == "LOOCV") {
    warning("value 'LOOCV' (leave one out cross-validation) for argument resampling_method is depreciated; please use value 'loocv' instead.")
    resampling_method <- "loocv"
  }
  if (resampling_method == "repeatedcv") {
     warning("value 'repeatedcv' (repeated k-fold cross-validation) for argument resampling_method is depreciated; please use value 'rep_kfold_cv' instead.")
    resampling_method <- "rep_kfold_cv"
  }

  # Perform calibration sampling -----------------------------------------------
  list_sampled <- split_data_q(
    spec_chem, split_method, ratio_val = ratio_val,
    ken_sto_pc = substitute(ken_sto_pc),
    evaluation_method = substitute(evaluation_method),
    invert = substitute(invert)
  )

  # Check on method for cross-validation to be used in caret model tuning ------
  if (resampling_method == "loocv") {
    # leave-one-out cross-validation
    tr_control <- control_train_loocv_q(x = list_sampled,
      response = substitute(response), env = env)
  } else if (resampling_method == "rep_kfold_cv") {
    # repeated k-fold cross-validation
    tr_control <- control_train_rcv_q(x = list_sampled,
      response = substitute(response),
      resampling_seed = substitute(resampling_seed), env = env)
  } else if (resampling_method == "none") {
    # no resampling; calls caret::train(..., method = "none");
    # fixed number of PLS components; tuning_method argument has also
    # to be set to "none"
    tr_control <- control_train_none_q(x = list_sampled,
      response = substitute(response),
      resampling_seed = substitute(resampling_seed), env = env)
  } else if (resampling_method == "kfold_cv") {
    # k-fold cross validation
    tr_control <- control_train_q(x = list_sampled,
      response = substitute(response),
      resampling_seed = substitute(resampling_seed), env = env)
  }
  # Fit a pls calibration model; pls object is output from caret::train() ------
  if (tuning_method == "resampling") {
    pls <- train_pls_q(x = list_sampled,
      evaluation_method = "test_set",
      response = substitute(response), tr_control = tr_control,
      center = center, scale = scale,
      pls_ncomp_max = substitute(pls_ncomp_max), env
      )
  } else if (tuning_method == "none") {
    pls <- train_pls_q(x = list_sampled,
      evaluation_method = "test_set",
      response = substitute(response), tr_control = tr_control,
      center = center, scale = scale, tuning_method = "none",
      ncomp_fixed = substitute(ncomp_fixed), env
      )
  }
  # Evaluate model accuracy (predicted vs. observed) ---------------------------
  stats <- evaluate_model_q(x = list_sampled, model = pls,
    response = substitute(response),
    evaluation_method = substitute(evaluation_method),
    tuning_method = substitute(tuning_method),
    resampling_method = substitute(resampling_method),
    env = parent.frame()
  )
  list(data = list_sampled, p_pc = list_sampled$p_pc,
    model = pls, stats = stats$stats, p_model = stats$p_model,
    predobs = stats$predobs)
}

## Old function name of `fit_pls`: `pls_ken_stone`

#' @rdname fit_pls
#' @export
pls_ken_stone <- fit_pls

## Random forest modeling in one function =======================

#' @title Calibration sampling, and random forest model tuning and evaluation
#' @description Perform calibration sampling and use selected
#' calibration set for model tuning
#' @param spec_chem Tibble that contains spectra, metadata and chemical
#' reference as list-columns. The tibble to be supplied to \code{spec_chem} can
#' be generated by the `join_chem_spc() function`
#' @param response Response variable as symbol or name
#' (without quotes, no character string). The provided response symbol needs to be
#' a column name in the \code{spec_chem} tibble.
#' @param variable Depreciated and replaced by `response`
#' @param evaluation_method Character string stating evaluation method.
#' Either \code{"test_set"} (default) or \code{"resampling"}. \code{"test_set"}
#' will split the data into a calibration (training) and validation (test) set,
#' and evaluate the final model by predicting on the validation set.
#' If \code{"resampling"}, the finally selected model will be evaluated based
#' on the cross-validation hold-out predictions.
#' @param validation Depreciated and replaced by \code{evaluation_method}.
#' Default is \code{TRUE}.
#' @param split_method Method how to to split the data into a independent test
#' set. Default is \code{"ken_sto"}, which will select samples for calibration
#' based on Kennard-Stone sampling algorithm of preprocessed spectra. The
#' proportion of validation to the total number of samples can be specified
#' in the argument \code{ratio_val}.
#' \code{split_method = "random"} will create a single random split.
#' @param ratio_val Ratio of validation (test) samples to
#' total number of samples (calibration (training) and validation (test)).
#' @param ken_sto_pc Number of component used
#' for calculating mahalanobsis distance on PCA scores for computing
#' Kennard-Stone algorithm.
#' Default is \code{ken_sto_pc = 2}, which will use the first two PCA
#' components.
#' @param pc Depreciated; renamed argument is `ken_sto_pc`.
#' @param invert Logical
#' @param tuning_method Character specifying tuning method. Tuning method
#' affects how caret selects a final tuning value set from a list of candidate
#' values. Possible values are \code{"resampling"}, which will use a
#' specified resampling method such as repeated k-fold cross-validation (see
#' argument \code{resampling_method}) and the generated performance profile
#' based on the hold-out predictions to decide on the final tuning values
#' that lead to optimal model performance. The value \code{"none"} will force
#' caret to compute a final model for a predefined canditate PLS tuning
#' parameter number of PLS components. In this case, the value
#' supplied by \code{ncomp_fixed}` is used to set model complexity at
#' a fixed number of components.
#' @param resampling_seed Random seed (integer) that will be used for generating
#' resampling indices, which will be supplied to \code{caret::trainControl}.
#' This makes sure that modeling results are constant when re-fitting.
#' Default is \code{resampling_seed = 123}.
#' @param cv Depreciated. Use \code{resampling_method} instead.
#' @param ntree_max Maximum random forest trees
#' by caret::train. Caret will aggregate a performance profile using resampling
#' for an integer sequence from 1 to \code{ntree_max} trees.
#' @param print Logical expression whether model evaluation graphs shall be
#' printed
#' @param env Environment where function is evaluated. Default is
#' \code{parent.frame}.
#' @export
# Note: check non standard evaluation, argument passing...
fit_rf <- function(spec_chem,
    response, variable = NULL, # variable depreciated, will not work
    evaluation_method = "test_set", validation = NULL, # Validation is depreciated
    split_method = "ken_stone",
    ratio_val,
    ken_sto_pc = 2, pc = NULL,
    invert = TRUE, # only if split_method = "ken_stone
    tuning_method = "resampling",
    resampling_seed = 123,
    cv = NULL, # cv depreciated
    ntree_max = 500,
    print = TRUE,
    env = parent.frame()) {

  calibration <- NULL

  # Warning messages and reassignment for depreciated arguments ----------------
  # Replace validation = TRUE or FALSE with
  # new argument evaluation_method = "test_set" or "resampling"
  if (!missing(validation)) {
    warning("argument validation is deprecated; please use evaluation_method instead.",
      call. = FALSE)
    evaluation_method <- validation
  }
  if (!missing(variable)) {
    warning("argument variable is deprecated; please use response instead.",
      call. = FALSE)
    response <- variable
  }
  # Depreciate argument pc, use more consistent and verbose argument ken_sto_pc
  if (!missing(pc)) {
    warning("argument pc is depreciated; please use ken_sto_pc instead.",
      call. = FALSE)
    ken_sto_pc <- pc
  }

  # Calibration sampling -------------------------------------------------------
  list_sampled <- split_data_q(
    spec_chem, split_method, ratio_val = ratio_val,
    ken_sto_pc = substitute(ken_sto_pc),
    evaluation_method = substitute(evaluation_method),
    invert = substitute(invert)
  )
  # Control parameters for caret::train ----------------------------------------
  tr_control <- control_train_q(x = list_sampled,
    response = substitute(response),
    resampling_seed = substitute(resampling_seed), env = env)
  # Train random forest model (model tuning) -----------------------------------
  rf <- train_rf_q(x = list_sampled,
    evaluation_method = "test_set",
    response = substitute(response), tr_control = tr_control, env,
    ntree_max = substitute(ntree_max)
  )
  # Evaluate finally chosen random forest model --------------------------------
  stats <- evaluate_model_q(x = list_sampled, model = rf,
    response = substitute(response),
    evaluation_method = substitute(evaluation_method),
    tuning_method = substitute(tuning_method),
    resampling_method = substitute(resampling_method),
    env = parent.frame()
  )
  # Return list with results ---------------------------------------------------
  list(data = list_sampled, p_pc = list_sampled$p_pc,
    rf_model = rf, stats = stats$stats, p_model = stats$p_model)
}
