#' @title Assess multiple pairs of measured and predicted values
#' @description Return performance metrics for test set predictions and
#' measured values, e.g. for different model outcome variables.
#' @param data Data frame with all measured (observed) and predicted variables.
#' @param ... Multiple arguments with observed (measured)-predicted pairs,
#' specified with \code{dplyr::vars(o = <column_name>, p = <column_name>)}.
#' Column names can strings or symbols. The arguments in `...` need to be named.
#' @param .metrics Character vector with package used for metrics calculation.
#' Default is \code{"simplerspec"}, which uses
#' \code{simplerspec::evaluate_model()}.
#' @param .model_name String with name for the new column that specifies the
#' model or the outcome variable. Default is \code{"model"}.
#'
#' @return Data frame with with summary statistics for measured values and
#' performance metrics for the pairs of measured and predicted values.
#' @importFrom purrr modify_depth imap
#' @export
assess_multimodels <- function(data,
                               ...,
                               .metrics = c("simplerspec", "yardstick"),
                               .model_name = "model") {
  args <- rlang::enquos(...)
  args_tidy <- map(args, rlang::eval_tidy)
  stopifnot(
    all(map_int(names(args), nchar)),
    all(map_lgl(args_tidy, is.list))
  )

  vars_names <- modify_depth(args_tidy, .depth = 1, names)
  vars_names_ok <- map_lgl(vars_names, ~ all(c("o", "p") %in% .x))
  if (!all(vars_names_ok)) {
    stop("Assessment variables supplied in `vars()` need to be named with
      'o' (observed) and 'p' (predicted)")}

  metrics <- match.arg(.metrics)
  # pb 2018-05-09: todo: support yardstick metrics
  assessment <- switch (metrics,
    "simplerspec" = map(
      .x = args_tidy,
      ~ evaluate_model(data = data,
        obs = !!.x[["o"]], pred = !!.x[["p"]]))
  )
  assessment_models <- imap(assessment,
    ~ tibble::add_column(.x, !!.model_name := .y, .before = 1))

  dplyr::bind_rows(assessment_models)
}


#' @title Calculate model evaluation metrics
#' @description Calculate observed summary statistics and model evaluation
#' statistics for assessing agreement between observed (`obs`) and
#' predicted (`pred`) values.
#' @param data `data.frame`` with predicted and observed data in columns.
#' @param obs Column that contains observed values, `symbol`/`name` or
#' `character` (wrapped in "").
#' @param pred Column that contains predicted values, `symbol`/`name` or
#' `character` (wrapped in "").
#' @importFrom e1071 kurtosis
#' @export

# Note that coefficient of determination (r2) derived from a linear regression
# of observed values on the prediction
# solely describes what proportion of variance in the measured data is
# simulated by the model. A linear regression line between observed (x) and
# predicted (y) does not provide a measure of model error!

# Mean squared error (MSE) can be decomposed into squared bias/deviation (SE^2)
# and mean squared variation (MSV) (:= SDE^2 := "squared standard deviation of
# the error"; see e.g. Kobayashi and Salam (2000)
# Gauch et al. (2003) propose a more sophisticated additive partitioning of the
# MSE that is more informative about the sources of error and link to the
# regression parameters; namely, these are squared bias (SB),
# non-unity slope (NU) and lack of correlation (LC)

evaluate_model <- function(data, obs, pred) {
  # Implement quasiquotation to first quote the `obs` and `pred` arguments, and
  # unquote/evaluate in the context of the `data` data.frame;
  # obs` and `pred` can be # both symbols or character;
  # `obs = obs` or `obs = "obs"` when the column that contains observed values
  # is named "obs"
  obs <- rlang::enquo(obs)
  pred <- rlang::enquo(pred)
  obs <- dplyr::pull(data, !!obs)
  pred <- dplyr::pull(data, !!pred)

  tibble::tibble(
    ## Compute descriptive statistics of the observations/measurements
    n = length(obs),
    min = min(obs, ra.rm = TRUE),
    max = max(obs, na.rm = TRUE),
    mean = mean(obs, na.rm = TRUE),
    median = median(obs, na.rm = TRUE),
    sdev = sd(obs, na.rm = TRUE),
    cv = sd(obs, na.rm = TRUE) / mean(obs, na.rm = TRUE),
    skewness_b1 = e1071::skewness(obs, na.rm = TRUE, type = 3),
    kurtosis = e1071::kurtosis(obs, na.rm = TRUE),

    ## Compute model evaluation measures to address different aspects
    ## of how well predictions correspond to observations/measurements
    # Root mean squared error
    rmse = mean((obs - pred)^2, na.rm = TRUE)^.5,
    # Mean squared error; mse^2 = me^2 + msv = me^2 + sde^2
    mse = mean((obs - pred)^2, na.rm = TRUE),
    # Terms mean error (ME) and bias are equivalent
    me = mean(obs - pred, na.rm = TRUE),
    bias = mean(obs - pred, na.rm = TRUE),
    # Mean squared variation (of the error); difference between the simulation
    # and the measurement with respect to the deviation from the means
    msv = mean(((mean(obs, na.rm = TRUE) - obs)
      - (mean(pred, na.rm = TRUE) - pred))^2),
    # Standard deviation of the error := SDE = MSV^0.5
    sde = mean(((mean(obs, na.rm = TRUE) - obs)
      - (mean(pred, na.rm = TRUE) - pred))^2)^0.5,
    # Mean absolute error
    mae = mean(abs(obs - pred), na.rm = TRUE),
    r2  = cor(obs, pred, use = "pairwise.complete.obs")^2,
    b = lm(obs ~ pred)$coefficients[2],
    # Ratio of performance to deviation
    rpd = sd(obs, na.rm = TRUE) / mean((obs - pred)^2, na.rm = TRUE)^.5,
    # Ratio of performance to interquartile range
    rpiq = (quantile(obs, .75, na.rm = TRUE) - quantile(obs, .25, na.rm = TRUE))
     / mean((obs - pred)^2, na.rm = TRUE)^.5,
    # See Gauch et. al., 2003) for MSD decomposition into SB, NU and LC;
    # Squared bias
    SB = (mean(obs - pred, na.rm = TRUE))^2,
    # Non-unity slope
    NU = mean((pred - mean(pred))^2) * (1 - lm(obs ~ pred)$coefficients[2])^2,
    # Lack of correlation
    LC = mean((obs - mean(obs))^2)
      * (1 - cor(obs, pred, use = "pairwise.complete.obs")^2),
    # Proportional contributions of SB, NU and LC to MSE in percent
    SB_prop = round((mean(obs - pred, na.rm = TRUE))^2
      / mean((pred - obs)^2) * 100, 0),
    NU_prop = round(mean((pred - mean(pred))^2)
      * (1 - lm(obs ~ pred)$coefficients[2])^2 / mean((pred - obs)^2) * 100, 0),
    LC_prop = round(mean((obs - mean(obs))^2)
      * (1 - cor(obs, pred, use = "pairwise.complete.obs")^2)
      / mean((pred - obs)^2) * 100, 0)
  )
}

# Wrapper function to ensure compatibility with old summary_df() function
# !!! note that y are predicted values and x are observed values;
# this is only for backward compatibility;
# This deviates from the principles that observed (O) should be denoted in the
# in the y-axis vs. predicted (P) in the x-axis ((OP) regressions), according
# Pineiro et al. (2008)
#' @rdname evaluate_model
#' @export
summary_df <- function(df, x, y) {
  evaluate_model(data = df, obs = x, pred = y)
}

# Function to calculate standard error of the mean
sem_ci <- function(x) {
  qt(0.975, df = length(na.omit(x)) - 1) *
    sqrt(var(x, na.rm = TRUE) / length(na.omit(x)))
}

# Calculate standard error
se <- function(x) {
  sqrt(var(x, na.rm = TRUE) / length(na.omit(x)))
}
