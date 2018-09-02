#' @title Calculate model evaluation metrics
#' @description Calculates model statistics for predicted (y)
#' vs. observed (y) values
#' @param df data.frame with predicted and observed data
#' @param x column with observed values
#' @param y column with predicted values
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

summary_df <- function(df, x, y) {
  # !!! note that y are predicted values and x are observed values
  x <- rlang::enquo(x)
  y <- rlang::enquo(y)
  x <- dplyr::pull(df, !!x)
  y <- dplyr::pull(df, !!y)

  tibble::tibble(
    ## Compute descriptive statistics of the observations/measurements
    n = length(x),
    min = min(x, ra.rm = TRUE),
    max = max(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sdev = sd(x, na.rm = TRUE),
    cv = sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE),
    skewness_b1 = e1071::skewness(x, na.rm = TRUE, type = 3),
    kurtosis = e1071::kurtosis(x, na.rm = TRUE),

    ## Compute model evaluation measures to address different aspects
    ## of how well predictions correspond to observations/measurements
    # Root mean squared error
    rmse = mean((y - x)^2, na.rm = TRUE)^.5,
    # Mean squared error; mse^2 = me^2 + msv = me^2 + sde^2
    mse = mean((y - x)^2, na.rm = TRUE),
    # Terms mean error (ME) and bias are equivalent
    me = mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE),
    bias = mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE),
    # Mean squared variation (of the error); difference between the simulation
    # and the measurement with respect to the deviation from the means
    msv <- mean(((mean(y, na.rm = TRUE) - y) - (mean(x, na.rm = TRUE) - x))^2),
    # Standard deviation of the error := SDE = MSV^0.5
    sde <- mean(((mean(y, na.rm = TRUE) - y)
      - (mean(x, na.rm = TRUE) - x))^2)^0.5,
    # Mean absolute error
    mae = mean(abs(y - x), na.rm = TRUE),
    r2  = cor(x, y, use = "pairwise.complete.obs")^2,
    b = lm(x ~ y)$coefficients[2],
    # Ratio of performance to deviation
    rpd = sd(x, na.rm = TRUE) /
      sqrt(sum((y - x)^2, na.rm = TRUE) / (length(x) - 1)),
    # Ratio of performance to interquartile range
    rpiq = (quantile(x, .75, na.rm = TRUE) - quantile(x, .25, na.rm = TRUE)) /
      sqrt(sum((x - y)^2, na.rm = TRUE) / (length(x) - 1)),
    # See Gauch et. al., 2003) for MSD decomposition into SB, NU and LC;
    # Squared bias
    SB = (mean(x - y, na.rm = TRUE))^2,
    # Non-unity slope
    NU = mean((y - mean(y))^2) * (1 - lm(x ~ y)$coefficients[2])^2,
    # Lack of correlation
    LC = mean((x - mean(x))^2)
      * (1 - cor(x, y, use = "pairwise.complete.obs")^2),
    # Proportional contributions of SB, NU and LC to MSE in percent
    SB_prop = round((mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE))^2
      / mean((y - x)^2) * 100, 0),
    NU_prop = round(mean((y - mean(y))^2) * (1 - lm(x ~ y)$coefficients[2])^2
      / mean((y - x)^2) * 100, 0),
    LC_prop = round(mean((x - mean(x))^2)
      * (1 - cor(x, y, use = "pairwise.complete.obs")^2)
      / mean((y - x)^2) * 100, 0)
  )
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
