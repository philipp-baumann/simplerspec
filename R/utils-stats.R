#' @title Calculate model evaluation metrics
#' @description Calculates model statistics for predicted (y)
#' vs. observed (y) values
#' @param df data.frame with predicted and observed data
#' @param x column with observed values
#' @param y column with predicted values
#' @export
summary_df <- function(df, x, y) {
  # !!! note that y are predicted values and x are observed values
  x <- rlang::enquo(x)
  y <- rlang::enquo(y)
  x <- dplyr::pull(df, !!x)
  y <- dplyr::pull(df, !!y)

  tibble::tibble(
    n = length(x),
    min = min(x, ra.rm = TRUE),
    max = max(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sdev = sd(x, na.rm = TRUE),
    rmse = mean((y - x)^2, na.rm = TRUE)^.5,
    mse = mean((y - x)^2, na.rm = TRUE),
    r2  = cor(x, y, use = "pairwise.complete.obs")^2,
    b = lm(x ~ y)$coefficients[2],
    rpd = sd(x, na.rm = TRUE) /
      sqrt(sum((y - x)^2, na.rm = TRUE) / (length(x) - 1)),
    rpiq = (quantile(x, .75, na.rm = TRUE) - quantile(x, .25, na.rm = TRUE)) /
      sqrt(sum((x - y)^2, na.rm = TRUE) / (length(x) - 1)),
    bias  = mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE),
    SB = (mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE))^2,
    NU = mean((y - mean(y))^2) * (1 - lm(x ~ y)$coefficients[2])^2,
    LC = mean((x - mean(x))^2)
      * (1 - cor(x, y, use = "pairwise.complete.obs")^2),
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
