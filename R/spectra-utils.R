#' @title Calculate model statistics
#' @description Calculates model statistics for predicted (y)
#' vs. observed (y) values
#' @param df data.frame with predicted and observed data
#' @param x column with observed values
#' @param y column with predicted values
#' @export
summary_df <- function(df, x, y){
  # !!! note that y are predicted values and x are observed values
  x <- df[, x]
  y <- df[, y]
  b <- lm(x ~ y)$coefficients[2]
  data.frame(
    rmse = sqrt(sum((x - y)^2, na.rm = T) / (length(x)-1)),
    rmsd = mean((y - x)^2)^.5,
    msd = mean((y - x)^2),
    sdev = sd(x, na.rm = T),
    rpd =  sd(x,na.rm = T) /
      sqrt(sum((y - x)^2, na.rm = T) / (length(x) - 1)),
    rpiq = (quantile(x, .75, na.rm = T) - quantile(x, .25, na.rm = T)) /
      sqrt(sum((x - y)^2, na.rm = T) / (length(x) - 1)),
    r2  = cor(x, y, use = "pairwise.complete.obs")^2,
    bias  = mean(y, na.rm = T) - mean(x, na.rm = T),
    SB = (mean(y, na.rm = T) - mean(x, na.rm = T))^2,
    NU = mean((y - mean(y))^2) * (1 - lm(x ~ y)$coefficients[2])^2,
    LC = mean((x - mean(x))^2) * (1 - cor(x, y, use = "pairwise.complete.obs")^2),
    SB_prop = round((mean(y, na.rm = T) - mean(x, na.rm = T))^2 / mean((y - x)^2) * 100, 0),
    NU_prop = round(mean((y - mean(y))^2) * (1 - lm(x ~ y)$coefficients[2])^2 / mean((y - x)^2) * 100, 0),
    LC_prop = round(mean((x - mean(x))^2) * (1 - cor(x, y, use = "pairwise.complete.obs")^2) / mean((y - x)^2) * 100, 0),
    n = length(x),
    b = b
  )
}

#' @export
# Function to calculate standard error of the mean
sem_ci <- function(x) {
  qt(0.975, df = length(na.omit(x)) - 1) *
    sqrt(var(x, na.rm = TRUE) / length(na.omit(x)))
}
#' @export
# Calculate standard error
se <- function(x) {
  sqrt(var(x, na.rm = TRUE) / length(na.omit(x)))
}
