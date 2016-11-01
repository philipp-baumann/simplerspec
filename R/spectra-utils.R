#' @title Calculate model statistics
#' @description Calculates model statistics for predicted (y)
#' vs. observed (y) values
#' @param df data.frame with predicted and observed data
#' @param x column with observed values
#' @param y column with predicted values
#' @export
summary_df <- function(df, x, y){
  x <- df[, x]
  y <- df[, y]
  data.frame(rmse = sqrt(sum((x - y)^2, na.rm = T) / (length(x)-1)),
    rmsd = mean((x - y)^2)^.5,
    msd = mean((x - y)^2),
    sdev = sd(x, na.rm = T),
    rpd =  sd(x,na.rm = T) /
      sqrt(sum((x - y)^2, na.rm = T) / (length(x) - 1)),
    rpiq = (quantile(x, .75, na.rm = T) - quantile(x, .25, na.rm = T)) /
      sqrt(sum((x - y)^2, na.rm = T) / (length(x) - 1)),
    r2  = cor(x, y, use = "pairwise.complete.obs")^2,
    bias  = mean(x, na.rm = T) - mean(y, na.rm = T),
    SB = (mean(x, na.rm = T) - mean(y, na.rm = T))^2,
    NU = var(x, na.rm = T) * (1 - lm(y ~ x)$coefficients[2])^2,
    LC = var(y, na.rm = T) *
      (1 - cor(x, y, use = "pairwise.complete.obs")^2),
    SB_prop = (mean(x, na.rm = T) - mean(y, na.rm = T))^2 / mean((x - y)^2),
    NU_prop = var(x, na.rm = T) * (1 - lm(y ~ x)$coefficients[2])^2 / mean((x - y)^2),
    LC_prop = var(y, na.rm = T) *
      (1 - cor(x, y, use = "pairwise.complete.obs")^2) / mean((x - y)^2),
    n = length(x)
  )
}
