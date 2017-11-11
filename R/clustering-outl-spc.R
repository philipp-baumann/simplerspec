## Robust multivariate outlier detection method based on semi-robust principal
## components ==================================================================

# Apply pcout for a nested list of matrices at a certain depth
pcout_depth <- function(x, depth = 2) {
  purrr::modify_depth(.x = x, .depth = depth,
    .f = ~ mvoutlier::pcout(x = ., makeplot = FALSE))
}

# Extracts mvoutlier::pcout $wscat elements in nested list (default at depth 2)
wscat_depth <- function(x, depth = 2) {
  purrr::modify_depth(.x = x, .depth = 2, .f = "wscat")
}

# Returns nested list containing logicals from test which scattering weights
# are zero
which_wscat0_depth <- function(x) {
  purrr::modify_depth(wscat_depth(x), 2, ~ . == 0)
}

# Extract $wfinal01 elements
wfinal01_depth <- function(x) {
  purrr::modify_depth(x, 2, "wfinal01")
}

# Returns nested list containing logicals from test which final 0/1 weights
# are zero
which_wfinal0_depth <- function(x) {
  purrr::modify_depth(wfinal01_depth(x), 2, ~ . == 0)
}
