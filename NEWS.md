# simplerspec 0.1.0.9000

* `model_list` can now also be a list of "train" objects and not only a simplerspec `fit_pls()` output. This saves memory when doing predictions
* Omit message "Joining..." by explicitly providing `by` argument in `dplyr::inner_join`
* Ungroup when slicing tibble so that return data frame does not show "Groups:   sample_id [?]" when printing
* subsequent simplerspec processing functions return (list of) data.table's -> consistency
* newly planned json import/export API based on jsonlite package can passes column names as fields_id only for data.frames and not matrices


# simplerspec 0.1.0

* Added a `NEWS.md` file to track changes to the package

# simplerspec 0.1.0.1

* `read_opus_bin_univ()`: Add support for Bruker files that have undefined `PLF` value (`:= NULL`)
