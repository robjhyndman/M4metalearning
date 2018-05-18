#' Calculate features from Talagala, Hyndman, Athanaspoulos and add them to the dataset
#'
#'  For each series in \code{dataset}, the feature set used in
#'   (Talagala, Hyndman and Athanasopoulos, 2018) is calculated
#'    and added to the \code{forec_err_dataset}. A \code{tibble} named \code{THA_features}
#'     is added.
#' @param dataset A list the elements having a \code{ts} object with the name \code{x}
#' @param n.cores The number of cores to be used. \code{n.cores > 1} means parallel processing.
#'
#' @examples
#' processed <- generate_THA_feature_dataset(Mcomp::M3[c(1:3)], n.cores=1)
#' processed[[1]]$THA_features
#'
#' @export
generate_THA_feature_dataset <-
  function(dataset, n.cores=1) {
    list_process_fun <- lapply
    cl = -1
    require(tsfeatures)
    if (n.cores > 1) {
      cl <- parallel::makeCluster(n.cores)
      parallel::clusterExport(cl, varlist="dataset", envir=environment())
      parallel::clusterCall(cl, function() library(tsfeatures)) #required to find functions within tsfeatures
      list_process_fun <- function(my_list, ...) {
        parallel::parLapplyLB(cl, my_list, ...)
      }
    }

    dataset_feat <- list_process_fun(dataset,
                               function (serdat) {
                                 tryCatch({
                                   #additional features from Talagala, Hyndman, Athanasopoulos 2018
                                   featrow <-
                                     tsfeatures::tsfeatures(
                                       serdat$x,
                                       features = c(
                                         "acf_features",
                                         "arch_stat",
                                         "crossing_points",
                                         "entropy",
                                         "flat_spots",
                                         "heterogeneity",
                                         "holt_parameters",
                                         "hurst",
                                         "lumpiness",
                                         "nonlinearity",
                                         "pacf_features",
                                         "stl_features",
                                         "stability",
                                         "hw_parameters",
                                         "unitroot_kpss",
                                         "unitroot_pp"
                                       )
                                     )


                                   #additional features
                                   series_length <- length(serdat$x)

                                   featrow <- tibble::add_column(
                                     featrow,
                                     "series_length" = series_length)

                                   featrow[is.na(featrow)] <-
                                     0 #SET NAs TO 0 ?


                                   #adding dummy variables for non seasonal series
                                   #that are not output by tsfeatures
                                    if (length(featrow) == 37) {
                                      featrow <- tibble::add_column(featrow, "seas_acf1" = 0, .before = 7)
                                      featrow <- tibble::add_column(featrow, "seas_pacf" =
                                                                     0, .before = 24)
                                      featrow = tibble::add_column(
                                        featrow,
                                      "seasonal_strength" = 0,
                                      "peak" = 0,
                                      "trough" = 0,
                                      .before=33)
                                    }
                                    serdat$THA_features <- featrow
                                    serdat
                                 }, error = function(e) {
                                   print(e)
                                   return(e)
                                 })
                               })

    if (n.cores > 1) {
      parallel::stopCluster(cl)
    }

    dataset_feat
  }


