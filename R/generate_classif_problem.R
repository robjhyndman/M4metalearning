#' Calculate features from Talagala, Hyndman, Athanaspoulos and add them to the dataset
#'
#'  For each series in \code{dataset}, the feature set used in
#'   (Talagala, Hyndman and Athanasopoulos, 2018) is calculated
#'    and added to the \code{forec_err_dataset}. A \code{tibble} named \code{features}
#'     is added.
#' @param dataset A list the elements having a \code{ts} object with the name \code{x}
#' @param n.cores The number of cores to be used. \code{n.cores > 1} means parallel processing.
#'
#' @examples
#' processed <- generate_THA_features(Mcomp::M3[c(1:3)], n.cores=1)
#' processed[[1]]$features
#'
#' @export
THA_features <-
  function(dataset, n.cores=1) {
    list_process_fun <- lapply
    cl = -1
    require(tsfeatures)
    if (n.cores > 1) {
      cl <- parallel::makeCluster(n.cores)
      #parallel::clusterExport(cl, varlist="dataset", envir=environment())
      parallel::clusterExport(cl, varlist=ls(), envir=environment())
      parallel::clusterExport(cl, varlist=ls(envir=environment(THA_features)),
                              envir = environment(THA_features))
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
                                         "heterogeneity_tsfeat_workaround",
                                         "holt_parameters",
                                         "hurst",
                                         "lumpiness",
                                         "nonlinearity",
                                         "pacf_features",
                                         "stl_features",
                                         "stability",
                                         "hw_parameters_tsfeat_workaround",
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
                                    serdat$features <- featrow
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

#' @export
heterogeneity_tsfeat_workaround <- function(x) {
  output <- c(arch_acf =0, garch_acf=0, arch_r2=0, garch_r2=0)
  try( output <- tsfeatures::heterogeneity(x) )
  output
}

#' @export
hw_parameters_tsfeat_workaround <- function(x) {
  hw_fit <- NULL
  hw_fit$par <- c(NA, NA, NA)
  try(hw_fit <- forecast::ets(x, model=c("AAA")), silent=TRUE)
  names(hw_fit$par) <- c("hw_alpha", "hw_beta" , "hw_gamma")
  hw_fit$par[1:3]
}
