#calculate the MASE and SMAPE errors for all forecasts in serielement$ff
#' @export
calc_mase_smape_errors <- function(serielement) {

  insample <- serielement$x

  ff <- serielement$ff

  frq <- frq <- stats::frequency(insample)
  insample <- as.numeric(insample)
  outsample <- as.numeric(serielement$xx)
  masep <- mean(abs(utils::head(insample,-frq) - utils::tail(insample,-frq)))


  repoutsample <- matrix(
    rep(outsample, each=nrow(ff)),
    nrow=nrow(ff))

  smape_err <- 200*abs(ff - repoutsample) / (abs(ff) + abs(repoutsample))

  if (anyNA(smape_err)) {
    warning("Invalid values when calculating SMAPE error, setting them to the mean error")
    smape_err[is.na(smape_err)] <- mean(smape_err, na.rm=TRUE)
  }
  mase_err <- abs(ff - repoutsample) / masep

  serielement$smape_err <- smape_err
  serielement$mase_err <- mase_err
  serielement
}

#this function goes over the dataset
#' @export
process_naive2_avg_errors <- function(dataset) {
  stopifnot("xx" %in% names(dataset[[1]]))
  stopifnot("x" %in% names(dataset[[1]]))

  naive2_errors <- rowMeans(sapply(dataset, function(lentry) {
    lentry <- calc_forecasts(lentry, list("naive2_forec"))
    lentry <- calc_mase_smape_errors(lentry)
    c( mean(lentry$mase_err) , mean(lentry$smape_err))
  }))

  if (anyNA(naive2_errors)) {
    stop(paste("Invalid values when calculating naive2 errors"))
  }

  avg_naive2_errors <- list(avg_mase=naive2_errors[1],
                            avg_smape=naive2_errors[2])
  attr(dataset, "avg_naive2_errors") <- avg_naive2_errors
  dataset
}

#' @export
process_owa_errors <- function(dataset) {
  stopifnot("xx" %in% names(dataset[[1]]))
  stopifnot("x" %in% names(dataset[[1]]))
  stopifnot("mase_err" %in% names(dataset[[1]]))
  stopifnot("smape_err" %in% names(dataset[[1]]))

  avg_naive2_errors <- attr(dataset, "avg_naive2_errors")
  if (is.null(avg_naive2_errors)) {
    message("Calulating Average Naive2 errors for OWA...")
    dataset <- process_naive2_avg_errors(dataset)
    avg_naive2_errors <- attr(dataset, "avg_naive2_errors")
  }

  for (i in 1:length(dataset)) {
    dataset[[i]] <- calc_mase_smape_errors(dataset[[i]])
  }

  for (i in 1:length(dataset)) {
    dataset[[i]]$errors <- 0.5*(rowMeans(dataset[[i]]$mase_err)/avg_naive2_errors$avg_mase +
                        rowMeans(dataset[[i]]$smape_err)/avg_naive2_errors$avg_smape)
  }

  attr(dataset, "avg_naive2_errors") <- avg_naive2_errors
  dataset
}


############################################################################
####### CALULATE THE FORECASTS AND ERRORS FOR A GIVEN FORECAST METHOD ######
############################################################################

# Given a list of methods, it is applied to the element in the input dataset
# and the element of the output dataset is generated
# the list of R functions is used for easy application of many methods
# the name of the methods in the output entry forec.methods is taken from the functions
# as strings.
# These R functions should take as input a parameter x (the time series) and
# h (the forecast horizon)
# and output only a vector of h elements with the forecasts.
# This way any method in the forecast package can be easily added to the list, and also other
#custom methods



#processes forecast methods on a series
#given a series component (the series and the required horizon)
#and the list of forecast methods to apply
#' @export
calc_forecasts <- function(seriesdata, methods_list) {

  #process each method in methods_list to produce the forecasts and the errors
  ff <- t(sapply(methods_list, function (mentry) {
    method_name <- mentry
    method_fun <- get(mentry)
    forecasts <- tryCatch( method_fun(x=seriesdata$x, h=seriesdata$h),
                           error=function(error) {
                             print(error)
                             print(paste("ERROR processing series: ", seriesdata$st))
                             print(paste("The forecast method that produced the error is:",
                                         method_name))
                             print("Returning snaive forecasts instead")
                             snaive_forec(seriesdata$x, seriesdata$h)
                           })
  }))
  rownames(ff) <- unlist(methods_list)
  seriesdata$ff <- ff
  seriesdata
}
