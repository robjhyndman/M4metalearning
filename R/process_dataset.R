##########################################################
## Given a dataset with the correct format              ##
## forecast are applied to all series in the dataset    ##
## given their respective forecast horizons             ##
## and a new dataset is created with based on the input ##
## storing the individual forecasts, forecast errors,   ##
## and which method is best                             ##
##########################################################


# The input dataset is in the following format (names borrowed from the Mcomp package)
# a list with elements of the following structure
#  x   : The series
#  h   : The number of time steps required to forecast
#  xx  : The true future series (of length h)

# The output dataset will have the input structure, additionally,
# new entries in the structure will be added
# ff              :  A matrix with F rows and h columns, with F being the number of
#                   forecast methods.
# errors          : a vector of F elements containing the errors of the methods





#####################################################################
####### ERROR CALCULATION, TAKEN FROM THE M4 Competition GitHub  ####
####### https://github.com/M4Competition/M4-methods #################
#####################################################################


smape_cal <- function(outsample, forecasts) {
  #Used to estimate sMAPE
  outsample <- as.numeric(outsample) ; forecasts<-as.numeric(forecasts)
  smape <- (abs(outsample-forecasts)*200)/(abs(outsample)+abs(forecasts))
  return(smape)
}

mase_cal <- function(insample, outsample, forecasts) {
  stopifnot(stats::is.ts(insample))
  #Used to estimate MASE
  frq <- stats::frequency(insample)
  forecastsNaiveSD <- rep(NA,frq)
  for (j in (frq+1):length(insample)){
    forecastsNaiveSD <- c(forecastsNaiveSD, insample[j-frq])
  }
  masep<-mean(abs(insample-forecastsNaiveSD),na.rm = TRUE)

  outsample <- as.numeric(outsample) ; forecasts <- as.numeric(forecasts)
  mase <- (abs(outsample-forecasts))/masep
  return(mase)
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



#calculate forecast predictions
calculate_forecast_preds <- function(insample, h, forec.method) {
  forec.method(x=insample, h=h)
}


#calc SMAPE and MASE errors for a given model and parameters
calculate_errors <- function(insample, outsample, forecasts) {
  SMAPE <- smape_cal(outsample, forecasts)
  MASE <- mase_cal(insample, outsample, forecasts)
  c(mean(SMAPE) , mean(MASE))
}

#output only the owi error
calculate_owi <- function(insample, outsample, snaive_errors, forecasts) {
  errors <- calculate_errors(insample, outsample, forecasts)
  0.5*( (errors[1] / snaive_errors[1]) +  (errors[2]/snaive_errors[2]))
}


#processes forecast methods on a series, outputting the error
#given a series component
#and list of forecast methods
process_forecast_methods <- function(seriesdata, methods_list) {

  #calc the snaive errors first, they will be used for the owi error of all the methods
  forecast_snaive <- snaive_forec(seriesdata$x, h=seriesdata$h)
  snaive_errors <- calculate_errors(seriesdata$x, seriesdata$xx, forecast_snaive)

  #clamp snaive errors to a very small number
  snaive_errors[snaive_errors < 0.0001] <- 0.0001


  #process each method in methods_list to produce the forecasts and the errors
  lapply(methods_list, function (mentry) {
    method_name <- mentry
    method_fun <- get(mentry)
    forecasts <- method_fun(x=seriesdata$x, h=seriesdata$h)
    owi_error <- calculate_owi(seriesdata$x, seriesdata$xx, snaive_errors, forecasts)
    list(error=owi_error, forecasts=forecasts, method_name=method_name)
  })
}


#' Generate Forecast and Errors for a time series dataset
#'
#' For each series in \code{dataset}, forecasts and
#' owi errors are generated for all methods in \code{methods_list}.
#'
#' \code{dataset} must be with each element having the following format:
#' \describe{
#'   \item{x}{A time series object \code{ts} with the historical data.}
#'   \item{h}{The number of required forecasts.}
#'   \item{xx}{A time series of length code{h} with the true future data.}
#' }
#' \code{methods_list} is a list of strings with the names of the functions that generate the
#' forecasts. The functions must exist and take as parameters (\code{x}, \code{h}), with
#' \code{x} being the \code{ts} object with the input series and \code{h} the number of required
#' forecasts (after the last observation of \code{x}). The output of these functions must be
#' a vector or \code{ts} object of length \code{h} with the produced forecast.
#' No additional parameters are required in the functions.
#'
#' @param dataset The list containing the series. See details for the required format.
#' @param methods_list A list of strings with the names of the functions that generate
#' the forecasts.
#' @param n.cores The number of cores to be used. \code{n.cores > 1} means parallel processing.
#'
#' @return A list with the elements having the following structure
#' \describe{
#'   \item{x}{A time series object \code{ts} with the historical data.}
#'   \item{h}{The number of required forecasts.}
#'   \item{xx}{A time series of length \code{h} with the true future data.}
#'   \item{ff}{A matrix with F rows and \code{h} columns. Each row contains
#'   the forecasts of each method in \code{methods_list}}
#'   \item{errors}{A vector of F elements containing the OWI errors produced by each of the
#'   methods in \code{methods_list}}
#'   }
#'
#' @examples
#' auto_arima_forec <- function(x, h) {
#'   model <- forecast::auto.arima(x, stepwise=FALSE, approximation=FALSE)
#'   forecast::forecast(model, h=h)$mean
#' }
#'
#' snaive_forec <- function(x,h) {
#'   model <- forecast::snaive(x, h=length(x))
#'   forecast::forecast(model, h=h)$mean
#' }
#' rw_drift_forec <- function(x, h) {
#'   model <- forecast::rwf(x, drift=TRUE, h=length(x))
#'  forecast::forecast(model, h=h)$mean
#' }
#'
#' create_example_list <- function() {
#'   methods_list <- list("auto_arima_forec")
#'   methods_list <- append(methods_list, "snaive_forec")
#'   methods_list <- append(methods_list, "rw_drift_forec")
#'   methods_list
#' }
#' methods <- create_example_list()
#' forec_results <- process_forecast_dataset(Mcomp::M3[1:4], methods, n.cores=1)
#'
#' @export
process_forecast_dataset <- function(dataset, methods_list, n.cores=1) {
  list_process_fun <- lapply
  cl = -1

  if (n.cores > 1) {
    cl <- parallel::makeCluster(n.cores)
    parallel::clusterExport(cl, varlist=ls(), envir=environment())
    parallel::clusterExport(cl, varlist=ls(envir=environment(process_forecast_methods)),
                            envir = environment(process_forecast_methods))
    list_process_fun <- function(my_list, ...) {
      parallel::parLapplyLB(cl, my_list, ...)
    }
  }

  ret_list <- list_process_fun(dataset, function (seriesdata) {
    results <- process_forecast_methods(seriesdata, methods_list)
    ff <- t(sapply(results, function (resentry) resentry$forecasts))
    method_names <- sapply(results, function (resentry) resentry$method_name)
    row.names(ff) <- method_names
    errors <- sapply(results, function (resentry) resentry$error)
    names(errors) <- method_names
    list(x = seriesdata$x,
         xx = seriesdata$xx,
         h = seriesdata$h,
         ff = ff,
         errors = errors)
  })

  if (n.cores > 1) {
    parallel::stopCluster(cl)
  }

  ret_list
}

