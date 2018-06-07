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


##-#'   \item{xx}{A time series of length \code{h} with the true future data.}
##-#
##-#'   \item{errors}{A vector of F elements containing the OWI errors produced by each of the
##-#'   methods in \code{methods}}
##-#'   }

#' @export
calc_errors <- function(dataset) {

  total_snaive_errors <- c(0,0)
  for (i in 1:length(dataset)) {
    tryCatch({
    lentry <- dataset[[i]]
    insample <- lentry$x

    #extrac forecasts and attach the snaive for completion
    ff <- lentry$ff
    ff <- rbind(ff, snaive_forec(insample, lentry$h))

    frq <- frq <- stats::frequency(insample)
    insample <- as.numeric(insample)
    outsample <- as.numeric(lentry$xx)
    masep <- mean(abs(utils::head(insample,-frq) - utils::tail(insample,-frq)))


    repoutsample <- matrix(
      rep(outsample, each=nrow(ff)),
      nrow=nrow(ff))

    smape_err <- 200*abs(ff - repoutsample) / (abs(ff) + abs(repoutsample))

    mase_err <- abs(ff - repoutsample) / masep

    lentry$snaive_mase <- mase_err[nrow(mase_err), ]
    lentry$snaive_smape <- smape_err[nrow(smape_err),]

    lentry$mase_err <- mase_err[-nrow(mase_err),]
    lentry$smape_err <- smape_err[-nrow(smape_err),]
    dataset[[i]] <- lentry
    total_snaive_errors <- total_snaive_errors + c(mean(lentry$snaive_mase),
                                                   mean(lentry$snaive_smape))
    } , error = function (e) {
      print(paste("Error when processing OWIs in series: ", i))
      print(e)
      e
    })
  }
  total_snaive_errors = total_snaive_errors / length(dataset)
  avg_snaive_errors <- list(avg_mase=total_snaive_errors[1],
                            avg_smape=total_snaive_errors[2])


  for (i in 1:length(dataset)) {
    lentry <- dataset[[i]]
    dataset[[i]]$errors <- 0.5*(rowMeans(lentry$mase_err)/avg_snaive_errors$avg_mase +
                                  rowMeans(lentry$smape_err)/avg_snaive_errors$avg_smape)
    #dataset[[i]]$errors <- rowMeans(lentry$smape_err)
  }
  attr(dataset, "avg_snaive_errors") <- avg_snaive_errors
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


#processes forecast methods on a series
#given a series component
#and list of forecast methods
process_forecast_methods <- function(seriesdata, methods_list) {

  #process each method in methods_list to produce the forecasts and the errors
  lapply(methods_list, function (mentry) {
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
    list( forecasts=forecasts, method_name=method_name)
  })
}


#' Generate Forecasts for a Time Series Dataset
#'
#' For each series in \code{dataset}, forecasts
#' are generated for all methods in \code{methods}.
#'
#' \code{dataset} must be a list with each element having the following format:
#' \describe{
#'   \item{x}{A time series object \code{ts} with the historical data.}
#'   \item{h}{The number of required forecasts.}
#' }
#'
#' \code{methods} is a list of strings with the names of the functions that generate the
#' forecasts. The functions must exist and take as parameters (\code{x}, \code{h}), with
#' \code{x} being the \code{ts} object with the input series and \code{h} the number of required
#' forecasts (after the last observation of \code{x}). The output of these functions must be
#' a vector or \code{ts} object of length \code{h} with the produced forecast.
#' No additional parameters are required in the functions.
#'
#' @param dataset The list containing the series. See details for the required format.
#' @param methods A list of strings with the names of the functions that generate
#' the forecasts.
#' @param n.cores The number of cores to be used. \code{n.cores > 1} means parallel processing.
#'
#' @return A list with the elements having the following structure
#' \describe{
#'   \item{x}{A time series object \code{ts} with the historical data.}
#'   \item{h}{The number of required forecasts.}

#'   \item{ff}{A matrix with F rows and \code{h} columns. Each row contains
#'   the forecasts of each method in \code{methods} }
#' }
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
#'   methods <- list("auto_arima_forec")
#'   methods <- append(methods, "snaive_forec")
#'   methods <- append(methods, "rw_drift_forec")
#'   methods
#' }
#' methods <- create_example_list()
#' forec_results <- calc_forecasts(Mcomp::M3[1:4], methods, n.cores=1)
#'
#' @export
calc_forecasts <- function(dataset, methods, n.cores=1) {
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
    results <- process_forecast_methods(seriesdata, methods)
    ff <- t(sapply(results, function (resentry) resentry$forecasts))
    method_names <- sapply(results, function (resentry) resentry$method_name)
    row.names(ff) <- method_names
    seriesdata$ff <- ff
    seriesdata
  })

  if (n.cores > 1) {
    parallel::stopCluster(cl)
  }

  ret_list
}

