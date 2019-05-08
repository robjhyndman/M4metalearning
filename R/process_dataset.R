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



#calculate the MASE and SMAPE errors for all forecasts in seriesentry$ff
calc_err <- function(serielement) {

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

  mase_err <- abs(ff - repoutsample) / masep

  serielement$smape_err <- smape_err
  serielement$mase_err <- mase_err
  serielement
}

proc.err <- function(dataset, n.cores=1) {
  (futurlapplyfy(calc_err, n.cores))(dataset)
}

process_errors <- function(dataset, n.cores=1) {

  if (is.null(attr(dataset, "avg_naive2_errors"))) {
    message("Naive2 errors for OWA not found, calculating them...")

    datanaive2 <- process_forecasts(dataset, list(naive2_forec()), n.cores)
    datanaive2 <- proc.err(datanaive2)
    #get the average of the two errors
    err_naive2 <- sapply(datanaive2, function(ll) c(mean(ll$mase_err), mean(ll$smape_err)))
    total_naive2_errors <- colMeans(err_naive2)

    avg_naive2_errors <- list(avg_mase=total_naive2_errors[1],
                              avg_smape=total_naive2_errors[2])

    attr(dataset, "avg_naive2_errors") <- avg_naive2_errors
  }

  dataset <- proc.err(dataset)
  avg_naive2_errors <- attr(dataset, "avg_naive2_errors")
  lapply(dataset, function(ll) {
    ll$errors <- 0.5*(rowMeans(ll$mase_err)/avg_naive2_errors$avg_mase +
           rowMeans(ll$smape_err)/avg_naive2_errors$avg_smape)
  })

}

##-#'   \item{xx}{A time series of length \code{h} with the true future data.}
##-#
##-#'   \item{errors}{A vector of F elements containing the OWI errors produced by each of the
##-#'   methods in \code{methods}}
##-#'   }

#' @export
calc_errors <- function(dataset, use.precalc.naive2 = FALSE) {
  #sanity checking
  stopifnot("xx" %in% names(dataset[[1]]))
  stopifnot("ff" %in% names(dataset[[1]]))
  stopifnot("x" %in% names(dataset[[1]]))

  if (use.precalc.naive2) {
    stopifnot(!is.null(attr(dataset, "avg_naive2_errors") ))
    #stopifnot("naive2_mase" %in% names(dataset[[1]]))
    #stopifnot("naive2_smape" %in% names(dataset[[1]]))
  }

  total_naive2_errors <- c(0,0) #accumulate the errors of the naive2 to get its average over the dataset
  for (i in 1:length(dataset)) {
    tryCatch({
    lentry <- dataset[[i]]
    insample <- lentry$x

    #extrac forecasts and attach the naive2 for completion
    ff <- lentry$ff
    if (!use.precalc.naive2) {
      ff <- rbind(ff, naive2_forec(insample, lentry$h))
    }

    frq <- frq <- stats::frequency(insample)
    insample <- as.numeric(insample)
    outsample <- as.numeric(lentry$xx)
    masep <- mean(abs(utils::head(insample,-frq) - utils::tail(insample,-frq)))


    repoutsample <- matrix(
      rep(outsample, each=nrow(ff)),
      nrow=nrow(ff))

    smape_err <- 200*abs(ff - repoutsample) / (abs(ff) + abs(repoutsample))

    mase_err <- abs(ff - repoutsample) / masep

    if (!use.precalc.naive2) {
      naive2_mase <- mase_err[nrow(mase_err), ]
      naive2_smape <- smape_err[nrow(smape_err),]
      total_naive2_errors <- total_naive2_errors + c(mean(naive2_mase),
                                                     mean(naive2_smape))
      if (anyNA(total_naive2_errors)) {
        stop(paste("Invalid error values when processing series: ",i))
      }
      #remove the naive2 calculation from the output forecasts
      lentry$mase_err <- mase_err[-nrow(mase_err),]
      lentry$smape_err <- smape_err[-nrow(smape_err),]
    } else {
      lentry$mase_err <- mase_err
      lentry$smape_err <- smape_err
    }

    dataset[[i]] <- lentry

    } , error = function (e) {
      print(paste("Error when processing OWAs in series: ", i))
      print(e)
      e
    })
  }

  if (!use.precalc.naive2) {
    total_naive2_errors <- total_naive2_errors / length(dataset)
    avg_naive2_errors <- list(avg_mase=total_naive2_errors[1],
                            avg_smape=total_naive2_errors[2])
  } else {
    avg_naive2_errors <- attr(dataset, "avg_naive2_errors")
  }


  for (i in 1:length(dataset)) {
    lentry <- dataset[[i]]
    dataset[[i]]$errors <- 0.5*(rowMeans(lentry$mase_err)/avg_naive2_errors$avg_mase +
                                  rowMeans(lentry$smape_err)/avg_naive2_errors$avg_smape)
    #dataset[[i]]$errors <- rowMeans(lentry$smape_err)
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


#
# #calculate forecast predictions
# calculate_forecast_preds <- function(insample, h, forec.method) {
#   forec.method(x=insample, h=h)
# }
#
#
# #calc SMAPE and MASE errors for a given model and parameters
# calculate_errors <- function(insample, outsample, forecasts) {
#   SMAPE <- smape_cal(outsample, forecasts)
#   MASE <- mase_cal(insample, outsample, forecasts)
#   c(mean(SMAPE) , mean(MASE))
# }
#
# #output only the owi error
# calculate_owi <- function(insample, outsample, snaive_errors, forecasts) {
#   errors <- calculate_errors(insample, outsample, forecasts)
#   0.5*( (errors[1] / snaive_errors[1]) +  (errors[2]/snaive_errors[2]))
# }


#processes forecast methods on a series
#given a series component (the series and the required horizon)
#and the list of forecast methods to apply
calc_forecast_methods <- function(seriesdata, methods_list) {

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
#' @param forecast_methods A list of strings with the names of the functions that generate
#' the forecasts.
#' @param n.cores The number of cores to be used. \code{n.cores > 1} means parallel processing.
#' #' @param n.cores The number of cores to be used. \code{n.cores > 1} means parallel processing.
#'
#' @return A list with the same structure as \code{dataset} but with the added entry
#' \describe{
#'   \item{ff}{A matrix with F rows and \code{h} columns. Each row contains
#'   the forecasts of each method in \code{methods} }
#' } to each element on a list
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
#' forec_results <- process_forecasts(Mcomp::M3[1:4], methods, n.cores=1)
#'
#' @export
process_forecasts <- function(dataset, forecast_methods, n.cores=1, chunk_size=0, do_shuffle=0,
                              save_checkpoint_filename,
                              load_checkpoint_filename) {

  calcforec_fun <- function (seriesdata) {
    results <- calc_forecast_methods(seriesdata, forecast_methods)
    ff <- do.call("rbind", lapply(results, function (resentry) resentry$forecasts))
    method_names <- sapply(results, function (resentry) resentry$method_name)
    row.names(ff) <- method_names
    seriesdata$ff <- ff
    seriesdata
  }

  parchunk_calcforec <- chunkparfy( calcforec_fun,
                                    n.cores, chunk_size, do_shuffle,
                                    save_checkpoint_filename,
                                    load_checkpoint_filename)

  parchunk_calcforec(dataset)
}




#transform a function that works on a list, into the same function, but with saving/resuming capabilities
#' @export
chunkify <- function( myFUN, chunk_size, do_shuffle=TRUE,
                      save_checkpoint_filename=NULL,
                      load_checkpoint_filename=NULL) {

  function(dataset) {
    #no chunk_size, return the original function
    if (chunk_size==0) {
      return(myFUN(dataset))
    }

    temp_dataset <- NULL

    if (do_shuffle) {
      shuffling <- sample(length(dataset))
    } else {
      shuffling <- 1:length(dataset)
    }


    data_length <- length(dataset)

    chunk_index <- seq(1, data_length, chunk_size)
    chunk_index <- c(chunk_index, data_length+1) #add a final DUMMY chunk
    start_chunk <- 1

    if (!is.null(load_checkpoint_filename)) {
      message(paste("Resuming processing of the file:", load_checkpoint_filename))
      tmp_chunk_size <- chunk_size
      tryCatch( load(load_checkpoint_filename),
                error= function(e) { stop("Resuming filename not found!!");e})
      if (tmp_chunk_size != chunk_size) {
        stop(paste("Function call chunk_size and resumed chunk_size do not coincide, please check!! Resuming chunk_size was: ", chunk_size))
      }
      resume_ind <- length(temp_dataset)+1
      start_chunk <- which(chunk_index == resume_ind)

      if (length(start_chunk) != 1) {
        stop("Error when calculating resuming starting chunk")
      }

      if (start_chunk == length(chunk_index)) {
        message("No more processing is required, returning previous computations")
        return(temp_dataset[shuffling])
      }
    }

    #shuffle the dataset to even the processing cost of each chunk
    #after we have loaded the previous shuffling if we are resuming the computation
    dataset <- dataset[shuffling]


    for (i in start_chunk:(length(chunk_index)-1)) {
      start_ind = chunk_index[i]
      end_ind = chunk_index[i+1]-1
      temp_dataset <- append(temp_dataset, myFUN(dataset[start_ind:end_ind]))
      message(paste("From ", start_ind, " to", end_ind,
                    ", ", round(100*(start_ind-1) / length(dataset),2),
                    "% of the dataset processed"))
      if (!is.null(save_checkpoint_filename)) {
        save(temp_dataset, shuffling, chunk_size, file = save_checkpoint_filename)
      }
    }

    #return the processed dataset
    temp_dataset[order(shuffling)]

  }
}


### get a function that works on an element and
### return a lapply version of the function, but parallel
#' @export
futurlapplyfy <- function(myFUN, n.cores=1) {
  function (dataset) {
    list_process_fun <- lapply
    aaaa <- n.cores
    if (n.cores > 1) {
      oldplan <- future::plan()
      newstrat <- future::tweak(future::multiprocess, workers = n.cores)
      future::plan(newstrat)
      on.exit(future::plan(oldplan), add = TRUE)
      list_process_fun <- function(dataset, ...) {
        future.apply::future_lapply(dataset, ...)
      }
    }
    list_process_fun(dataset, myFUN)
  }
}



chunkparfy <- function(myFUN, n.cores, chunk_size, do_shuffle=TRUE,
                       save_checkpoint_filename=NULL,
                       load_checkpoint_filename=NULL) {

  myparfun <- futurlapplyfy(myFUN, n.cores)
  mychunkparfun <- chunkify(myparfun, chunk_size,
                                            do_shuffle,
                                            save_checkpoint_filename,
                                            load_checkpoint_filename)
  mychunkparfun

}
