#' Forecasting methods list
#' A list of the forecasting methods for use in the metalearnig process
#' The list follows the format described in the parameter \code{method_list}
#' of \code{\link{process_forecast_dataset}}
#' @seealso \code{\link{process_forecast_dataset}}
#' @export
create_forec_method_list <- function() {
  methods_list <- list("auto_arima_forec")
  methods_list <- append(methods_list, "ets_forec")
  methods_list <- append(methods_list, "nnetar_forec")
  methods_list <- append(methods_list, "tbats_forec")
  methods_list <- append(methods_list, "stlm_ar_forec")
  methods_list <- append(methods_list, "rw_drift_forec")
  methods_list <- append(methods_list, "naive_forec")
  methods_list <- append(methods_list, "snaive_forec")
  methods_list
}

#' @describeIn create_forec_method_list forecast::snaive
#' @param x A \code{ts} object with the input time series
#' @param h The amount of future time steps to forecast
#' @export
snaive_forec <- function(x,h) {
  model <- forecast::snaive(x, h=length(x))
  forecast::forecast(model, h=h)$mean
}

#' @describeIn create_forec_method_list forecast::naive
#' @export
naive_forec <- function(x,h) {
  model <- forecast::naive(x, h=length(x))
  forecast::forecast(model, h=h)$mean
}

#' @describeIn create_forec_method_list forecast::auto.arima
#' @export
auto_arima_forec <- function(x, h) {
  model <- forecast::auto.arima(x, stepwise=FALSE, approximation=FALSE)
  forecast::forecast(model, h=h)$mean
}

#' @describeIn create_forec_method_list forecast::ets
#' @export create_forec_method_list
ets_forec <- function(x, h) {
  model <- forecast::ets(x)
  forecast::forecast(model, h=h)$mean
}

#' @describeIn create_forec_method_list forecast::nnetar
#' @export
nnetar_forec <- function(x, h) {
  model <- forecast::nnetar(x)
  forecast::forecast(model, h=h)$mean
}

#' @describeIn create_forec_method_list forecast::tbats
#' @export
tbats_forec <- function(x, h) {
  model <- forecast::tbats(x)
  forecast::forecast(model, h=h)$mean
}


#' @describeIn create_forec_method_list forecast::stlm with ar modelfunction
#' @export
stlm_ar_forec <- function(x, h) {
  model <- tryCatch({
    forecast::stlm(x, modelfunction = stats::ar)
  }, error = function(e) forecast::auto.arima(x, d=0,D=0))
  forecast::forecast(model, h=h)$mean
}

#' @describeIn create_forec_method_list forecast::rwf
#' @export
rw_drift_forec <- function(x, h) {
  model <- forecast::rwf(x, drift=TRUE, h=length(x))
  forecast::forecast(model, h=h)$mean
}


