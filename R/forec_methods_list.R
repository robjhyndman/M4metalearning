

#' @export
snaive_forec <- function(x,h) {
  model <- forecast::snaive(x, h=length(x))
  forecast::forecast(model, h=h)$mean
}

#' @export
auto_arima_forec <- function(x, h) {
  model <- forecast::auto.arima(x, stepwise=FALSE, approximation=FALSE)
  forecast::forecast(model, h=h)$mean
}

#' @export
ets_forec <- function(x, h) {
  model <- forecast::ets(x)
  forecast::forecast(model, h=h)$mean
}

#' @export
nnetar_forec <- function(x, h) {
  model <- forecast::nnetar(x)
  forecast::forecast(model, h=h)$mean
}

#' @export
tbats_forec <- function(x, h) {
  model <- forecast::tbats(x)
  forecast::forecast(model, h=h)$mean
}


#' @export
stlm_ar_forec <- function(x, h) {
  model <- tryCatch({
    forecast::stlm(x, modelfunction = stats::ar)
  }, error = function(e) forecast::auto.arima(x, d=0,D=0))
  forecast::forecast(model, h=h)$mean
}

#' @export
rw_drift_forec <- function(x, h) {
  model <- forecast::rwf(x, drift=TRUE, h=length(x))
  forecast::forecast(model, h=h)$mean
}

#' @export
create_seas_method_list <- function() {
  methods_list <- list("auto_arima_forec")
  methods_list <- append(methods_list, "ets_forec")
  methods_list <- append(methods_list, "nnetar_forec")
  methods_list <- append(methods_list, "tbats_forec")
  methods_list <- append(methods_list, "stlm_ar_forec")
  methods_list <- append(methods_list, "rw_drift_forec")
  methods_list <- append(methods_list, "snaive_forec")
}
