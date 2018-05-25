check_customxgboost_version <- function() {
  #check if the custom xgboost version is installed
  if ( !requireNamespace("xgboost", quietly = TRUE)
       || (utils::packageVersion("xgboost") != '666.6.4.1') ) {
    warning("Xgboost CUSTOM version is required!")
    warning("Installing it from github pmontman/customxgboost")
    devtools::install_github("pmontman/customxgboost")
  }
}

# FIX FOR THE CUSTOM MULTICLASS OBJECTIVE : https://github.com/dmlc/xgboost/issues/2776



#' Softmax Transform
#' @param x A numeric vector.
#' @export
softmax_transform <- function(x) {
  exp(x) / sum(exp(x))
}

# user defined objective function for xgboost
# minimizes de class probabilities * owi_errors
#' @export
error_softmax_obj <- function(preds, dtrain) {
  labels <- xgboost::getinfo(dtrain, "label")
  errors <- attr(dtrain, "errors")

  preds <- exp(preds)
  sp <- rowSums(preds)
  preds <- preds / replicate(ncol(preds), sp)
  rowsumerrors <- replicate( ncol(preds), rowSums(preds * errors))

  grad <- preds*(errors - rowsumerrors)
  hess <- errors*preds*(1.0-preds) - grad*preds
  #hess <- grad*(1.0 - 2.0*preds)
  #hess <- pmax(hess, 1e-16)
  #the true hessian should be grad*(1.0 - 2.0*preds) but it produces numerical problems
  #what we use here is a upper bound

  #print(mean(rowSums(preds*errors)))
  return(list(grad = t(grad), hess = t(hess)))
}



#prepare the time series dataset with extracted features and pose it as a
#custom classification problem
#TO-DO: when there are draws in the errors, which class to pick as label?

#' Create a classification problem from a forecasting-processed time series dataset
#'
#' @param dataset A list with each element having a \code{THA_features} and a \code{errors} fields.
#'     See \code{generate_THA_feature_dataset} and \code{process_forecast_dataset} for more information.
#'
#'@return \code{create_feat_classif_problem} returns a list with the entries:
#' \describe{
#'   \item{data}{The features extracted from the series.}
#'   \item{errors}{The errors produced by the forecasting method.}
#'   \item{labels}{The target classification problem, created by selecting the method that produces.
#'       Integer from 0 to (nmethods-1).}
#'   }
#' @export
create_feat_classif_problem <- function(dataset) {
  stopifnot("THA_features" %in% names(dataset[[1]]))
  extracted <- t(sapply(dataset, function (lentry) {
    seriesdata <- c(as.numeric(lentry$THA_features), which.min(lentry$errors) -1,
      lentry$errors
      )
    names(seriesdata) <- c( names(lentry$THA_features), "best_method", names(lentry$errors))
    seriesdata
  }))

  return_data <- list(data = extracted[, 1:length(dataset[[1]]$THA_features)],
       labels = extracted[, length(dataset[[1]]$THA_features) +1],
       errors = extracted[, -(1:(length(dataset[[1]]$THA_features) +1))]
       )

  return_data
}

#' @describeIn metatemp_train Train a method-selecting ensemble that minimizes forecasting error
#'
#' @param data A matrix with the input features data (extracted from the series).
#'     One observation (the features from the original series) per row.
#' @param errors A matrix with the errors produced by each of the forecasting methods.
#'     Each row is a vector with the errors of the forecasting methods.
#' @param labels A numeric vector from 0 to (nclass -1) with the targe labels for classification.
#'
#' @export
train_selection_ensemble <- function(data, errors, param=NULL) {

  check_customxgboost_version()

  dtrain <- xgboost::xgb.DMatrix(data)
  attr(dtrain, "errors") <- errors

  if (is.null(param)) {
    param <- list(max_depth=10, eta=0.4, nthread = 3, silent=1,
                  objective=error_softmax_obj,
                  num_class=ncol(errors),
                  subsample=0.9,
                  colsample_bytree=0.6
    )
  }

  bst <- xgboost::xgb.train(param, dtrain, 200)
  bst
}

#' @describeIn metatemp_train Produces predictions probabilities for the selection ensemble.
#' @param model The xgboost model
#' @param newdata The feature matrix, one row per series
#' @export
predict_selection_ensemble <- function(model, newdata) {
  pred <- stats::predict(model, newdata, outputmargin = TRUE, reshape=TRUE)
  pred <- t(apply( pred, 1, softmax_transform))
  pred
}


#' @describeIn metatemp_train Analysis of the predictions
#' @param predictions A NXM matrix with N the number of observations(time series)
#'                    and M the number of methods.
#' @param errors The NXM matrix with the erros per series and per method
#' @param labels Integer vector The true labels of the would be classification problem.
#'               Possible values are from 0 to M-1
#' @param dataset The list with the meta information, if given additional details will be provided
#' @param print.summary Boolean indicating wheter to print the information
#'
#' @export
summary_performance <- function(predictions, dataset, print.summary = TRUE) {

  labels <- sapply(dataset, function(lentry) which.min(lentry$errors) - 1)
  max_predictions <- apply(predictions, 1, which.max) - 1
  class_error <- 1 - mean(max_predictions == labels)
  n_methods <- length(dataset[[1]]$errors)

  #selected_error <- mean( sapply(1:nrow(errors),
  #                               function (i) errors[i,max_predictions[i] + 1]) )
  #oracle_error <- mean( sapply(1:nrow(errors),
  #                             function (i) errors[i,labels[i] + 1]) )
  #single_error <- min(colMeans(errors))
  #average_error <- mean(errors)

  #calculate the weighted prediction

    for (i in 1:length(dataset)) {
      weighted_ff <- t(predictions[i,]) %*% dataset[[i]]$ff
      naive_combi_ff <- colMeans(dataset[[i]]$ff)
      selected_ff <- dataset[[i]]$ff[max_predictions[i]+1,]
      oracle_ff <- dataset[[i]]$ff[labels[i]+1,]
      dataset[[i]]$ff <- rbind(dataset[[i]]$ff,
                               weighted_ff,
                               naive_combi_ff,
                               selected_ff,
                               oracle_ff)
    }

    dataset <- fast_errors_dataset(dataset)
    all_errors <- sapply(dataset, function (lentry) {
      lentry$errors})

    all_errors <- rowMeans(all_errors)
    weighted_error <- all_errors[n_methods + 1]
    naive_weight_error  <- all_errors[n_methods + 2]
    selected_error  <- all_errors[n_methods + 3]
    oracle_error <-  all_errors[n_methods + 4]
    single_error <- min(all_errors[1:n_methods])
    average_error <- mean(all_errors[1:n_methods])



  if (print.summary) {
    print(paste("Classification error: ", round(class_error,4)))
    print(paste("Selected OWI : ", round(selected_error,4)))
    if (!is.null(weighted_error)) {
      print(paste("Weighted OWI : ", round(weighted_error,4)))
      print(paste("Naive Weighted OWI : ", round(naive_weight_error,4)))
    }
    print(paste("Oracle OWI: ", round(oracle_error,4)))
    print(paste("Single method OWI: ", round(single_error,3)))
    print(paste("Average OWI: ", round(average_error,3)))

  }
  list(selected_error = selected_error, weighted_error = weighted_error)
}






#' Create Temporal Crossvalidated Dataset
#'
#' Creates a datasets of time series for forecasting and metalearning
#' by extracting the final observations of each series
#'
#' The final \code{h} observation are extrated an posed as the true future values
#' in the entry \code{xx} of the output list. At least 7 observations are kept as the
#' observable time series \code{x}.
#'
#' @param dataset A list with each element having at least the following
#' \describe{
#'   \item{x}{A time series object \code{ts} with the historical data.}
#'   \item{h}{The number of required forecasts.}
#' }
#'
#' @return A list with the same structure as the input,
#'  the following entries may be added or overwritten if existing
#' \describe{
#'   \item{x}{A time series object \code{ts} with the historical data.}
#'   \item{xx}{A time series with the true future data. Has length \code{h}
#'   unless the remaining \code{x} would be too short.}
#' }
#' @export
create_tempcv_dataset <- function(dataset) {
  lapply(dataset, function(seriesentry) {
    frq <- stats::frequency(seriesentry$x)
    if (length(seriesentry$x) - seriesentry$h < max(2 * frq +1, 7)) {
      length_to_keep <- max(2 * stats::frequency(seriesentry$x) +1, 7)
      seriesentry$h <- length(seriesentry$x) - length_to_keep
      if (seriesentry$h < 2) {
        warning( paste( "cannot subset series by",
                        2 - seriesentry$h,
                        " observations, adding a mean constant") )

        seriesentry$x <- stats::ts(c(seriesentry$x, rep(mean(seriesentry$x),2 - seriesentry$h )),
                          frequency = frq)
        seriesentry$h <- 2
      }
    }
    #note: we get first the tail, if we subset first, problems will arise (a temp variable for x should be used)
    seriesentry$xx <- utils::tail(seriesentry$x, seriesentry$h)
    seriesentry$x <- stats::ts( utils::head(seriesentry$x, -seriesentry$h),
                                frequency = frq)
    if (!is.null(seriesentry$n)) {
      seriesentry$n <- length(seriesentry$x)
    }
    seriesentry
  })
}
