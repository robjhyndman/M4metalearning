check_customxgboost_version <- function() {
  #check if the custom xgboost version is installed
  if ( !requireNamespace("xgboost", quietly = TRUE)
       || (utils::packageVersion("xgboost") != '666.6.4.1') ) {
    stop("Xgboost CUSTOM version is required!")

    #warning("Installing it from github pmontman/customxgboost")
    #devtools::install_github("pmontman/customxgboost")
  }
}

# FIX FOR THE CUSTOM MULTICLASS OBJECTIVE : https://github.com/dmlc/xgboost/issues/2776


.onLoad <- function(libname, pkgname) {
  message("Loading metalearning: Checking custom xgboost version ")
  check_customxgboost_version()
}


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

#' Create a classification-like problem from a forecasting-processed time series dataset
#'
#' @param dataset A list with each element having a \code{features} and a optionally \code{errors} fields.
#'     See \code{THA_feature} and \code{process_forecasts} for more information.
#'
#'@return \code{create_feat_classif_problem} returns a list with the entries:
#' \describe{
#'   \item{data}{The features extracted from the series.}
#'   \item{errors}{The errors produced by the forecasting method. If the air availabe in the dataset.}}
#' @export
create_feat_classif_problem <- function(dataset) {
  stopifnot("features" %in% names(dataset[[1]]))

  data <- t(sapply(dataset, function (lentry) as.numeric(lentry$features)))

  if ("errors" %in% names(dataset[[1]]) ) {
    errors <-  t(sapply(dataset, function (lentry) as.numeric(lentry$errors)))
    return(list(data=data, errors=errors))
  } else {
    return(list(data=data))
  }

}

#' Train a method-selecting ensemble that minimizes forecasting error
#'
#' @param data A matrix with the input features data (extracted from the series).
#'     One observation (the features from the original series) per row.
#' @param errors A matrix with the errors produced by each of the forecasting methods.
#'     Each row is a vector with the errors of the forecasting methods.
#' @param params A list containing thr speficic parameters to be passed to the xgboost::xgb.train function
#' @param nrounds nrounds param in xgb.train
#'
#' @export
train_selection_ensemble <- function(data, errors, params, nrounds) {

  dtrain <- xgboost::xgb.DMatrix(data)
  attr(dtrain, "errors") <- errors

  params$objective <- error_softmax_obj
  params$silent <- 0
  params$num_class= length(errors[1,])

  bst <- xgboost::xgb.train(params=params, dtrain, nrounds)
  bst
}



#' @describeIn train_selection_ensemble Produces predictions probabilities for the selection ensemble.
#' @param model The xgboost model
#' @param newdata The feature matrix, one row per series
#' @export
predict_selection_ensemble <- function(model, newdata) {
  pred <- stats::predict(model, newdata, outputmargin = TRUE, reshape=TRUE)
  pred <- t(apply( pred, 1, softmax_transform))
  pred
}



#on one hand we should have the model as first argument following the "predict" R convention
#but this is also a function to be used with *apply type functions
#' @export
predict_weights_meta <- function(seriesentry, model) {
  pred <- stats::predict(model, as.matrix(seriesentry$features), outputmargin = TRUE, reshape=TRUE)
  pred <- t(apply( pred, 1, softmax_transform))
  seriesentry$weights <- pred
  seriesentry
}

#' @export
ensemble_meta <- function(seriesentry) {
  seriesentry$ff_meta_sel <- seriesentry$ff[ which.max(seriesentry$weights), ]
  seriesentry$ff_meta_avg <- seriesentry$weights %*% seriesentry$ff
  seriesentry
}




#' @describeIn train_selection_ensemble Analysis of the predictions, the weighted error and the selection error, along with extra information
#' @param predictions A NXM matrix with N the number of observations(time series)
#'                    and M the number of methods. Each row contains the weights assigned to the methods for the series
#' @param dataset The list with the meta information, forecasts of each method...MUST contain the precalculated errors with process_errors()!
#' @param print.summary Boolean indicating wheter to print the information
#' @param use.precalc.naive2 Boolean indicating wheter the naive2 errors are already contained in \code{dataset} for skip its calculation
#'
#' @export
summary_performance <- function(predictions, dataset, print.summary = TRUE) {
  stopifnot("xx" %in% names(dataset[[1]]))

  #requires precalculated average errors of the dataset
    if (is.null(attr(dataset, "avg_naive2_errors") )) {
      stop("summary_performance requires precalculated avg_naive2_errors, please run process_OWA_errors on dataset")
    }

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

    dataset <- process_owa_errors(dataset)



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
    print(paste("Selected OWA : ", round(selected_error,4)))
    if (!is.null(weighted_error)) {
      print(paste("Weighted OWA : ", round(weighted_error,4)))
      print(paste("Naive Weighted OWA : ", round(naive_weight_error,4)))
    }
    print(paste("Oracle OWA: ", round(oracle_error,4)))
    print(paste("Single method OWA: ", round(single_error,3)))
    print(paste("Average OWA: ", round(average_error,3)))
  }

  list(selected_error = selected_error,
       weighted_error = weighted_error)
}


#' @export
summary_meta <- function(dataset, print.summary = TRUE) {
  stopifnot("xx" %in% names(dataset[[1]]))
  stopifnot("ff" %in% names(dataset[[1]]))
  stopifnot("weights" %in% names(dataset[[1]]))

  weights <- t(sapply(dataset, function (ll) ll$weights))
  summary_performance(weights, dataset, print.summary)
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
#' @import forecast
temp_holdout <- function(dataset) {
  lapply(dataset, temporal_holdout)
}
#' @export
temporal_holdout <- function(seriesentry) {
  frq <- stats::frequency(seriesentry$x)
  if (length(seriesentry$x) - seriesentry$h < max(2 * frq +1, 7)) {
    length_to_keep <- max(2 * stats::frequency(seriesentry$x) +1, 7)
    seriesentry$h <- length(seriesentry$x) - length_to_keep
    if (seriesentry$h < 2) {
      warning( paste( "cannot subset series by",
                      2 - seriesentry$h,
                      " observations, adding a mean constant") )

      seriesentry$x <- stats::ts(c(seriesentry$x, rep(mean(seriesentry$x),2 - seriesentry$h )),
                                 frequency = frq,
                                 start=stats::start(seriesentry$x))
      seriesentry$h <- 2
    }
  }
  #note: we get first the tail, if we subset first, problems will arise (a temp variable for x should be used)
  stopifnot(class(seriesentry$x)=="ts")
  seriesentry$xx <- subset(seriesentry$x, subset=NULL, start = length(seriesentry$x) - seriesentry$h + 1)
  seriesentry$x <- subset(seriesentry$x, subset=NULL, end = length(seriesentry$x) - seriesentry$h )
  if (!is.null(seriesentry$n)) {
    seriesentry$n <- length(seriesentry$x)
  }
  seriesentry
}

