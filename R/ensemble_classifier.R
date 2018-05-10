# FIX FOR THE CUSTOM MULTICLASS OBJECTIVE : https://github.com/dmlc/xgboost/issues/2776



#' Softmax Transform
#' @param x A numeric vector.
#' @export
softmax_transform <- function(x) {
  exp(x) / sum(exp(x))
}

# user define objective function, given prediction, return gradient and second order gradient
error_softmax_obj <- function(preds, dtrain) {
  labels <- xgboost::getinfo(dtrain, "label")
  errors <- attr(dtrain, "errors")

  preds <- exp(preds)
  sp <- rowSums(preds)
  preds <- preds / replicate(ncol(preds), sp)
  rowsumerrors <- replicate( ncol(preds), rowSums(preds * errors))

  grad <- preds*(errors - rowsumerrors)
  hess <- errors*preds*(1.0-preds) - grad*preds
  #the true hessian should be grad*(1.0 - 2.0*preds) but it produces numerical problems
  #what we use here is a upper bound

  #print(mean(rowSums(preds*errors)))
  return(list(grad = t(grad), hess = t(hess)))
}




#prepare the time series dataset with extracted features and pose it as a
#custom classification problem

#TO-DO: when draws in the errors, which class to pick?

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

#' Train a method-selecting ensemble that minimizes forecasting error
#'
#' @param data A matrix with the input features data (extracted from the series).
#'     One observation (the features from the original series) per row.
#' @param errors A matrix with the errors produced by each of the forecasting methods.
#'     Each row is a vector with the errors of the forecasting methods.
#' @param labels A numeric vector from 0 to (nclass -1) with the targe labels for classification.
#'     ACTUALLY, THIS IS ALMOST IGNORED, MAY DISSAPEAR IN THE FUTURE.
#'
#' @export
train_selection_ensemble <- function(data, errors, labels) {

  #check if the custom xgboost version is installed
  if ( !requireNamespace("xgboost", quietly = TRUE)
       || (utils::packageVersion("xgboost") != '666.6.4.1') ) {
    warning("Xgboost CUSTOM version is required!")
    warning("Installing it from github pmontman/customxgboost")
    devtools::install_github("pmontman/customxgboost")
  }

  dtrain <- xgboost::xgb.DMatrix(data,
                        label = labels)
  attr(dtrain, "errors") <- errors

  param <- list(max_depth=10, eta=0.1, nthread = 2, silent=1,
                objective=error_softmax_obj,
                num_class=ncol(errors),
                subsample=0.6,
                colsample_bytree=0.6)

  bst <- xgboost::xgb.train(param, dtrain, 100)
  bst
}

#' @describeIn train_selection_ensemble Produces predictions probabilities for the selection ensemble.
#' @export
predict_selection_ensemble <- function(model, newdata) {
  pred <- stats::predict(model, newdata, outputmargin = TRUE, reshape=TRUE)
  pred <- t(apply( pred, 1, softmax_transform))
  pred
}


#' @describeIn train_selection_ensemble Analysis of the predictions
#' @export
summary_performance <- function(predictions, errors, labels, dataset=NULL, print.summary = TRUE) {

  max_predictions <- apply(predictions, 1, which.max) - 1
  class_error <- 1 - mean(max_predictions == labels)
  selected_error <- mean( sapply(1:nrow(errors),
                            function (x) errors[x,max_predictions[x] + 1]) )
  oracle_error <- mean( sapply(1:nrow(errors),
                               function (x) errors[x,labels[x] + 1]) )
  single_error <- min(colMeans(errors))
  average_error <- mean(errors)

  #calculate the weighted prediction
  weighted_error <- NULL
  if (!is.null(dataset)) {
    if ("snaive_forec" %in% rownames(dataset[[1]]$ff) ) {
      snaive_index <- which("snaive_forec" == rownames(dataset[[1]]$ff))
      weighted_error <- sapply(1:nrow(errors), function (i) {
        weighted_forecast <- t(predictions[i,]) %*% dataset[[i]]$ff
        snaive_errors <- calculate_errors(dataset[[i]]$x, dataset[[i]]$xx,
                                          dataset[[i]]$ff[snaive_index,])
        calculate_owi(dataset[[i]]$x, dataset[[i]]$xx,
                      snaive_errors,
                      weighted_forecast)
      })
      weighted_error <- mean(weighted_error)
    }
  }

  if (print.summary) {
    print(paste("Classification error: ", round(class_error,4)))
    print(paste("Selected OWI : ", round(selected_error,4)))
    if (!is.null(weighted_error)) {
      print(paste("Weighted OWI : ", round(weighted_error,4)))
    }
    print(paste("Oracle OWI: ", round(oracle_error,4)))
    print(paste("Single method OWI: ", round(single_error,3)))
    print(paste("Average OWI: ", round(average_error,3)))

  }
}


#' @export
create_tempcv_dataset <- function(dataset) {
  lapply(dataset, function(seriesentry) {
    if (length(seriesentry$x) - seriesentry$h < max(2 * stats::frequency(seriesentry$x) +1, 7)) {
      length_to_keep <- max(2 * stats::frequency(seriesentry$x) +1, 7)
      seriesentry$h <- length(seriesentry$x) - length_to_keep
      if (seriesentry$h < 2) {
        warning( paste( "cannot subset series by",
                        2 - seriesentry$h,
                        " observations, adding a mean constant") )
        seriesentry$x <- stats::ts(c(seriesentry$x, rep(mean(seriesentry$x),2 - seriesentry$h )),
                          frequency = stats::frequency(seriesentry$x))
      }
    }
    #note: we get first the tail, if we subset first, problems will arise (a temp variable for x should be used)
    seriesentry$xx <- utils::tail(seriesentry$x, seriesentry$h)
    seriesentry$x <- utils::head(seriesentry$x, -seriesentry$h)
    if (!is.null(seriesentry$n)) {
      seriesentry$n <- length(seriesentry$x)
    }
    seriesentry
  })
}

if (0) {
train_data <- create_feat_classif_problem(feat_forec_M3)
test_data <- create_feat_classif_problem(feat_forec_M1)

require(xgboost)

dtrain <- xgboost::xgb.DMatrix(train_data$data,
                      label = train_data$labels)
attr(dtrain, "errors") <- train_data$errors
dtest <- xgboost::xgb.DMatrix(test_data$data,
                      label = test_data$labels)
attr(dtest, "errors") <- test_data$errors


num_round <- 2000

param <- list(max_depth=10, eta=0.1, nthread = 2, silent=1,
              objective=error_softmax_obj,
              num_class=ncol(train_data$errors),
              subsample=0.3,
              colsample_bytree=0.3)

bst <- xgboost::xgb.train(param, dtrain, 100)

a <- predict(bst, train_data$data)

pred <- predict(bst, train_data$data, outputmargin = TRUE, reshape=TRUE)
pred = apply(pred, 1, which.max)
mean((pred - 1) == train_data$labels)
owi_error <- mean( sapply(1:nrow(train_data$errors),
                          function (x) train_data$errors[x,pred[x] ]) )
owi_error

#analysis of the train
pred <- predict(bst, test_data$data, outputmargin = TRUE, reshape=TRUE)
pred = apply(pred, 1, which.max) - 1
class_error <- 1 - mean(pred == as.factor(xgboost::getinfo(dtest, "label")))
owi_error <- mean( sapply(1:nrow(test_data$errors),
                    function (x) test_data$errors[x,pred[x] + 1]) )
oracle_error <- mean( sapply(1:nrow(test_data$errors),
                       function (x) test_data$errors[x,test_data$labels[x] + 1]) )
print(paste("Classification error: ", round(class_error,4)))
print(paste("OWI_error: ", round(owi_error,4)))
print(paste("oracle_error: ", round(oracle_error,4)))
print(paste("Single method OWI: ", round(min(colMeans(test_data$errors)),3)))
print(paste("Average OWI: ", round(mean(test_data$errors),3)))

importance <- xgb.importance(model = bst)
head(importance)

forest <- randomForest::randomForest(x = as.matrix(train_data$data),
                                                   y = as.factor(train_data$labels))

testm <- as.matrix(test_data$data)
colnames(testm) <- NULL
pred <- predict(forest, newdata = testm)

1 - mean(pred == as.factor(getinfo(dtrain, "label")))
mean( sapply(1:nrow(test_data$errors),
             function (x) test_data$errors[x,as.numeric(pred[x])]) )
}
