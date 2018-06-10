
##
# put here everything to check in the datasets
# same number of methods, and ids of them
# errors
# fields of the list, etc
sanity_check_function <- function() {

}


#squared loss
combi_forec_square_obj <- function(preds, dtrain) {

  ff <- attr(dtrain, "ff")
  xx <- attr(dtrain, "xx")
  ew <- attr(dtrain, "ew")

  SE <-  ew*(rowSums(preds * ff) - xx)
  grad <- 2*ew *ff* SE
  hess <- 2* ew*ew*ff * ff

  print(mean(abs(SE)))
  return(list(grad = t(grad), hess = t(hess)))
}

combi_softmax_square <- function(preds, dtrain) {
  ff <- attr(dtrain, "ff")
  xx <- attr(dtrain, "xx")
  ew <- attr(dtrain, "ew")

  preds <- exp(preds)
  sp <- rowSums(preds)
  preds <- preds / replicate(ncol(preds), sp)

  S <- rowSums(preds * ff)
  Sxx <- ew*(S - xx)
  GradSxx <- ew*preds*(ff - S)
  grad <- 2*Sxx*GradSxx #derivative of squaredSxx
  hess <- 2* GradSxx*( GradSxx + Sxx*(1.0 - 2.0*preds))

  return(list(grad = t(grad), hess = t(hess)))
}

sgm <- function(x) {
  1.0 / (1.0 + exp(-x))
}


sigmoid_clamp = function(x) {
  x <- pmax(pmin(x, 9), -9)
  1.0 / (1.0 + exp(-x))
}

smooth_sign <- function(x) {
  2.0 * sigmoid_clamp(x) - 1
}


combi_forec_absolute_obj <- function(preds, dtrain) {

  ff <- attr(dtrain, "ff")
  xx <- attr(dtrain, "xx")
  ew <- attr(dtrain, "ew")

  SE <-  ew*(rowSums(preds * ff) - xx)
  SS <- smooth_sign(SE)
  grad <- ew*SS*ff
  hess <- ew*ff * 2*sigmoid_clamp(SE)*(1-sigmoid_clamp(SE))*ff*ew

  return(list(grad = t(grad), hess = t(hess)))
}


DEPRECATED_old_combi_softmax_abs <- function(preds, dtrain) {
  ff <- attr(dtrain, "ff")
  xx <- attr(dtrain, "xx")
  ew <- attr(dtrain, "ew")

  preds <- exp(preds)
  sp <- rowSums(preds)
  preds <- preds / replicate(ncol(preds), sp)

  S <- rowSums(preds * ff)
  Sxx <-  ew*(S - xx)
  GradSxx <- ew*(preds*(ff - S))
  grad = smooth_sign(Sxx)*GradSxx
  hess = 2*sigmoid_clamp(Sxx)*(1-sigmoid_clamp(Sxx))*GradSxx #+
   #0* smooth_sign(Sxx)*ew*(ff - preds*(ff-0*S))

  print( mean(abs(Sxx)))

  return(list(grad = t(grad), hess = t(hess)))
}

combi_softmax_abs <- function(preds, dtrain) {
  ff <- attr(dtrain, "ff")
  xx <- attr(dtrain, "xx")
  ew <- attr(dtrain, "ew")
  #lambda <- attr(dtrain, "lambda")

  a = preds

  preds <- exp(preds)
  sp <- rowSums(preds)
  preds <- preds / replicate(ncol(preds), sp)

  S <- rowSums(preds * ff)
  Sxx <-  ew*(S - xx)
  GradSxx <- ew*(preds*(ff - S))
  grad = smooth_sign( Sxx )*GradSxx
  HesSxx = ew*( preds*(1-preds)*ff - preds*(1-preds)*S - preds*GradSxx/ew)

  hess = 2 * sigmoid_clamp(Sxx)*(1 - sigmoid_clamp(Sxx))*GradSxx*GradSxx +
    smooth_sign(Sxx)*( HesSxx)
  hess = pmax(hess, 1e-16)
  #hess[(hess > -1e-13) & (hess < 0)] <- -1e-13
  #hess[(hess < 1e-16) ] <- 1e-16


  #print( mean(abs(Sxx)))

  return(list(grad = t(grad), hess = t(hess)))
}


#' @export
combi_softmax_owi <- function(preds, dtrain) {
  ff <- attr(dtrain, "ff")
  xx <- attr(dtrain, "xx")
  ew <- attr(dtrain, "ew")
  eh <- attr(dtrain, "eh")
  avg_mase <- attr(dtrain, "avg_mase")
  avg_smape <- attr(dtrain, "avg_smape")

  preds <- exp(preds)
  sp <- rowSums(preds)
  preds <- preds / replicate(ncol(preds), sp)

  S <- rowSums(preds * ff)
  Sxx <-  (S - xx)
  GradSxx = preds*(ff - S)
  D = abs(xx) + abs(S)
  GradD = GradSxx   #smooth_sign(S)* GradSxx #always positive

  gradMASE = smooth_sign( Sxx )*GradSxx
  gradSMAPE =  (gradMASE * D  - abs(Sxx)*GradD) / D^2

  HesSxx = GradSxx*(1 - 2*preds)
  #note, S could alway be positive, and maybe the derivative could be set to 1 instead of sigmooth_sign
  #the second derivative to 0 instead of sclamp(1-sclamp)
  HesD = HesSxx

  hessMASE = 2 * sigmoid_clamp(Sxx)*(1 - sigmoid_clamp(Sxx))*GradSxx*GradSxx +
    smooth_sign(Sxx)*( HesSxx)


  hterm1= (hessMASE*D - smooth_sign(Sxx)*GradSxx*GradD) / D^2

  ht2 = smooth_sign(Sxx)*GradSxx*GradD + abs(Sxx)*HesD

  hterm2 = - ( ht2*D^2 - abs(Sxx)*GradD*2*D*GradD) / D^4


  hessSMAPE = hterm1 + hterm2

  grad = 0.5*(ew*gradMASE / avg_mase + eh*gradSMAPE/avg_smape)
  hess = 0.5*(ew*hessMASE / avg_mase + eh*hessSMAPE/avg_smape)

  hess <- pmax(hess, 1e-16)
  return(list(grad = t(grad), hess = t(hess)))
}




fair_obj <- function(preds, dtrain) {
  ff <- attr(dtrain, "ff")
  labels <- attr(dtrain, "xx")
  c <- 2  #the lower the "slower/smoother" the loss is. Cross-Validate.
  x <-  rowSums(preds * ff)-labels
  grad <- ff * (c*x / (abs(x)+c))
  hess <- ff*ff*(c^2 / (abs(x)+c)^2)
  print(sum(abs(x)))
  return(list(grad = grad, hess = hess))
}

ln_cosh_obj <- function(preds, dtrain) {
  ff <- attr(dtrain, "ff")
  labels <- attr(dtrain, "xx")
  ew <- attr(dtrain, "ew")
  norm <- abs(mean(labels))
  x <- (rowSums(preds * ff) - labels) * ew
  grad <- tanh(x)*ff
  hess <- ff*(1-tanh(x)^2)*ff
  print(sum(abs(x)))
  return(list(grad = grad, hess = hess))
}

#this function calculates all the parts that can be
#taken outside of the gradient-hess calculation,
#and added later as a multiplication:
#normalizing via the horizon steps, the denom in the mase, calc, etc..
calc_external_weight <- function(insample, h) {
  frq <- stats::frequency(insample)
  forecastsNaiveSD <- rep(NA,frq)
  for (j in (frq+1):length(insample)){
    forecastsNaiveSD <- c(forecastsNaiveSD, insample[j-frq])
  }
  masep<-mean(abs(insample-forecastsNaiveSD),na.rm = TRUE)
  c( 1 / (masep *h ), 1/h)
}


#' @export
create_combi_info <- function(dataset) {
  #transform the dataset into matrix format
  #with one forecast horizon per row, that will be properly weighted
  #instead of one series per row, to simplify the objective function
  forec_true_w <- lapply(dataset, function (seriesentry) {
    external_weight <- calc_external_weight(seriesentry$x, seriesentry$h)
    t(sapply(1:seriesentry$h, function(i) {
      c(t(seriesentry$ff[, i]), seriesentry$xx[i], external_weight)
    }))
  })
  forec_true_w <- do.call(rbind, forec_true_w)

  data <- lapply(dataset, function (lentry) {
    seriesdata <- t(replicate(lentry$h,
                              as.numeric(lentry$features)))
    colnames(seriesdata) <- names(lentry$features)
    seriesdata
  })
  data <- do.call(rbind, data)

  ff <- forec_true_w[, 1:nrow(dataset[[1]]$ff)]
  xx <- forec_true_w[, nrow(dataset[[1]]$ff) + 1]
  ew <- forec_true_w[, ncol(forec_true_w) - 1]
  eh <- forec_true_w[, ncol(forec_true_w) ]
  list(data=data, ff=ff, xx=xx, ew=ew, eh=eh)
}



train_combination_ensemble <- function(dataset, params=NULL) {

  check_customxgboost_version()

  train_data <- create_combi_info(dataset)

  dtrain <- xgboost::xgb.DMatrix(train_data$data)
  attr(dtrain, "ff") <- train_data$ff
  attr(dtrain, "xx") <- train_data$xx
  attr(dtrain, "ew") <- train_data$ew


  # if (is.null(test_dataset)) {
  #   test_dataset = dataset
  # }
  # test_data <- create_combi_info(test_dataset)
  # dtest <- xgboost::xgb.DMatrix(test_data$data)
  # attr(dtest, "ff") <- test_data$ff
  # attr(dtest, "xx") <- test_data$xx
  # attr(dtest, "ew") <- test_data$ew


  param <- list(max_depth=4, eta=0.019*6, nthread = 2, silent=0,
                objective=combi_softmax_square,
                num_class=nrow(dataset[[1]]$ff),
                subsample=0.9,
                colsample_bytree=0.6)

  bst <- xgboost::xgb.train(param, dtrain, 150)
  bst
}


#' Train Metalearning Models
#' Trains a metalearning model by performing search on the hyperparameters
#'
#'
#' @param train_dataset A list with elements in the metalearning format.
#' E.g. the output of a combination of \code{process_forecast_dataset} and \code{generate_THA_feature_dataset}
#' \describe{
#'   \item{x}{A time series object \code{ts} with the historical data.}
#'   \item{h}{The number of required forecasts.}
#'   \item{xx}{The number of required forecasts.}
#'   \item{THA_features}{The number of required forecasts.}
#'   \item{ff}{The number of required forecasts.}
#'   \item{errors}{The number of required forecasts.}
#' }
#' @param eval_dataset A list in the format of \code{train_dataset}, used for evaluating the error.
#' Can be the same as \code{train_dataset} to show training error
#' @param obj_fun The objective loss function that the metalearning minimizes
#' @param filename Name of the file used for saving the metalearning process
#' @param verbose Boolean indicating whether training progress messages may be printed
#'
#' @return
#' \describe{
#'   \item{model}{The \code{xgboost} model found by the metalearning}
#'   \item{eval_log}{The log of hyperparametes tested with their produced errors}
#' }
#'
#' @export
metatemp_train <- function(train_dataset,
                           eval_dataset,
                           obj_fun = c("select",
                                       "combi:smax_abs",
                                       "combi:smax_sq",
                                       "combi:square"),
                           filename = "meta_results.RData",
                           verbose = FALSE,
                           n.cores=3) {

  check_customxgboost_version()

  type <- match.arg(obj_fun,
                    choices = c("combi:smax_abs",
                                "combi:smax_sq",
                                "combi:square",
                                "select"))
  train_data <- NULL
  dtrain <- NULL
  objective_fun <- NULL

  apply_softmax <- FALSE

  if (type == "select") {
    objective_fun <- error_softmax_obj
    apply_softmax <- TRUE
    train_data <- create_feat_classif_problem(train_dataset)
    dtrain <- xgboost::xgb.DMatrix(train_data$data)
    attr(dtrain, "errors") <- train_data$errors
  } else {
    if (type=="combi:smax_abs") {
      apply_softmax = TRUE
      objective_fun <- combi_softmax_abs
    }
    if (type=="combi:smax_sq") {
      apply_softmax = TRUE
      objective_fun <- combi_softmax_square
    }
    if (type=="combi:square") {
      apply_softmax = FALSE
      objective_fun <- combi_forec_square_obj
    }
    train_data <- create_combi_info(train_dataset)
    dtrain <- xgboost::xgb.DMatrix(train_data$data)
    attr(dtrain, "ff") <- train_data$ff
    attr(dtrain, "xx") <- train_data$xx
    attr(dtrain, "ew") <- train_data$ew
  }


  eval_data <- create_feat_classif_problem(eval_dataset)

  max_depth_grid <- c(4,6,8,10)
  nrounds_grid <- c(20,50,100,200)
  eta_grid <- 2*c(0.001, 0.01, 0.05, 0.1, 0.2)
  colsample_bytree_grid <- c(0.6, 1.0)
  subsample_grid <- c(0.6, 0.9)

  params_grid <- expand.grid(max_depth_grid,
                             nrounds_grid,
                             eta_grid,
                             colsample_bytree_grid,
                             subsample_grid)

  #randomize grid search
  params_grid <- params_grid[sample(nrow(params_grid)),]

  results_grid <- cbind(params_grid, -666, -666)
  colnames(results_grid) <- c("max_depth", "nrounds",
                              "eta", "colsample_bytree",
                               "subsample", "selected_owi_error",
                              "weighted_owi_error")



  best_owi <- 9999999999.9
  meta_results <- NULL
  best_model <- NULL
  for (i in 1:nrow(params_grid)) {
    max_depth <- params_grid[i,1]
    nrounds <- params_grid[i,2]
    eta <- params_grid[i,3]
    colsample_bytree <- params_grid[i,4]
    subsample <- params_grid[i,5]

    if (verbose)  print(paste("Training with: ",
                "max_depth=", max_depth,
                "nrounds=",nrounds,
                "eta=",eta,
                "colsample_bytree", colsample_bytree,
                "subsample", subsample))


    param <- list(max_depth=max_depth,
                  eta=eta, nthread = n.cores,
                  silent=0,
                  objective=objective_fun,
                  num_class=nrow(train_dataset[[1]]$ff),
                  subsample=subsample,
                  colsample_bytree=colsample_bytree)

    bst <- xgboost::xgb.train(param, dtrain, nrounds)

    preds <- stats::predict(bst, eval_data$data, outputmargin = TRUE, reshape=TRUE)
    if (apply_softmax) {
      preds <- t(apply( preds, 1, softmax_transform))
    }
    errors <- summary_performance(preds,
                                     eval_dataset,
                                     print.summary = verbose)
    owi_error <- min(unlist(errors))
    if (is.nan(owi_error)) {
      owi_error <- 999999999.9
    }
    results_grid[i, (ncol(results_grid)-1):ncol(results_grid)] <- unlist(errors)

    if (verbose) print(paste("Iter: ", i, " OWI: ", round(owi_error,3)))

    if (owi_error < best_owi) {
      if (verbose) print("Improved the owi!")
      best_owi <- owi_error
      best_model <- bst
    }
    meta_results <- list(model = best_model,
                         eval_log= results_grid[1:i,])
    save(meta_results, file=filename)
  }
  meta_results
}
