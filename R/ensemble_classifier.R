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


#CANDIDATE TO REMOVE, ZERO CLAMPING....
#' @export
ensemble_forecast <- function(predictions, dataset, clamp_zero=TRUE) {
  for (i in 1:length(dataset)) {
    weighted_ff <- as.vector(t(predictions[i,]) %*% dataset[[i]]$ff)
    if (clamp_zero) {
      weighted_ff[weighted_ff < 0] <- 0
    }
    dataset[[i]]$y_hat <- weighted_ff
  }
  dataset
}



#' @describeIn train_selection_ensemble Analysis of the predictions, the weighted error and the selection error, along with extra information
#' @param predictions A NXM matrix with N the number of observations(time series)
#'                    and M the number of methods. Each row contains the weights assigned to the methods for the series
#' @param dataset The list with the meta information, forecasts of each method...MUST contain the precalculated errors with process_errors()!
#' @param print.summary Boolean indicating wheter to print the information
#' @param use.precalc.naive2 Boolean indicating wheter the naive2 errors are already contained in \code{dataset} for skip its calculation
#'
#' @export
summary_performance <- function(predictions, dataset, print.summary = TRUE, use.precalc.naive2=FALSE) {
  stopifnot("xx" %in% names(dataset[[1]]))

  #requires precalculated average errors of the dataset
    if (is.null(attr(dataset, "avg_naive2_errors") )) {
      stop("summary_performance requires precalculated avg_naive2_errors, please run process_errors on dataset")
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

    dataset <- process_errors(dataset)



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
  })
}

#' @export
forecast_meta_M4 <- function(meta_model, x, h) {

  ele <- list(list(x=x, h=h))
  ele <- process_THA_features(ele)

  ff <- process_forecasts(ele[[1]], forec_methods())
  ff <- t(sapply(ff, function (lentry) lentry$forecast))
  rownames(ff) <- unlist(forec_methods())

  preds <- predict_selection_ensemble(meta_model, as.matrix(ele[[1]]$features))

  y_hat <- preds %*% ff

  y_hat <- as.vector(y_hat)

  #for the interval forecast we use only thetaf, naive and snaive
  thetamod <- forecast::thetaf(x, h)
  radius_theta <- thetamod$upper[,2] - thetamod$mean

  naivemod <-forecast::naive(x, h)
  radius_naive <- naivemod$upper[,2] - naivemod$mean

  snaivemod <- forecast::snaive(x, h)
  radius_snaive <- snaivemod$upper[,2] - snaivemod$mean

  ff_radius <- rbind(radius_theta, radius_naive, radius_snaive)

  radius <- rep(0, h)
  if (h > 48) {
    stop("Maximum forecasting horizon is 48")
  }
  for (i in 1:h) {
    radius[i] <- sum(M4_interval_weights[,i] * ff_radius[,i])
  }

  upper <- as.vector(y_hat + radius)
  lower <- as.vector(y_hat - radius)

  y_hat[y_hat < 0] <- 0
  upper[upper < 0] <- 0
  lower[lower < 0] <- 0

  list(mean=y_hat, upper=upper, lower=lower)
}

#' @export
get_M4_interval_weights <- function() {
  M4_interval_weights
}

#weights used for calu
M4_interval_weights <- structure(c(1.20434805479858, 0.0349832013212199, -0.0135553803216437,
            1.6348667017876, -0.0291391104169728, -0.137106578797595, 1.89343593467985,
            -0.0556309211578089, -0.193624560527705, 1.97820794828226, -0.0710718154535777,
            -0.123100424872484, 2.20852069562932, -0.0710600158690502, -0.211316149359757,
            2.00533194439932, -0.0735768842959382, 0.00639025637846198, 1.54942520065088,
            -0.057520179590224, 0.0192717604264258, 1.66779655654594, -0.0577404254880085,
            -0.0301536402316914, 1.57793111309566, -0.0571366696960703, 0.00074454324835706,
            1.58014785718421, -0.0546758191445934, 0.00023514471533616, 1.51051833028159,
            -0.0502489429193765, 0.0306686761759747, 1.49667215498576, -0.0504760468509232,
            0.0633292476660604, 1.46230478866238, -0.0532688561557852, 0.0966347659478812,
            1.61065833858158, -0.0586474808085979, -0.010800707477267, 1.50576342267903,
            -0.0520210844283552, 0.0144436883743055, 1.54921239534018, -0.053533338280372,
            -0.0108667128746493, 1.53417971494702, -0.0579026110677012, 0.000297052007093576,
            1.55732972883077, -0.0553048336745351, 0.0367638357882623, 0.116813998810538,
            0.323446912667969, 0.630743390147047, 0.73284831139316, -0.0328228748038042,
            0.614008235063802, 1.76483588802515, -0.0543503687837068, -0.0194824287225068,
            0.40893801004056, -0.0372774403169106, 1.03983187771261, 0.874914242855032,
            0.0520531705411014, 0.265317677176088, 1.43631066965087, -0.0197534066247038,
            -0.173619599584651, 0.898313936819893, -0.0705817867062166, 0.325840179586548,
            0.625354817708426, 0.0730184383944889, 0.327004041690249, 0.311911360827102,
            -0.045054781310374, 1.04704513182606, 0.141901550600215, -0.0382311974845755,
            1.34096763681398, 0.491822729685341, -0.0374825159401713, 0.915188301573958,
            0.742763540929105, -0.0735719158347089, 0.906734816018297, 0.997546107480632,
            -0.0835548323138982, 0.544652116459324, 0.687991299960411, -0.0463052425732432,
            0.645709731320379, 0.564295324221999, -0.0634237302002972, 1.09921867894254,
            -0.0931400534770518, -0.0191455415920918, 1.42725756848195, 0.425542670313752,
            -0.0525723014099128, 1.11315125793807, 2.37441008580272, -0.182316365736858,
            0.208272797333258, 1.37681109893514, -0.109037211688484, 0.322339972304184,
            1.02158407919349, -0.0702467319013755, 0.547900910728034, 0.939473970980379,
            -0.0831346087836634, 0.710478207012984, 0.088098458991328, -0.0279187923322283,
            1.24491828863518, 1.03485500909627, -0.0845414892639553, 0.418196261825132,
            0.432858831338303, -0.0486860688171931, 0.960107743238668, 0.0505708749888397,
            -0.0195777073302953, 0.957162564794512, 0.520298978200399, 0.0250325108758569,
            0.50637198342853, 0.502726050210653, -0.0498970101052844, 0.7843393684093,
            1.19955393367803, -0.0995989563120128, 0.592116095974745, 0.0150399012793929,
            -0.0310974825919329, 1.8790134248291, 1.42083348006252, -0.30878943417154,
            1.97721144553016), .Dim = c(3L, 48L))

