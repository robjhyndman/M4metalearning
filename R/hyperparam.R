
#https://towardsdatascience.com/doing-xgboost-hyper-parameter-tuning-the-smart-way-part-1-of-2-f6d255a45dde
#https://www.analyticsvidhya.com/blog/2016/03/complete-guide-parameter-tuning-xgboost-with-codes-python/
#https://towardsdatascience.com/fine-tuning-xgboost-in-python-like-a-boss-b4543ed8b1e
#https://www.kaggle.com/general/17120
#https://www.analyticsvidhya.com/blog/2016/01/xgboost-algorithm-easy-steps/
#https://www.slideshare.net/ShangxuanZhang/kaggle-winning-solution-xgboost-algorithm-let-us-learn-from-its-author
#https://www.slideshare.net/odsc/owen-zhangopen-sourcetoolsanddscompetitions1

#' @export
hyperparameter_search <- function(dataset, objective= c("averaging", "selection"),
                                  rand_points=100,
                                  n_iter=1000, n.cores=1,
                                  save_filename="bayes_hyper_search.RData",
                                  resume_filename=NULL) {

  if (is.null(save_filename)) {
    stop("hyperparameter_search requires a file to store temp results, please set the save_filename parameter")
  }

  #weird check to remove a save file already exisiting, when we are not resuming, to avoid continoing from it
  if (!is.null(save_filename)) {
    if (is.null(resume_filename) || (save_filename!=resume_filename)) {
      if (file.exists(save_filename)) {
        message("Rewriting saving file...")
        file.remove(save_filename)
      }
    }
  }


  if (length(dataset) < 10) {
    stop("Not enough data to do the crossvalidation!")
  }

  if (is.null(attr(dataset, "avg_naive2_errors"))) {
    stop("Need to calculate the average naive2 OWA errors with process_errors")
  }

  type_objective <- match.arg(objective)

  N_THREAD = n.cores
  whole_dataset <- dataset
  #prepare the folds
  folds <- rBayesianOptimization::KFold(1:length(whole_dataset), nfolds=5, seed=31-05-2018)

  train_ds <- NULL
  test_ds <- NULL
  train_feat <- NULL
  test_feat <- NULL

  for (i in 1:length(folds)) {
    train_ds[[i]] <- whole_dataset[ -folds[[i]] ]
    train_feat[[i]] <- create_feat_classif_problem(train_ds[[i]])

    test_ds[[i]] <- whole_dataset[ folds[[i]] ]
    test_feat[[i]] <- create_feat_classif_problem(test_ds[[i]])
  }

  bay_results <- NULL


  bayes_xgb <- function(max_depth, eta, gamma, min_child_weight,
                        subsample, colsample_bytree, nrounds) {


    param_bay <- list(max_depth=max_depth, eta=eta,
                  gamma=gamma,
                  min_child_weight=min_child_weight,
                  subsample=subsample,
                  colsample_bytree=colsample_bytree,
                  nrounds=nrounds)

    final_error = NULL
    final_preds = NULL
    for (i in 1:1) {

      bst <- .train_data_from_bayes_res(train_feat[[i]], param_bay, N_THREAD)
      preds <- M4metalearning::predict_selection_ensemble(bst, test_feat[[i]]$data)

      attr(test_ds[[i]], "avg_naive2_errors") <- attr(dataset, "avg_naive2_errors")
      er <- M4metalearning::summary_performance(preds,
                                                test_ds[[i]],
                                                print.summary = FALSE, use.precalc.naive2 = TRUE)
      #maybe improve this a bit to avoid calculating both errors always
      er <- switch(type_objective,
        selection = er$selected_error,
        averaging = er$weighted_error)

      final_error <- c(final_error, er)
      final_preds <- rbind(final_preds, preds)
    }

    bay_results <- rbind(NULL, c(max_depth, eta, gamma, min_child_weight,
                     subsample, colsample_bytree, nrounds, mean(final_error)))

    try({colnames(bay_results) <- c("max_depth", "eta", "gamma", "min_child_weight",
                                    "subsample", "colsample_bytree", "nrounds", "Value")})
    bay_results <- data.frame(bay_results)

    if (!is.null(save_filename)) {
      oldres <- NULL
      try(oldres<-readRDS(save_filename))
      bay_results <- rbind(oldres, bay_results)
      saveRDS(bay_results, file=save_filename)
    }
    list(Score=-mean(final_error), Pred=0)
  }


  precalc_grid <- NULL
  if (!is.null(resume_filename)) {
    message("Resuming hyperparameter search")
    precalc_grid <- readRDS(resume_filename)
    bay_results <- precalc_grid

    n_iter <- n_iter +  min(rand_points - nrow(bay_results), 0)#when resuming calc how many iterations are left
    rand_points <- max(rand_points - nrow(bay_results), 0) #we are resuming, we assume the first points come from random
  }

  if (rand_points + n_iter > 0) {
  k=2
  rBayesianOptimization::BayesianOptimization(bayes_xgb, bounds=list(max_depth=c(2L,50L),
                                                         eta=c(0.001, 1.0),
                                                         gamma=c(0.00001, 2.0),
                                                         min_child_weight=c(0.00001, 5.0),
                                                         subsample=c(0.5,1.0),
                                                         colsample_bytree=c(0.5,1.0),
                                                         nrounds=c(1L,500L)),
                                  init_grid_dt = precalc_grid,
                                  init_points= rand_points,
                                  kappa = 2.576,
                                  n_iter=n_iter,
                                  kernel=list(type = "matern", nu=(2*k+1)/2))
  }

  bay_results <- readRDS(save_filename)
  bay_results
}

#from a bay_results row, train the model

.train_data_from_bayes_res <- function (data, bayes_res, n.cores) {
  param <- as.list(bayes_res)
  param$nthread <- n.cores

  nrounds <- param$nrounds
  param$nrounds <- NULL
  param$Value <- NULL


  train_selection_ensemble(data$data, data$errors, param, nrounds)
}

#' @export
.train_from_bayes_res <- function (dataset, bayes_res, n.cores) {
  data <- create_feat_classif_problem(dataset)
  .train_data_from_bayes_res(data, bayes_res, n.cores)
}
