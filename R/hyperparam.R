

#' @export
hyperparameter_search <- function(dataset, filename="meta_hyper.RData", n_iter=1000, n.cores=1) {

  if (length(dataset) < 10) {
    stop("Not enough data to do the crossvalidation!")
  }

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



  bayes_xgb <- function(max_depth, eta, subsample, colsample_bytree, nrounds) {

    n_class <- nrow(train_ds[[1]][[1]]$ff)
    bay_results <- NULL

    param <- list(max_depth=max_depth, eta=eta, nthread = N_THREAD, silent=0,
                  objective=error_softmax_obj,
                  num_class=n_class,
                  subsample=subsample,
                  colsample_bytree=colsample_bytree)
    final_error = NULL
    final_preds = NULL
    for (i in 1:1) {

      dtrain <- xgboost::xgb.DMatrix(train_feat[[i]]$data)
      attr(dtrain, "errors") <- train_feat[[i]]$errors

      bst <- xgboost::xgb.train(param, dtrain, nrounds)
      preds <- M4metalearning::predict_selection_ensemble(bst, test_feat[[i]]$data)

      attr(test_ds[[i]], "avg_naive2_errors") <- attr(dataset, "avg_naive2_errors")
      er <- M4metalearning::summary_performance(preds,
                                                test_ds[[i]],
                                                print.summary = FALSE, use.precalc.naive2 = TRUE)

      final_error <- c(final_error, er$weighted_error)
      final_preds <- rbind(final_preds, preds)
    }

    try({load(filename)})
    bay_results <- rbind(bay_results,
                         c(max_depth, eta, subsample, colsample_bytree, nrounds, mean(final_error))
    )
    try({colnames(bay_results) <- c("max_depth", "eta", "subsample", "colsample_bytree", "nrounds", "combi_OWA")})
    bay_results <- as.data.frame(bay_results)
    save(bay_results, file=filename)
    list(Score=-mean(final_error), Pred=final_preds)
  }

  prefound_grid <- list(max_depth=10L,
                        eta=0.4,
                        subsample=0.9,
                        colsample_bytree=0.6,
                        nrounds=200)

  k=2
  bay_res <- rBayesianOptimization::BayesianOptimization(bayes_xgb, bounds=list(max_depth=c(2L,30L),
                                                         eta=c(0.001, 1.0),
                                                         subsample=c(0.5,1.0),
                                                         colsample_bytree=c(0.5,1.0),
                                                         nrounds=c(1L,350L)),
                                  init_grid_dt = prefound_grid,
                                  init_points= 5,
                                  kappa = 2.576,
                                  n_iter=n_iter,
                                  kernel=list(type = "matern", nu=(2*k+1)/2))
  bay_res
}
