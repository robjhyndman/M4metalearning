##get the number of workers from the future plan, used in auto chunk_size calculation and to launch xgboost
autodetect_num_workers <- function() {
  strat <- future::plan()
  num_workers <- as.list(args(strat))$workers
  library(future) #ugly namespace thing
  num_workers <- eval(num_workers)
  if (is.null(num_workers)) {
    num_workers = 1
  }
  num_workers
}


calculate_chunk_size <- function(dataset, num_workers) {
  N <- length(dataset)
  data_size <- object.size(dataset) / (1024*1024)


  #maximum of each chunk is ... in MegaBytes
  MAX_MEMORY_PERCHUNK = 45
  #try to get to the largest chunk size which is multiple of the workers and under MAX_MEMORY_PERCHUNK
  max_chunk_size <- MAX_MEMORY_PERCHUNK / (data_size / N )
  chunk_size <- num_workers*floor(max_chunk_size / num_workers)
  if (chunk_size > N/2) {
    chunk_size <- floor(N/2)
  }
  message(paste("Auto calculated chunk size: ", chunk_size))
  chunk_size
}

#quick function to append a suffix to a filename so that it goes before the "." extension of the filename
append_suffix <- function(fname_string, suffix) {
  if (is.null(fname_string)) {
    return(fname_string)
  }
  tmp_split <- unlist(strsplit(fname_string, split="\\."))
  paste(c(tmp_split[1:(length(tmp_split) -1)], suffix, tmp_split[length(tmp_split)])
        , collapse="")
}


#' @export
train_metalearning <- function(dataset, forec_methods = M4_forec_methods(),
                               objective = "averaging",
                               chunk_size=NULL,
                               save_filename="tmp_train_meta.RData",
                               resume_filename=NULL) {

  bayes_results <- NULL
  #all the function call parameters plus the state variable are saved on checkpoints
  save_training <- function() {
    save(dataset, state, objective, save_filename, chunk_size, forec_methods, bayes_results,
         file = save_filename)
  }



  state <- "STARTING"




  #the intermediate files to store chunkified processes and hyperparameter search
  proc_save_filename <- append_suffix(save_filename, suffix="_proc.")
  proc_resume_filename <-  append_suffix(resume_filename, suffix="_proc.")


  #when we are the resuming, check the state and go to the proper state, good old spaghetti code
  if (!is.null(resume_filename)) {
    load(resume_filename)
  }

  save_training()

  #state machine, to guide which processing is further required
  do_holdout <- state == "STARTING"
  do_forecast <- do_holdout || state == "HOLDOUT_DONE"
  do_features <- do_forecast || state == "FORECAST_DONE"
  do_errors <- do_features || state == "FEATURES_DONE"
  do_hyper <- do_errors || state == "CALCERRORS_DONE"

  ##get the number of workers from the future plan, used in auto chunk_size calculation and to launch xgboost

  num_workers <- autodetect_num_workers()

  if (do_holdout) {
    message("Calculating holdout data...")
    ##check if we have true future values in $xx, if we do, skip the holdout, if we dont, use the horizon $h for temporal holdout
    ##if we do not have horizon, set $h to a period of frequency, if frequency is 1 then set $h to 6
    dataset <- lapply(dataset, function (ll) {
      #we can add here some sanity checks

      #check if we have prediction horizon, if not
      if (is.null(ll$xx)) {
        if (is.null(ll$h)) {
          ll$h <- frequency(ll$x)
          if (ll$h == 1) {
            ll$h <- 4
          }
        }
      }
      ll})

    dataset <- M4metalearning::temp_holdout(dataset)
    gc()

  }




  #process the forecasts
  #check if we are resuming calculations

  #calculate chunks and everything
  if (is.null(chunk_size)) {
    chunk_size <- calculate_chunk_size(dataset,   num_workers)
  }


  if (do_forecast) {
    if (state != "HOLDOUT_DONE") {
      message("Holdout calculated, saving...")
      state <- "HOLDOUT_DONE"
      proc_resume_filename <- NULL #we are doing forecast for the first time, we cannot resume
      save_training()
    }
    message("Processing forecasts...")
    dataset <- M4metalearning::process_forecasts(dataset, forec_methods,
                                                 chunk_size = chunk_size, do_shuffle = TRUE,
                                                 save_checkpoint_filename = proc_save_filename,
                                                 load_checkpoint_filename = proc_resume_filename)
    gc()
  }

  if (do_features) {
    if (state != "FORECAST_DONE") {
      #save the state, set it to forecast_done
      message("Forecasts calculated, saving...")
      state <- "FORECAST_DONE"
      proc_resume_filename <- NULL  #we are doing features for the first time, we cannot resume
      save_training()
    }
    #go for the features
    message("Calculating features...")
    dataset <- M4metalearning::process_THA_features(dataset, chunk_size,
                                                    do_shuffle = TRUE,
                                                    save_checkpoint_filename = proc_save_filename,
                                                    load_checkpoint_filename = proc_resume_filename)
    gc()
    #save the state, set it to features_done
    state <- "FEATURES_DONE"
    message("Features calculated, saving...")
    save_training()
  }

  if (do_errors) {
    #calculate the errors
    dataset <- M4metalearning::process_errors(dataset)
    gc()
  }

  #now do the hyperparameter search


  ### from this point on, we can already use the model!!!!




  if (do_hyper) {
    if (state != "CALCERRORS_DONE") {
      state <- "CALCERRORS_DONE"
      proc_resume_filename <- NULL
      message("Errors calculated, saving...")
      save_training()
    }
    bayes_results <- M4metalearning::hyperparameter_search(dataset, objective = objective,
                                                           n_iter=5, n.cores=num_workers,
                                                           rand_points = 4,
                                                           save_filename=proc_save_filename,
                                                           resume_filename=proc_resume_filename)
    state <- "HYPER_DONE"
    save_training()
  }


  best_params <- bayes_results[which.min(bayes_results[, ncol(bayes_results)]), ]

  meta_model <- .train_from_bayes_res(dataset, best_params, n.cores = num_workers)


  gc()


  list(dataset=dataset, meta_model=meta_model, forec_methods=forec_methods,
       objective=objective,
       bayes_results=bayes_results)

}


#forecast meta

#pass the object for training, and a dataset
#uses similar structure as the R predict

#returns a dataset with added component y hat and the individual forecasts

#last step is to add the summary performance, we save it for the output also
#we can do the comparison through crossvalidation if sample and test are a holdout version

#the



#' @export
forecast_meta <- function(model, new.dataset,
                          chunk_size=NULL,
                          save_filename="tmp_forec_meta.RData",
                          resume_filename=NULL) {

  num_workers <- autodetect_num_workers()
  if (is.null(chunk_size)) {
    chunk_size <- calculate_chunk_size(new.dataset, num_workers)
  }

  weights <- NULL
  state <- "STARTING"


  save_forecasting <- function() {
    save(new.dataset, state, save_filename, chunk_size, weights,
         file = save_filename)
  }


  #the intermediate files to store chunkified processes and hyperparameter search
  proc_save_filename <- append_suffix(save_filename, suffix="_proc.")
  proc_resume_filename <-  append_suffix(resume_filename, suffix="_proc.")



  #when we are the resuming, check the state and go to the proper state, good old spaghetti code
  if (!is.null(resume_filename)) {
    load(resume_filename)
  }

  save_forecasting()

  do_features <- state == "STARTING"
  do_weights <- do_features || state == "FEATURES_DONE"
  do_forecast <- do_weights || state == "WEIGHTS_DONE"



  #calculate the features
  if (do_features) {
    message("Calculating features...")
    new.dataset <- M4metalearning::process_THA_features(new.dataset, chunk_size,
                                                        do_shuffle = TRUE,
                                                        save_checkpoint_filename = proc_save_filename,
                                                        load_checkpoint_filename = proc_resume_filename)
    gc()

  }

  #calculate the weights, attach them to our dataset
  if (do_weights) {
    if (state != "FEATURES_DONE") {
      state <- "FEATURES_DONE"
      proc_resume_filename <- NULL
      message("Features calculated, saving...")
      save_forecasting()
    }
    new.eval_data <- M4metalearning::create_feat_classif_problem(new.dataset)

    weights <- M4metalearning::predict_selection_ensemble(model$meta_model, new.eval_data$data)
    state <- "WEIGHTS_DONE"
    message("Weights calculated, saving...")
    save_forecasting()
  }

  #calculate the forecasts from the individual models
  if (do_forecast) {
    message("Processing forecasts...")
    new.dataset <- M4metalearning::process_forecasts(new.dataset, model$forec_methods,
                                                     chunk_size = chunk_size, do_shuffle = TRUE,
                                                     save_checkpoint_filename = proc_save_filename,
                                                     load_checkpoint_filename = proc_resume_filename)
    gc()

    state <- "FORECAST_DONE"
    save_forecasting()
  }


  for (i in 1:length(new.dataset)) {
    new.dataset[[i]]$meta_ff_average <- weights[i,] %*% new.dataset[[i]]$ff
    new.dataset[[i]]$meta_ff_selection <- new.dataset[[i]]$ff[which.max(weights[i,] ),]
    new.dataset[[i]]$meta_weights <- weights[i,]
  }
  #attach them to the output dataset

  #performance summary if we have the xx in the input dataset
  errors <- NULL
  if ( !is.null(new.dataset[[1]]$xx)) {
    new.dataset <- M4metalearning::process_errors(new.dataset)
    errors <- M4metalearning::summary_performance(weights, new.dataset)
  }

  #attach the predictions to the output dataset
  retval <- list(dataset=new.dataset)
  if (!is.null(errors)) {
    retval$errors <- errors
  }

  retval
}

