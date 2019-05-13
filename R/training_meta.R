
#' @export
train_metalearning <- function(dataset, forec_methods,
                               objective = "averaging",
                               chunk_size=NULL,
                               save_filename="tmp_train_meta.RData",
                               resume_filename=NULL) {

  TMP_FOREC_METHODS <- M4metalearning::forec_methods()


  #RECOMMENDED FOR PARALLELISM, future plan with gc=TRUE



  # a one function without parameters

  #all the function call parameters plus the state variable are saved on checkpoints
  save_training <- function() {
    save(dataset, state, objective, save_filename, chunk_size,
         file = save_filename)
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


  state <- "STARTING"




  #the intermediate files to store chunkified processes and hyperparameter search
  proc_save_filename <- append_suffix(save_filename, suffix="_proc.")
  proc_resume_filename <-  append_suffix(resume_filename, suffix="_proc.")


  #when we are the resuming, check the state and go to the proper state, good old spaghetti code
  if (!is.null(resume_filename)) {
    load(resume_filename)
  }



  #state machine, to guide which processing is further required
  do_holdout <- state == "STARTING"
  do_forecast <- do_holdout || state == "HOLDOUT_DONE"
  do_features <- do_forecast || state == "FORECAST_DONE"
  do_errors <- do_features || state == "FEATURES_DONE"
  do_hyper <- do_errors || state == "CALCERRORS_DONE"

  ##get the number of workers from the future plan, used in auto chunk_size calculation and to launch xgboost
  strat <- future::plan()
  num_workers <- as.list(args(strat))$workers
  library(future) #ugly namespace thing
  num_workers <- eval(num_workers)
  if (is.null(num_workers)) {
    num_workers = 1
  }

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
  }


  if (do_forecast) {
    if (state != "HOLDOUT_DONE") {
      message("Holdout calculated, saving...")
      state <- "HOLDOUT_DONE"
      proc_resume_filename <- NULL #we are doing forecast for the first time, we cannot resume
      save_training()
    }
    message("Processing forecasts...")
    dataset <- M4metalearning::process_forecasts(dataset, TMP_FOREC_METHODS,
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

  }
  state <- "FINISHED"
  save_training()

  best_params <- bayes_results[which.min(bayes_results[, ncol(bayes_results)]), ]

  meta_model <- .train_from_bayes_res(dataset, best_params, n.cores = num_workers)


  gc()

  list(dataset=dataset, meta_model=meta_model, bayes_results=bayes_results)

}
