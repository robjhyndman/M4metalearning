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
train_metalearning <- function(train_dataset, forec_methods = M4_forec_methods(),
                               objective = "averaging",
                               chunk_size=NULL,
                               save_foldername=NULL) {

  num_workers <- autodetect_num_workers()
  if (is.null(chunk_size)) {
    chunk_size <- calculate_chunk_size(train_dataset, num_workers)
  }

  train_proc <- function(lentry, methods_list) {
    lentry <- temporal_holdout(lentry)
    lentry <- calc_forecasts(lentry, methods_list)
    lentry <- calc_features(lentry)
    lentry <- calc_mase_smape_errors(lentry)
  }

  train_dataset <- chunk_xapply(train_dataset, chunk_size, save_foldername, "train_call",
                                future.apply::future_lapply,
                                train_proc, forec_methods)

  train_dataset <- process_owa_errors(train_dataset)


  ##hyper search uses a file internally
  bayes_resume_filename <- NULL
  bayes_save_filename <- "meta_bayes_hypersearch.rds"
  if (!is.null(save_foldername)) {
    bayes_save_filename <- paste(save_foldername, "/", bayes_save_filename, collapse="", sep="")
    if (file.exists(bayes_save_filename)) {
      bayes_resume_filename <- bayes_save_filename
    }
  }

  bayes_results <- M4metalearning::hyperparameter_search(train_dataset, objective = objective,
                                                         n_iter=5, n.cores=num_workers,
                                                         rand_points = 4,
                                                         save_filename=bayes_save_filename,
                                                         resume_filename=bayes_resume_filename)

  best_params <- bayes_results[which.min(bayes_results[, ncol(bayes_results)]), ]

  meta_model <- .train_from_bayes_res(train_dataset, best_params, n.cores = num_workers)

  list(train_dataset=train_dataset, meta_model=meta_model, forec_methods=forec_methods,
       objective=objective,
       bayes_results=bayes_results)
}

#' @export
#' @import memoise
forecast_metalearning <- function(model, new_dataset,
                          chunk_size=NULL,
                          save_foldername=NULL) {

  num_workers <- autodetect_num_workers()
  if (is.null(chunk_size)) {
    chunk_size <- calculate_chunk_size(new_dataset, num_workers)
  }


  #all the processing steps for doing the forecasting
  forec_steps <- function (lentry, model) {
    lentry <- calc_forecasts(lentry, model$forec_methods)
    lentry <- calc_features(lentry)
    lentry <- predict_weights_meta(lentry, model$meta_model)
    lentry <- ensemble_meta(lentry)
    lentry <- calc_mase_smape_errors(lentry)
    lentry
  }

  new_dataset <- chunk_xapply(new_dataset, chunk_size, save_foldername, "forec_call",
               future.apply::future_lapply,
               forec_steps, model)

  owa_errors <- NULL
  if (!is.null(new_dataset[[1]]$xx)) {
    new_dataset <- process_owa_errors(new_dataset)
    owa_errors <- summary_meta(new_dataset)
  }
  list(dataset=new_dataset, owa_errors=owa_errors)
}

#' @import memoise
chunk_xapply <- function( .chunk_dataset, chunk_size, save_foldername, .idcall, .apply_FUN, ... ) {

  if (!is.null(save_foldername)) {
    message(paste("using cache in:", save_foldername, "for saving/resuming computations"))
  }
  temp_dataset <- NULL
  data_length <- length(.chunk_dataset)

  chunk_index <- seq(1, data_length, chunk_size)
  chunk_index <- c(chunk_index, data_length+1) #add a final DUMMY chunk
  start_chunk <- 1


  start_time = proc.time()
  for (i in start_chunk:(length(chunk_index)-1)) {
    start_ind = chunk_index[i]
    end_ind = chunk_index[i+1]-1

    if (!is.null(save_foldername)) {
      whether_memoize <- memoize
    } else {
      whether_memoize <- function(x, cache) x
    }
    memofun <- whether_memoize( function(call_id, chunk_id, myMeMoFun) {
      myMeMoFun(.chunk_dataset[start_ind:end_ind],...)
      },
                        cache=cache_filesystem(save_foldername) )
    temp_dataset <- c(temp_dataset, memofun(.idcall, start_ind, .apply_FUN))

    #remaining time calculations
    endchunk_time <- proc.time()

    message(paste("From ", start_ind, " to", end_ind,
                  ", ", round(100*(end_ind) / length(.chunk_dataset),2),
                  "% of the dataset processed, remaining time: ",
                  round( (endchunk_time - start_time)[3]* (length(.chunk_dataset) / (end_ind) -1), 2 ),
                  "seconds") )
  }
  temp_dataset
}





#############old forec

if (0) {


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

