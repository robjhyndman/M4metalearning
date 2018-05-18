dirty_calc_owi <- function(insample, outsample, forecast) {
  snaive_ff <- snaive_forec(insample, length(forecast))
  snaive_errors <- calculate_errors(insample, outsample, snaive_ff)
  calculate_owi(insample, outsample, snaive_errors, forecast)
}


#' Calculate the theoretical error for subsets of the pool of methods
#' Calculates the error for each possible subset of method in the processed dataset
#'
#' @param dataset A list of time series, with each element having at least
#'                the entry \code{errors}, a vector with the errors produced by
#'                the methods in the pool
#' @return A list of list, the first dimension is the size of of subset
#'                         and the second the results of each subset,
#'                        having two entries: \code{iset} the indices of the subset in the
#'                        \code{errors} vector and \code{error}, the minimum error produced
#'                        by that subset.
#' @export
calc_subset_methods <- function(dataset) {
  errors <- t(sapply(dataset, function (lentry) lentry$errors))
  cinfo <- create_combi_info(dataset)

  results_list <- replicate(ncol(errors),list())
  for (N in 1:(ncol(errors))) {
    comb <- utils::combn(ncol(errors), N)
    for (K in 1:ncol(comb)) {
      iset <- comb[,K]
      if (N==1) {
        iset <- c(iset,iset)
      }
      oracle_error_set <- mean(apply(errors[, iset], 1, min))
      naive_combi_error_set <- mean(abs(rowMeans(cinfo$ff[,iset]) - cinfo$xx))
      if (N==1) {
        iset <- iset[1]
      }
      results_list[[N]] <- append(results_list[[N]],
                                  list(list(iset = iset,
                                            selec_error = oracle_error_set,
                                            naivecomb_error = naive_combi_error_set)))
    }
    #save(results_list, file="pool_methods.RData")
  }
  results_list
}
