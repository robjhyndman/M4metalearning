Reproducibility: Combination of Forecast Methods by Feature-based Learning
================
Pablo Montero-Manso
2018-06-10

This page explains how to reproduce the results for our submission. The first part is the methodological description of our methods. [The methodology can be seen here](M4_methodology.html)

Authorship
==========

-   Pablo Montero-Manso
-   Thiyanga Talagala
-   Rob J Hyndman
-   George Athanasopoulos

Overview
========

The technical part of the reproducibility is divided in two sections. First, the tools for producing the forecast of our method and second, a way of replicating the training of our model. It is divided in two parts because replicating the training is a very computationaly expensive process, so a pretrained model is provided. Also for very long time series, the actual calculation of there forecast can also computationally expensive, because 9 different forecasting methods are calculated independently, and some of these methods do exhaustive search on its parameters, which is time consuming.

The model
=========

Two R packages are provided as part of our submission, the `M4metalearning` package contains the functions for training new models and making the predictions.

The `M4metaresults` package is just a 'data' package containing the trained model used in ou submission, as well as the training set, explained in our [methodology page](M4_methodology.md) as well as the 'test' or 'final' dataset, which is the original M4 dataset but with the forecasts of the individual methods and the features of the series already extracted.

We first show how to reproduce the exact results of our submission and also how to easily use our trained model and provided tools to produce new forecast for any time series.

Simple forecasting
------------------

The function `forecas_M4` is included in the M4metalearning package to simplify forecasting. `forecast_M4` uses the pretrained model of the `M4metaresults` package. This takes as input a time series and its required forecasting horizon and outputs the mean and interval forecasts. It may be computationally expensive, especially for long time series (length&gt;1000).

\*\*NOTE: The results may slightly vary since one of the individual forecasting methods, `nnetar` uses random initialization. Nevertheless, the results should vary less than 0.1%.\*

The following example shows how to produce forecasts for a example 'new' time series and serveral of the ones in the M4 dataset.

``` r
#install the package if not already installed
devtools::install_github("robjhyndman/M4metalearning")

devtools::install_github("carlanetto/M4comp2018")

#the M4metaresults package must be installed from sources due to its large size
#note that this package is about 1GB in size
install.packages("https://github.com/pmontman/M4metaresults/releases/download/v0.0.0.9000/M4metaresults_0.0.0.9000.tar.gz", 
                 repos = NULL, type="source")
```

``` r
library(M4metalearning)
library(M4metaresults)
#we will showcase how to forecast a syntethic time series
#first, we generate the full series
set.seed(10-06-2018)
truex = (rnorm(60)) + seq(60)/10

#we subtract the last 10 observations to use it as 'true future' values
#and keep the rest as the input series in our method
h = 10
x <- head(truex, -h)
x <- ts(x, frequency = 1)

#forecasting with our method using our pretrained model in one line of code
#just the input series and the desired forecasting horizon
forec_result <- forecast_meta_M4(model_M4, x, h=h)
#> Loading required package: tsfeatures
#show the output, the mean, upper and lower interval forecasts
print(forec_result, digits=2)
#> $mean
#>  [1] 5.3 5.3 5.4 5.5 5.5 5.6 5.6 5.7 5.8 5.8
#> 
#> $upper
#>  [1]  8.1  8.6  8.9  9.5  9.7 10.7  9.8  9.9 10.0 10.2
#> 
#> $lower
#>  [1] 2.37 2.01 1.84 1.39 1.37 0.48 1.51 1.53 1.52 1.47
```

``` r
plot(truex, ylim=range(c(truex,unlist(forec_result))), type="l")
lines(c(x,forec_result$mean), col="blue")
lines(c(x,forec_result$upper), col="blue")
lines(c(x,forec_result$lower), col="blue")
lines(truex)
```

![example forecast](/docs/example_forecast-1.png)

Any time series of the M4 dataset may be easily forecasted using our `forecast_meta_M4` and the `M4comp2018` package containing the original M4 series.

``` r
forec_M4_2018 <- forecast_meta_M4(model_M4,
                                  M4comp2018::M4[[2018]]$x,
                                  M4comp2018::M4[[2018]]$h)
print(forec_M4_2018, digits=2)
#> $mean
#>  [1] 12287 12279 12273 12267 12261 12255 12249 12244 12239 12234 12229
#> [12] 12225 12221 12217
#> 
#> $upper
#>  [1] 12623 12838 13033 13216 13403 13517 13312 13429 13450 13513 13541
#> [12] 13612 13662 13743
#> 
#> $lower
#>  [1] 11950 11721 11513 11317 11118 10993 11187 11059 11028 10955 10918
#> [12] 10837 10779 10690
```

Exact reproduction of the submitted forecasts
---------------------------------------------

Calculating the individual forecast methods for the 100000 time series in the M4 dataset is a very time consuming process (hours). One of the individual methods `nnetar` in the `forecast` R package produces random results which difficult reproducibiliy, for instance when running in a cluster. Training the learning model is also time consuming.

In order to facilitate the exact reproduction of the process, the individual forecast are already calculated in the (`M4metaresults`package)\[<https://github.com/pmontman/M4metaresults>\] and we show here how to apply our approach to the dataset of precalculated individual forecasts.

``` r
library(M4metalearning)
library(M4metaresults)
data <- create_feat_classif_problem(submission_M4)

preds <- predict_selection_ensemble(model_M4, data$data))

replication_M4 <- ensemble_forecast(preds, submission_M4)


interval_M4 <- predict_interval(submission_M4,
                                get_M4_interval_weights())
```

By examining the `y_hat` elements of the `replication_M4` list we can check the actual mean forecasts. By examining the `upper` and `lower` elements of the interval\_M4 list we can check the upper and lower prediction interval bounds. \*\*NOTE that these are the exact results submitted to the M4 competition, in the same series if forecats using the `forecats_meta_M4` approach of the previous subsection, the resulst will be slightly different due to the randomness of the `nnetar` method.

Reproducing the Training of the Model
=====================================

This section explains how to arrive to the trained metalearning model, for producing the mean forecast and the prediction intervals. Note that this is a very time consuming process, the pretrained models have been provided in the previous section.

``` r

#training
library(M4metalearning)

#create the training set using temporal holdout
set.seed(10-06-2018)
meta_M4 <- temp_holdout(M4comp2018::M4)


#calculate the forecasts of each method in the pool
#THIS WILL TAKE A LOT OF TIME (hours...)
meta_M4 <- calc_forecasts(meta_M4, forec_methods(), n.cores=3)
#calculate the OWA errors
meta_M4 <- calc_errors(meta_M4)
#extract the features
meta_M4 <- THA_features(meta_M4, n.cores=3)

#search for hyperparameters
hyperparameter_search(meta_M4, filename = "M4_hyper.RData", n_iter=150)

#get the best hyperparameter found
load("M4_hyper.RData")
best_hyper <- bay_results[ which.min(bay_results$combi_OWA), ]

#Train the metalearning model with the best hyperparameters found

train_data <- create_feat_classif_problem(meta_M4)
param <- list(max_depth=best_hyper$max_depth,
              eta=best_hyper$eta,
              nthread = 3,
              silent=1,
              objective=error_softmax_obj,
              num_class=ncol(train_data$errors), #the number of forecast methods used
              subsample=bay_results$subsample,
              colsample_bytree=bay_results$colsample_bytree)

meta_model <- train_selection_ensemble(train_data$data,
                                       train_data$errors,
                                       param=param)

## Now the model is trained, lest produce the predictions

final_M4 <- M4comp2018::M4

#just calculate the forecast and features
final_M4 <- calc_forecasts(final_M4, forec_methods())
final_M4 <- THA_features(final_M4)

#get the feature matrix
final_data <- create_feat_classif_problem(final_M4)
#calculate the predictions using our model
preds <- predict_selection_ensemble(meta_model, final_data$data)
#calculate the final mean forecasts
final_M4 <- ensemble_forecast(preds, final_M4)
#the combination predictions are in the field y_hat of each element in the list
#lets check one
final_M4[[1]]$y_hat
}
```

For the interval, a similar procedure is required.

``` r

#we need the mean predictions on the training as the centers of the intervals
#we will calculate them
preds <- predict_selection_ensemble(meta_model, train_data$data)
#this will add the mean forecasts to the dataset, not called interval_M4
interval_M4 <- ensemble_forecast(preds, meta_M4)
#we calculate the radius of the individual forecasts methods
interval_M4 <- calc_radius_interval(interval_M4)
#prepare the data to pose it as a minimization problem
info <- prepare_radius_info(interval_M4)

set.seed("10-06-2018")
#calculate the weights for the combinations of individual methods radius
res_opt <- train_interval_weights(interval_M4, 48)

#transform the results of the optimization process to the format used by predict interval
weights <- lapply(res_opt, function (lentry) lentry$opt$par)

#in res opt we have the 
train_dataset <- predict_interval(train_dataset, weights, TRUE)
```
