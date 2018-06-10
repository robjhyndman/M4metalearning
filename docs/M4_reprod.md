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
devtools::install_github("rbojhyndman/M4metadata")
devtools::install_github("pmontman/M4metaresults")
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
times =system.time(
preds <- predict_selection_ensemble(model_M4, data$data))
times

replication_M4 <- ensemble_forecast(preds, submission_M4)


interval_M4 <- predict_interval(submission_M4,
                                get_M4_interval_weights())
```

Reproducing the Training of the Model
=====================================
