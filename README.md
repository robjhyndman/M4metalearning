
<!-- README.md is generated from README.Rmd. Please edit that file -->
M4Metalearning
==============

The goal of M4Metalearning is to provide a forecasting metalearning tools, see the example for more information.

Installation
------------

You can install M4Metalearning from github with:

``` r
# install.packages("devtools")
devtools::install_github("robjhyndman/M4metalearning", auth_token = "AUTH_TOKEN_REQUIRED_WHILE_REPOSITORY_IS_PRIVATE")
```

You will need and auth token from github, because the repository is private. See how to create personal tokens here: <https://help.github.com/articles/creating-a-personal-access-token-for-the-command-line/>

Example
-------

The purpose of this package is to apply metalearning to time series forecasting problems. A pretrained metalearning model will be provided, one that can be applied directly to a given time series. Additionally, if we have a large time series dataset, we can also train a new metalearning model for it using this package.

We will see an example of the proposed metalearning method, applied to the M3 and M1 forecasting competition datasets. We will consider the M3 dataset as the 'training set' and fit an ensemble metalearning model to it. We will test this new model on the M1 dataset and see the differences in error between using our ensemble,using a single forecasting method, and using other ensemble method.

We start with a dataset that contains the time series and the desired true forecasts.

The first step is to appy the available forecasting methods to each series in the dataset and keep their forecast, as well as calculate the OWI error they produce. This package provides the function `process_forecast_dataset` that takes as input the trainig dataset with the time series, and a list of forecasting methods to apply to it. We can see in the documentation of the function that the input dataset requires the following format: A `list` with each element having the following structure:

-   `x` : Containing the series as a `ts` object.
-   `h` : An integer with the amount of future moments to forecast.
-   `xx` : The 'true' future time series.

This format is inspired by the one used in `Mcomp` and `M4comp2018` packages. Each element of the list may contain additional fields, that will be simply ignored by `process_forecast_dataset`.

On the other hand, the second parameter, `methods_list` is a list of strings, each string has to be the name of an existing **`R`** function that will produce the forecast. The functions pointed by `methods_list` must have the following structure:

    my_forect_method <- fuction(x, h) { ... }

With x being a `ts` object to forecast and `h` the amount of future moments to forecast. **The output of these functions must be a numeric vector of length `h` containing the forecasts of the method.** The input format is kept for simplicity, any parameters used by the methods must be set internally in these functions.

In this package, some forecasting functions with this format are provided, which are wraps of methods in the `forecast` package.

``` r
library(M4Metalearning)
# see an example of function, just a wrap of a method in the forecast pacakge
auto_arima_forec
#> function(x, h) {
#>   model <- forecast::auto.arima(x, stepwise=FALSE, approximation=FALSE)
#>   forecast::forecast(model, h=h)$mean
#> }
#> <environment: namespace:M4Metalearning>
```

A list of methods wrapped from the `forecast` package is also provided in the function `create_seas_method_list()`.

Without further ado, lets process the M3 dataset to generate all forecast and errors, to be used in the metalearning.

``` r
#this will take time, use the pregenerated result include as data in the package
#forec_M3 <- process_forecast_dataset(Mcomp::M3, create_seas_method_list(), n.cores=3)
data(forec_M3)
```

Now we have in `forec_M3` the series, the forecasts of each method for each series, and the errors of these forecasts. Specifically, `process_forecast_dataset` produces as output a list, with each element having the following structure:

-   x : The series to be forecasted.
-   h : The amount of moments into the future to forecast.
-   xx : The true value of the future of the series.
-   ff : The forecasts of each methods, a matrix with the forecast of each methods per row.
-   errors : A vector with the OWI errors produced by each method.

The idea behind `process_forecast_dataset` is to produce all relevant information that a metalearning method may required for the training.

### A metalearning method: tsfeatures + gradient boosted trees

The information in `forec_M3` is enough to train a metalearning method, but many approaches will required further processing, e.g. extracting features, zero padding the series so all have the same length, etc.

The metalearning method we will show uses features extracted from the series instead of the raw series. With these features, it trains a xgboost model that gives weights to the forecast methods in order to minimize the OWI error (roughly speaking).

The specific features extracted from the series are a key part of the process, and in this package, the function `generate_THA_feature_dataset` calculates such features. A rationale and description of the features may be seen in a --upcoming paper from (Talagala, Hyndman and Athanasopoulos, 2018)--.

`generate_THA_feature_dataset` processes once again a dataset of times series, a list of elements containing at least the field `x` with a `ts` object (in the spirit of the input and output datasets of `process_forecast_dataset` and the `Mcomp` pacakge format). We can use `forec_M3` directly to generate the features for each series.

``` r
#Will take some time, use the pregenerated file include in the package
#feat_forec_M3 <- generate_THA_feature_dataset(forec_M3, n.cores=3)
data(feat_forec_M3)
```

The output of `generate_THA_feature_dataset` is its input list, but to each element, the field `THA_features` has been added. `THA_features` is a `tibble` object with the extracted features. `THA_features` uses the package `tsfeatures` for extracting the features.

Now we have in `feat_forec_M3` the series, its extracted features, forecasts and errors.

The next step is training the ensemble using xgboost the extracted features and the error produced by eac method.

The `create_feat_classif_problem` auxiliary function is provided to reformat the list produced by the previously shown functions to the more common format of a feature matrix and target labels used in classification functions such as `randomForest`, `svm`, `xgboost`, ...

`create_feat_classif_problem` simply produces a list with the entries:

-   data : The features extracted from the series
-   errors : The errors produced by the forecasting method
-   labels : The target classification problem, created by selecting the method that produces the smallest error for each series.

``` r
train_data <- create_feat_classif_problem(feat_forec_M3)
head(train_data$data, n=3)
#>          x_acf1  x_acf10 diff1_acf1 diff1_acf10   diff2_acf1 diff2_acf10
#> N0001 0.7623182 1.504539  0.5974236   0.6308634 -0.004813322   0.1934310
#> N0002 0.7507872 1.456497  0.2399691   0.3201949 -0.398246929   0.3245830
#> N0003 0.7687310 1.555912  0.4461251   0.5683772 -0.211798893   0.3048418
#>       seas_acf1 ARCH.LM crossing_points   entropy flat_spots  arch_acf
#> N0001         0       1               1 0.7729350          3 0.7452971
#> N0002         0       1               3 0.8374974          3 0.4920190
#> N0003         0       1               3 0.8250352          2 0.4863309
#>       garch_acf arch_r2 garch_r2     alpha         beta     hurst
#> N0001 0.4068810       1        1 0.9709143 0.9709142237 0.9710509
#> N0002 0.5131037       1        1 0.9998999 0.0001000064 0.9473065
#> N0003 0.5253350       1        1 0.9998998 0.2635651859 0.9486339
#>       lumpiness nonlinearity   x_pacf5 diff1x_pacf5 diff2x_pacf5 seas_pacf
#> N0001         0     2.124405 0.6152347    0.5483426    0.2301945         0
#> N0002         0     1.998710 0.8093241    0.1565805    0.3074159         0
#> N0003         0     1.449664 0.9173062    0.3708305    0.1717048         0
#>       nperiods seasonal_period     trend        spike linearity curvature
#> N0001        0               1 0.9950394 2.373423e-07  3.583026  0.423830
#> N0002        0               1 0.8687934 1.787445e-04  2.053111 -2.084470
#> N0003        0               1 0.8648297 1.933079e-04  1.751751 -2.256782
#>          e_acf1   e_acf10 seasonal_strength peak trough stability hw_alpha
#> N0001 0.4124236 1.0452773                 0    0      0         0        0
#> N0002 0.3240316 0.4849437                 0    0      0         0        0
#> N0003 0.4571183 0.7723079                 0    0      0         0        0
#>       hw_beta hw_gamma unitroot_kpss unitroot_pp series_length
#> N0001       0        0     0.5757141    1.329299            14
#> N0002       0        0     0.2965073   -3.735398            14
#> N0003       0        0     0.2542136   -3.978590            14
head(train_data$errors, n=3)
#>       auto_arima_forec ets_forec nnetar_forec tbats_forec stlm_ar_forec
#> N0001         0.186718 0.1863118    0.2709279   0.1309829      0.322879
#> N0002         1.000000 0.5328959    0.1564174   0.8564553      1.253115
#> N0003         1.713569 0.9999183    6.5455414   0.8393168      1.255951
#>       rw_drift_forec naive_forec snaive_forec
#> N0001      0.5188943           1            1
#> N0002      0.3518504           1            1
#> N0003      3.3038686           1            1
head(train_data$labels, n=3)
#> N0001 N0002 N0003 
#>     3     2     3
```

The data in this format is easy to use with any classifier, as we will see. The next step is training the metalearning model: a 'classifier' that selects the best forecasting method for a time series, given its extracted features. One of the proposals in this package is `train_selection_ensemble`, that will create a `xgboost` classifier, but it trains it with a custom objective function that requires the whole errors information instead of only which method is the best for each series.

``` r
set.seed(1345) #set the seed because xgboost is random!
meta_model <- train_selection_ensemble(train_data$data, train_data$errors, train_data$labels)
```

Now we have in `meta_model` a `xgb.Booster` object that can be uses indepently, but easy to use functions for prediction and performance measurement are also provided in the package. It only remains now to test this model with the M1 competition dataset. We 'really' only need to extract the features from the test dataset, but we will generate the forecasts and errors for the performance analysis.

``` r
#Processing will take some time, use the pregenerated data in the package
#forec_M1 <- process_forecast_dataset(Mcomp::M1, create_seas_method_list(), n.cores=3)
#feat_forec_M1 <- generate_THA_feature_dataset(forec_M1, n.cores=3)
data("feat_forec_M1")
#pose it as a classification problem for ease of use
test_data <- create_feat_classif_problem(feat_forec_M1)
pred <- predict_selection_ensemble(meta_model, test_data$data)
head(pred)
#>           [,1]       [,2]       [,3]       [,4]        [,5]      [,6]
#> [1,] 0.2269338 0.12371762 0.03105979 0.04737823 0.003582103 0.3412190
#> [2,] 0.1002761 0.04629103 0.03126886 0.08758427 0.002947387 0.1507138
#> [3,] 0.1459839 0.06405401 0.02346983 0.05545320 0.001838499 0.3958520
#> [4,] 0.3362676 0.11876423 0.04596597 0.11015412 0.004861711 0.1601122
#> [5,] 0.1502393 0.10872702 0.02980672 0.04845675 0.004260490 0.4099174
#> [6,] 0.2278793 0.11696084 0.03838183 0.08817677 0.003558375 0.2836177
#>            [,7]       [,8]
#> [1,] 0.09368592 0.13242356
#> [2,] 0.45994448 0.12097406
#> [3,] 0.11652602 0.19682251
#> [4,] 0.15561923 0.06825489
#> [5,] 0.12468542 0.12390682
#> [6,] 0.11220496 0.12922026
```

The output of `predict_selection_ensemble` produces probabilities for each class instead of just the selected class (which would be the one with max probability, usually ;) ) To show the performance of the metalearning model, we have the function `summary_performance`, that requires as input the class probabilities, the errors and labels produced on the test set.

``` r
summary_performance(pred, test_data$errors, test_data$labels)
#> [1] "Classification error:  0.8422"
#> [1] "Selected OWI :  0.8992"
#> [1] "Oracle OWI:  0.5807"
#> [1] "Single method OWI:  0.934"
#> [1] "Average OWI:  1.101"
```

Before commenting on the results, lets compare with other 'state-of-the-art' classifier, to also showcase the ease of use. We will use a basic `randomForest` classifier for selecting the method.

``` r
#train
rforest <- randomForest::randomForest(train_data$data, y=as.factor(train_data$labels))
#test
pred_forest <- predict(rforest, newdata=test_data$data, type="prob")
#show performance
summary_performance(pred_forest, test_data$errors, test_data$labels)
#> [1] "Classification error:  0.7243"
#> [1] "Selected OWI :  0.9468"
#> [1] "Oracle OWI:  0.5807"
#> [1] "Single method OWI:  0.934"
#> [1] "Average OWI:  1.101"
```

The important output of `summary_performance` is 'Classification error', which is the error in the classification problem of selecting the best forecasting method, and 'Selected OWI' which is the average OWI error produced by the methods selected by the classifier **The important measure!**. 'Oracle OWI' shows the theoretical minimum error that a classifier would produce, 'Single method OWI' is the best method in our pool of forecasting methods, and 'Average OWI' would be the error produced by selecting methods at random from out pool of methods for each series.

We see that while randomForest produces way better classification error, the real forecasting error of the selected methods is worse than our proposal.

### The End!

<!-- Vignettes are long form documentation commonly included in packages. Because they are part of the distribution of the package, they need to be as compact as possible. The `html_vignette` output type provides a custom style sheet (and tweaks some options) to ensure that the resulting html is as small as possible. The `html_vignette` format: -->
<!-- - Never uses retina figures -->
<!-- - Has a smaller default figure size -->
<!-- - Uses a custom CSS stylesheet instead of the default Twitter Bootstrap style -->
<!-- ## Vignette Info -->
<!-- Note the various macros within the `vignette` section of the metadata block above. These are required in order to instruct R how to build the vignette. Note that you should change the `title` field and the `\VignetteIndexEntry` to match the title of your vignette. -->
<!-- ## Styles -->
<!-- The `html_vignette` template includes a basic CSS theme. To override this theme you can specify your own CSS in the document metadata as follows: -->
<!--     output:  -->
<!--       rmarkdown::html_vignette: -->
<!--         css: mystyles.css -->
<!-- ## Figures -->
<!-- The figure sizes have been customised so that you can easily put two images side-by-side.  -->
<!-- ```{r, fig.show='hold'} -->
<!-- plot(1:10) -->
<!-- plot(10:1) -->
<!-- ``` -->
<!-- You can enable figure captions by `fig_caption: yes` in YAML: -->
<!--     output: -->
<!--       rmarkdown::html_vignette: -->
<!--         fig_caption: yes -->
<!-- Then you can use the chunk option `fig.cap = "Your figure caption."` in **knitr**. -->
<!-- ## More Examples -->
<!-- You can write math expressions, e.g. $Y = X\beta + \epsilon$, footnotes^[A footnote here.], and tables, e.g. using `knitr::kable()`. -->
<!-- ```{r, echo=FALSE, results='asis'} -->
<!-- knitr::kable(head(mtcars, 10)) -->
<!-- ``` -->
<!-- Also a quote using `>`: -->
<!-- > "He who gives up [code] safety for [code] speed deserves neither." -->
<!-- ([via](https://twitter.com/hadleywickham/status/504368538874703872)) -->
