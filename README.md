
<!-- README.md is generated from README.Rmd. Please edit that file -->
This page contains the explanation of our forecast method for the M4 competition, authored by Pablo Montero-Manso, Thiyanga Talagala, Rob J Hyndman and George Athanasopoulos.

As part of our submission, we are producing the R package `M4metalearning`, see instalation instructions below. This package is intended to facilitate reproducing the results of our submission, but can also the used on its own to apply our approach to other datasets, either using the pretrained model submitted to the M4 competition or by training new models.

Additionaly, the authors have produced another R package `M4comp2018`, which facilitates the users the access to the M4 competition dataset. [The M4comp2018 package can be found here](https://github.com/carlanetto/M4comp2018)

The description is divided into three sections, the first t

1.  [Methodology](docs/M4_methodology.html)
2.  [Reproducing the results](docs/M4_reprod.md)
3.  [Usage example of the package](docs/metalearning_example.md)

M4metalearning
==============

The goal of the `M4Metalearning` package is to provide a forecasting metalearning tools based on our submission to the M4 competition.

Installation
------------

You can install `M4metalearning` from github with:

``` r
# install.packages("devtools")
devtools::install_github("robjhyndman/M4metalearning)
```

### Note

M4metalearning is using, for the time being, a slight modification of the `tsfeatures` package. Please install it from:

``` r
# install.packages("devtools")
devtools::install_github("pmontman/tsfeatures")
```

Also, a custom version of the `xgboost` package is required. It is installed automatically when calling the training functions that use it, since it does not break compatibility. It just supports customized multiclass objective functions. You may install it manually from:

``` r
# install.packages("devtools")
devtools::install_github("pmontman/customxgboost")
```

### Usage

For an example of the usage of the package [see this page](docs/metalearning_example.md)
