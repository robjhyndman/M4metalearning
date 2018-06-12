#' @export
opt_msis <- function(w, ff, radius, xx, ew) {
  w <- t(replicate(nrow(radius), w))
  combi_radius <- rowSums(w * radius)

  err_up <- xx - (ff + combi_radius)
  err_low <- (ff - combi_radius) - xx
  sum((2*combi_radius + 40*err_up*(err_up>0) + 40*err_low*(err_low>0))*ew)
}

#' @export
prepare_radius_info <- function(dataset) {
  ff <- lapply(dataset, function (l) l$y_hat)
  ff <- do.call(c,ff)

  radius <- lapply(dataset, function (l) t(l$radius))
  radius <- do.call(rbind,radius)

  info <- create_combi_info(dataset)
  info$ff <- ff
  info$radius <- radius
  info
}


#' @export
subinfo_h <- function(info, dataset, h) {
  subinfo <- NULL
  subinfo$ew <- NULL
  subinfo$ff <- NULL
  subinfo$xx <- NULL
  subinfo$radius <- NULL
  counth = 1
  for (i in 1:length(dataset)) {
    lh <- dataset[[i]]$h
    if (lh < max(h)) next;
    indiceh <- counth:(counth + lh -1)
    subinfo$ff <- c(subinfo$ff, info$ff[indiceh][h])
    subinfo$ew <- c(subinfo$ew, info$ew[indiceh][h])
    subinfo$xx <- c(subinfo$xx, info$xx[indiceh][h])
    subinfo$radius <- rbind(subinfo$radius, info$radius[indiceh,][h,])
    counth = counth + lh
  }
  subinfo
}

#' @export
train_interval_weights <- function(dataset, maxh=48) {

  train_dataset <-dataset
  info <- prepare_radius_info(train_dataset)

  res_perh <- NULL
  for (h in 1:maxh) {
    print(h)
    sinfo <- subinfo_h(info, train_dataset, h)

    w <- rnorm(3)

    opt <- optim(w, opt_msis, ff=sinfo$ff, radius=sinfo$radius,
                 xx=sinfo$xx, ew=sinfo$ew, method="CG", control=list(trace=0))

    res_perh <- append(res_perh, list(list(h=h, opt=opt)) )
  }
  res_perh
}

#' @export
predict_interval <- function(dataset, weights, clamp_zero=TRUE) {
    for (i in 1:length(dataset)) {
      radius = rep(0,dataset[[i]]$h)
      for (j in 1:dataset[[i]]$h) {
        radius[j] = sum(weights[,j] * dataset[[i]]$radius[,j])
      }
      upper <- dataset[[i]]$y_hat + radius
      lower <- dataset[[i]]$y_hat - radius
      dataset[[i]]$upper <- upper
      dataset[[i]]$lower <- lower
    }
    dataset
  }


#' @export
calc_radius_interval <- function(dataset) {
  for (i in 1:length(dataset)) {
    thetamod <- forecast::thetaf(dataset[[i]]$x, dataset[[i]]$h)
    radius_theta <- thetamod$upper[,2] - thetamod$mean

    snaivemod <- forecast::snaive(dataset[[i]]$x, dataset[[i]]$h)
    radius_snaive <- snaivemod$upper[,2] - snaivemod$mean

    naivemod <-forecast::naive(dataset[[i]]$x, dataset[[i]]$h)
    radius_naive <- naivemod$upper[,2] - naivemod$mean

    dataset[[i]]$radius <- rbind(radius_theta, radius_naive, radius_snaive)
  }
  dataset
}
