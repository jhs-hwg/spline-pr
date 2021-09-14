

#' Create bins for histogram
#'
#' @param x
#' @param y
#' @param x_min
#' @param x_max
#' @param bin_count
#' @param bin_yintercept
#' @param bin_length
#'
#' @return
#' @export
#'
#' @examples
#'

bin_segments <- function(x,
                         x_min = NULL,
                         x_max = NULL,
                         y = NULL,
                         by_y = FALSE,
                         bin_count,
                         bin_yintercept,
                         bin_length){

  bins <- seq(from = if(is.null(x_min)) min(x) else x_min,
              to   = if(is.null(x_max)) max(x) else x_max,
              length.out = bin_count + 1)

  if(!by_y){

    freqs <- table(cut(x, bins))
    freqs <- 0.1 * freqs / max(freqs)

    return(
      data.frame(
        x = bins[-(bin_count + 1)],
        y = bin_yintercept,
        xend = bins[-(bin_count + 1)],
        yend = as.numeric(bin_length * freqs + bin_yintercept)
      )
    )

  }

  if(is.null(y)) stop("y must be specified if by_y = TRUE")

  # separate frequency counts by event status
  f0	<- table(cut(x[y == 0], bins))
  f1	<- table(cut(x[y == 1], bins))

  j0	<- f0 > 0
  j1	<- f1 > 0

  bins0 <- (bins[-(bin_count + 1)])[j0]
  bins1 <- (bins[-(bin_count + 1)])[j1]

  f0	<- f0[j0]
  f1	<- f1[j1]

  maxf <- max(f0, f1)

  f0	<- (0.1 * f0) / maxf
  f1	<- (0.1 * f1) / maxf

  data_segments_below <- data.frame(
    x = bins0,
    y = bin_yintercept,
    xend = bins0,
    yend = as.numeric(-1 * bin_length * f0 + bin_yintercept),
    event_status = 0
  )

  data_segments_above <- data.frame(
    x = bins1,
    y = bin_yintercept,
    xend = bins1,
    yend = as.numeric(bin_length * f1 + bin_yintercept),
    event_status = 1
  )

  rbind(data_segments_above, data_segments_below)

}

#' Get index of pattern in coefficient names
#'
#' @param fit
#' @param pattern
#'
#' @return
#' @export
#'
#' @examples
get_index <- function(fit, pattern){
 which(stringr::str_detect(names(coef(fit)), pattern = pattern))
}



#' variance covariance for geeglm objects
#'
#' simple wrapper on summary.geeglm to extract covariance matrix.
#' A use for this function is pooling regression coefficients.
#'
#' @param object a geeglm object
#' @param ... currently not used
#'
#' @return covariance matrix for regression coefficients in `object`
vcov.geeglm <- function(object, ...){
 getElement(summary(object), 'cov.scaled')
}

#' get spline predictions
#'
#' @param fit
#' @param basis
#' @param pattern
#' @param x_min
#' @param x_max
#' @param x_ref
#'
#' @return
#' @export
#'
#' @examples
get_spline_preds <- function(fit,
                             basis,
                             pattern,
                             x_min,
                             x_max,
                             x_ref = NULL){

 spline_index <- get_index(fit, pattern = pattern)
 spline_coef <- coef(fit)[spline_index]
 spline_covariance <- vcov(fit)[spline_index, spline_index, drop = FALSE]

 x_grid <- seq(x_min, x_max, length.out = 1000)
 x_mat <- predict(basis, newx = x_grid)

 if(!is.null(x_ref)){
  ref_index <- which.min(abs(x_grid - x_ref))
  x_mat <- sweep(x_mat, MARGIN = 2, STATS = x_mat[ref_index, ])
 }

 xb <- x_mat %*% spline_coef
 xb_stderr <- sqrt(diag(x_mat %*% tcrossprod(spline_covariance, x_mat)))
 lwr_ci <- xb - qt(0.975, fit$df.residual+1) * xb_stderr
 upr_ci <- xb + qt(0.975, fit$df.residual+1) * xb_stderr

 if(requireNamespace('tibble', quietly = TRUE)){
  return(
   tibble::tibble(
    x      = x_grid,
    pred   = as.numeric(xb),
    se     = as.numeric(xb_stderr),
    ci_lwr = as.numeric(lwr_ci),
    ci_upr = as.numeric(upr_ci)
   )
  )
 }

 data.frame(
  x      = x_grid,
  pred   = as.numeric(xb),
  se     = as.numeric(xb_stderr),
  ci_lwr = as.numeric(lwr_ci),
  ci_upr = as.numeric(upr_ci)
 )

}
