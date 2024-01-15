#' Softmax Function
#'
#' The softmax function, also known as the normalized exponential function, is a generalization
#' of the logistic function that squashes a K-dimensional vector of arbitrary real values
#' to a K-dimensional vector of real values in the range (0, 1) that add up to 1.
#' It is commonly used in classification problems.
#'
#' @param x A numeric vector or array.
#' @return A numeric vector representing the softmax transformation of the input.
#' @export
softmax <- function(x) exp(x) / sum(exp(x))
