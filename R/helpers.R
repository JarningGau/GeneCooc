#' @export
softmax <- function(x) exp(x) / sum(exp(x))
