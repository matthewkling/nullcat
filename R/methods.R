
# Single source of truth for supported methods
NULLCAT_METHODS <- c("curvecat", "swapcat", "tswapcat", "r0cat", "c0cat")


#' Supported nullcat methods
#'
#' Return the character vector of supported categorical randomization methods.
#'
#' @return A character vector of method names.
#' @examples
#' nullcat_methods()
#'
#' @export
nullcat_methods <- function() NULLCAT_METHODS
