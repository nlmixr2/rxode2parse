.rxModelVarsLast <- NULL

#' Internal translation to get model variables list
#'
#' 
#' @param model Model (either file name or string)
#' @export
#' @useDynLib rxode2parse, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom qs qsave
#' @importFrom dparser dparse
#' @importFrom utils capture.output
#' @eval rxode2parseFuns()
#' @examples
#' rxode2parse("a=3")
rxode2parse <- function(model) {
  rxParseSuppressMsg()
  checkmate::assertCharacter(model, len=1, any.missing=FALSE)
  modelPrefix=""
  md5=""
  meCode=""
  fullPrint=FALSE
  if (file.exists(model)) {
    .isStr <- 0L
  } else {
    .isStr <- 1L
  }
  if (missing(md5)) {
    md5 <- "none"
  }
  .ret <- .Call(`_rxode2parse_trans`, model, modelPrefix, md5, .isStr,
    as.integer(crayon::has_color()),
    meCode, .parseFuns,
    fullPrint)
  .ret
}

rxode2parseFuns <- function() {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("this requires devtools", call.=FALSE)
  }
  message("rebuild parseFuns.R from rxode2")
  try(source(devtools::package_file("build/refresh.R")))
  message("done")
  ""
}
