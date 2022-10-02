.rxModelVarsLast <- NULL

#' Internal translation to get model variables list
#'
#' 
#' @param model Model (either file name or string)
#' @export
#' @useDynLib rxode2parse, .registration=TRUE
#' @examples
#' .rxTrans("a=3")
rxode2parse <- function(model) {
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
