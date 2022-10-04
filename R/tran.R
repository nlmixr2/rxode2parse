.rxModelVarsLast <- NULL

#' Internal translation to get model variables list
#'
#' 
#' @param model Model (either file name or string)
#' @return A rxModelVars object that has the model variables of a rxode2 syntax expression
#' @export
#' @useDynLib rxode2parse, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom qs qsave
#' @importFrom dparser dparse
#' @importFrom utils capture.output assignInMyNamespace
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

#' This assigns the c level linkages for a roxde2 model
#'
#' @param df data frame containing the character column names rxFun,
#'   fun, type, package, packageFun and the integer column names
#'   argMin and argMax
#' @return Nothing called for side effects
#' @author Matthew L. Fidler
#' @export 
#' @examples
#' 
#' rxode2parseAssignTranslation(rxode2parseGetTranslation())
#' 
rxode2parseAssignTranslation <- function(df) {
  .char <- c("rxFun", "fun", "type", "package", "packageFun")
  .int <- c("argMin", "argMax", "threadSafe")
  .df <- df[,c(.char, .int)]
  for (.c in .char) {
    .df[[.c]] <- as.character(.df[[.c]])
  }
  for (.i in .int) {
    .df[[.i]] <- as.integer(.df[[.i]])
  }
  assignInMyNamespace(".rxode2parseDf", .df)
  invisible(.df)
}

#' This function gets the currently assigned translations
#' 
#' @return The currently assigned translations
#' @author Matthew L. Fidler
#' @export 
#' @examples
#' rxode2parseGetTranslation()
rxode2parseGetTranslation <- function() {
  .rxode2parseDf
}
#' Get the MD5 hash of the current language revision
#'
#' @return md5 hash of language revision
#' @author Matthew L. Fidler
#' @export 
#' @examples
#' rxode2parseMd5()
rxode2parseMd5 <- function() {
  rxode2parse.md5  
}
