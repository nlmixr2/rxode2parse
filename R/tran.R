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
    meCode, .parseEnv$.parseFuns,
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
  assign(".rxode2parseDf", .df, envir=.parseEnv)
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
  .parseEnv$.rxode2parseDf
}

.parseEnv$.packagesToLoad <- c("rxode2ll", "rxode2parse", "rxode2random", "rxode2et")

#'@rdname rxode2parseAssignPackagesToLoad
#'@export
rxode2parseGetPackagesToLoad <- function() {
  .parseEnv$.packagesToLoad
}
  
#' Control the packages that are loaded when a `rxode2` model dll is loaded
#'
#' @param pkgs The packages to make sure are loaded every time you load an rxode2 model.
#' @return List of packages to load
#' @author Matthew Fidler
#' @examples
#'
#' rxode2parseGetPackagesToLoad()
#'
#' rxode2parseAssignPackagesToLoad(rxode2parseGetPackagesToLoad())
#' @export
rxode2parseAssignPackagesToLoad <- function(pkgs=rxode2parseGetPackagesToLoad()) {
  assign(".packagesToLoad", pkgs, envir=.parseEnv)
  pkgs
}


.parseEnv$.rxode2parsePointerAssignment <- "rxode2parse"

#' This function gets the currently assigned function pointer assignments
#' 
#' @return The currently assigned pointer assignments
#' @author Matthew L. Fidler
#' @export 
#' @examples
#' rxode2parseGetTranslation()
rxode2parseGetPointerAssignment <- function() {
  .parseEnv$.rxode2parsePointerAssignment
}

#' This sets function gets the currently assigned function pointer assignments
#'
#' @param var List of packages where pointer assignment will be called.
#' 
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @keywords internal
#' @export 
#' @examples
#' rxode2parseAssignPointerTranslation("rxode2parse")
rxode2parseAssignPointerTranslation <- function(var) {
  checkmate::assertCharacter(var)
  assign(".rxode2parsePointerAssignment", var, envir=.parseEnv)
  invisible()
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
