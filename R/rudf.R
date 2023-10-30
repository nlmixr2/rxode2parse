.udfEnv <- new.env(parent=emptyenv())
.udfEnv$fun <- list()
.udfEnv$udf <- integer(0)

#' While parsing or setting up the solving, get information about the
#' user defined function
#'
#' @param fun function (character) to get information about
#' @return A list with two elements
#'   - nargs = `NA` if the user function isn't supported, or the number of arguments suported
#'   - string = Error message when `NA` or function string
#' @noRd
#' @author Matthew L. Fidler
.getUdfInfo <- function(fun) {
  .fun <- try(get(fun, mode="function"), silent=TRUE)
  if (inherits(.fun, "try-error")) {
    return(list(nargs=NA_integer_,
                sprintf("function '%s' is not supported; user function not found",
                        fun)))
  }
  .formals <- formals(.fun)
  if (any(names(.formals) == "...")) {
    return(list(nargs=NA_integer_,
                "rxode2 user defined R cannot have '...' arguments"))
  }
  .nargs <- length(.formals)
  .udfEnv$fun[[fun]] <- list(.fun, environment(.fun))
  .udfEnv$udf <- c(.udfEnv$udf, setNames(.nargs, fun))
  return(list(nargs=.nargs,
              fun))
}

#' This function is run before starting a rxode2 solve to make sure
#' the R-based user functions are setup correctly.
#'
#' This function also resets the udf-based run-time errors
#'
#' @param iv Named Integer Vector with the names representing the
#'   functions and the integers representing the number of arguments
#'   that were present when the model was compiled
#' @return nothing, called for side effect
#' @noRd
#' @author Matthew L. Fidler
.setupUdf <- function(iv) {
  .n <- names(iv)
  lapply(.n,
         function(n) {
           .oldArg <- iv[n]
           .new <- .getUdfInfo(n)
           if (is.na(.new[[1]])) {
             stop(.new[[2]], call.=FALSE)
           } else if (.new[[1]] != .oldArg) {
             stop("'", n,
                  "' had ", .oldArg, " arguments when model was compiled, now it has ",
                  .new[[1]], " arguments",
                  call.=FALSE)
           }
           NULL
         })
}
#' Reset the tracking of user defined functions
#'
#' This is called during parsing reset
#'
#' @return Nothing, called for side effects
#' @noRd
#' @author Matthew L. Fidler
.udfReset <- function() {
  .udfEnv$udf <- integer(0)
}

#' This gets the user defined functions information for incorporation
#' in the model variables
#'
#' @return A integer vector; The values are the number of arguments;
#'   the names are the function names
#' @author Matthew L. Fidler
#' @noRd
.udfInfo <- function() {
  .udfEnv$udf
}
#' Get the function name with the current arguments as a string
#'
#' @param fun function name
#' @param args  arguments
#' @return string of the form 'fun(arg1, arg2)':
#' @export
#' @author Matthew L. Fidler
#' @examples
.udfCallFunArg <- function(fun, args) {
  paste0("'", fun, "(",
         paste(vapply(seq_along(args),
                function(i) {
                  as.character(args[[i]])
                }, character(1), USE.NAMES=FALSE), collapse=", "),
         ")': ")
}
#' This is the function that is always called for every user function in rxode2
#'
#' @param fun A character vector representing the function
#' @param args A list of double numbers that will be used as the
#'   function arguments
#' @return A double numeric value, including `NA_real` when the
#'   function isn't working as expected
#' @noRd
#' @author Matthew L. Fidler
.udfCall <- function(fun, args) {
  .info <- .udfEnv$fun[[fun]]
  .fun <- .info[[1]]
  .envir <- .info[[2]]
  .env <- new.env(parent=.envir)
  .env$.fun <- .fun
  .env$.args <- args
  .ret <- try(with(.env, do.call(.fun, .args)), silent=TRUE)
  if (inherits(.ret, "try-error")) {
    .msg <- try(attr(.ret, "condition")$message, silent=TRUE)
    if (inherits(.msg, "try-error")) .msg <- "Unknown Error"
    # This can error since it isn't threaded
    stop(paste0(.udfCallFunArg(fun, args), .msg), call.=FALSE)
  }
  if (length(.ret) != 1L) {
    # This can error since it isn't threaded
    stop(paste0(.udfCallFunArg(fun, args), "needs to return a length 1 numeric"),
         call.=FALSE)
  }
  .ret <- try(as.double(.ret), silent=TRUE)
  if (inherits(.ret, "try-error")) {
    .msg <- try(attr(.ret, "condition")$message, silent=TRUE)
    if (inherits(.msg, "try-error")) .msg <- "Unknown Error"
    stop(paste0(.udfCallFunArg(fun, args), .msg), call.=FALSE)
  }
  .ret
}
