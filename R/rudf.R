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
    return(NA_real_)
  }
  if (length(.ret) != 1L) {
    return(NA_real_)
  }
  .ret <- try(as.double(.ret), silent=TRUE)
  if (inherits(.ret, "try-error")) {
    return(NA_real_)
  }
  .ret
}
