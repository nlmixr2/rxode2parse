.udfEnv <- new.env(parent=emptyenv())
.udfEnv$fun <- list()
.udfEnv$udf <- integer(0)
.udfEnv$envir <- new.env(parent=emptyenv())
.udfEnv$lockedEnvir <- FALSE
.udfEnv$rxSEeqUsr <- NULL
.udfEnv$rxCcode <- NULL
.udfEnv$symengineFs <- new.env(parent = emptyenv())
.udfEnv$extraCnow <- ""

#' Get the udf strings for creating model md5
#'
#' @return string vector
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
.udfMd5Info <- function() {
  .tmp <- ls(.udfEnv$symengineFs)
  .env <- new.env(parent=emptyenv())
  .env$found <- FALSE
  .ret <- vapply(.tmp, function(x) {
    .cur <- .udfEnv$fun[[x]]
    if (!is.null(.cur)) {
      .env$found <- TRUE
    }
    x
  }, character(1), USE.NAMES = FALSE)
  if (.env$found) {
    .ret <- c(.ret, data.table::address(.udfEnv$envir))
  }
  .ret
}

#' Generate extraC information for rxode2 models
#'
#' @param extraC Additional extraC from rxode2 compile optioioins
#' @return Nothing, called for side effects
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
.extraC <- function(extraC = NULL) {
  if (!is.null(extraC)) {
    if (file.exists(extraC)) {
      .ret <- sprintf("#include \"%s\"\n", extraC)
    } else {
      .ret <- paste(extraC, collapse = "\n")
    }
  } else {
    .ret <- ""
  }
  if (length(.udfEnv$rxCcode) > 0L) {
    .ret <- sprintf("%s\n%s\n", .ret, paste(.udfEnv$rxCcode, collapse = "\n"))
  }
  .udfEnv$extraCnow <- .ret
  return(invisible())
}
#' Get the extraCnow for compiling
#'
#'
#' @return string of extraC information
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
.extraCnow <- function() {
  .udfEnv$extraCnow
}

#' Add user function to rxode2
#'
#' This adds a user function to rxode2 that can be called.  If needed,
#' these functions can be differentiated by numerical differences or
#' by adding the derivatives to rxode2's internal derivative table
#' with rxode2's `rxD` function
#'
#' @param name This gives the name of the user function
#' @param args This gives the arguments of the user function
#' @param cCode This is the C-code for the new function
#' @return nothing
#' @author Matthew L. Fidler
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
rxFunParse <- function(name, args, cCode) {
  if (!is.character(name) || length(name) != 1L) {
    stop("name argument must be a length-one character vector", call. = FALSE)
  }
  if (missing(cCode)) stop("a new function requires a C function so it can be used in rxode2", call. = FALSE)
  if (any(name == names(.udfEnv$rxSEeqUsr))) {
    stop("already defined user function '", name, "', remove it fist ('rxRmFun')",
         call. = FALSE
         )
  }
  suppressWarnings(rxRmFunParse(name))
  .udfEnv$rxSEeqUsr <- c(.udfEnv$rxSEeqUsr, setNames(length(args), name))
  .udfEnv$rxCcode <- c(.udfEnv$rxCcode, setNames(cCode, name))
  assign(name, symengine::Function(name), envir = .udfEnv$symengineFs)
  return(invisible())
}
#' Return the equivalents symengine user functions from C
#'
#' @return equivalent symengine user functions
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
.rxSEeqUsr <- function() {
  .udfEnv$rxSEeqUsr
}

#' Return symengineFs from user functions
#'
#' @return symengineFs from user functions
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
.symengineFs <- function() {
  .udfEnv$symengineFs
}

#' @rdname rxFunParse
#' @export
rxRmFunParse <- function(name) {
  if (!is.character(name) || length(name) != 1L) {
    stop("name argument must be a length-one character vector",
         call. = FALSE)
  }
  if (!any(name == names(.udfEnv$rxSEeqUsr))) {
    warning("no user function '", name, "' to remove", call. = FALSE)
  }
  .w <- which(name == names(.udfEnv$rxSEeqUsr))
  if (length(.w) == 1L) {
    .udfEnv$rxSEeqUsr <- .udfEnv$rxSEeqUsr[-.w]
  }
  .w <- which(name == names(.udfEnv$rxCcode))
  if (length(.w) == 1L) {
    .udfEnv$rxCcode <- .udfEnv$rxCcode[-.w]
  }
  .rxD <- rxode2parse::rxode2parseD()
  if (exists(name, envir = .rxD)) {
    rm(list = name, envir = .rxD)
  }
  if (exists(name, envir = .udfEnv$symengineFs)) {
    rm(list = name, envir = .udfEnv$symengineFs)
  }
  return(invisible())
}
#' Setup the UDF environment (for querying user defined funtions)
#'
#' @param env environment where user defined functions are queried
#' @return nothing called for side effects
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
.udfEnvSet <- function(env) {
  if (.udfEnv$lockedEnvir) return(invisible())
  if (is.environment(env)) {
    .udfEnv$envir <- env
    return(invisible())
  }
  stop("'env' needs to be an environment")
}
#' Lock/Unlock environment for getting R user functions
#'
#'
#' @param lock logical to see if environment to look for user defined
#'   functions is locked.  If it is locked then environments are not
#'   assigned.  When NULL returns lock status.
#' @return lock status
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
.udfEnvLock <- function(lock=TRUE) {
  if (is.null(lock)) return(invisible(.udfEnv$lockedEnvir))
  .udfEnv$lockedEnvir <- lock
  if (!lock) {
    .udfEnv$fun <- list()
  }
  invisible(.udfEnv$lockedEnvir)
}
#' Lock the UDF function if the object exits inside of it
#'
#' @param obj object to check to see if it exists
#' @param envir When non-nil, look for object in environment and
#'   parent environments
#' @return logical saying if the environment was locked
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
.udfEnvLockIfExists <- function(obj, envir=NULL) {
  if (.udfEnv$lockedEnvir) return(invisible(FALSE))
  if (is.null(envir)) {
    if (any(vapply(ls(.udfEnv$envir, all=TRUE),
                   function(v) {
                     identical(obj, get(v, envir=.udfEnv$envir))
                   }, logical(1), USE.NAMES = FALSE))) {
      .udfEnvLock(lock=TRUE)
      return(invisible(TRUE))
    }
    return(invisible(FALSE))
  } else if (is.environment(envir)) {
    .env <- envir
    while(TRUE) {
      if (any(vapply(ls(.env, all=TRUE),
                     function(v) {
                       identical(obj, get(v, envir=.env))
                     }, logical(1), USE.NAMES = FALSE))) {
        .udfEnvSet(.env)
        .udfEnvLock(lock=TRUE)
        return(invisible(TRUE))
      }
      .env <- parent.env(.env)
      if (identical(.env, globalenv())) {
        if (any(vapply(ls(.env, all=TRUE),
                       function(v) {
                         identical(obj, get(v, envir=.env))
                       }, logical(1), USE.NAMES = FALSE))) {
          .udfEnvSet(.env)
          .udfEnvLock(lock=TRUE)
          return(invisible(TRUE))
        } else {
          return(invisible(FALSE))
        }
      }
      if (identical(.env, emptyenv())) break;
    }
  }
  invisible(FALSE)
}

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
  .fun <- try(get(fun, mode="function", envir=.udfEnv$envir), silent=TRUE)
  if (inherits(.fun, "try-error")) {
    .msg <- try(attr(.fun, "condition")$message, silent=TRUE)
    if (inherits(.msg, "try-error")){
      .msg <- sprintf("function '%s' is not supported; user function not found",
                      fun)
    }
    return(list(nargs=NA_integer_, .msg))
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
  .env <- new.env(parent=emptyenv())
  .env$needRecompile <- FALSE
  lapply(.n,
         function(n) {
           .oldArg <- iv[n]
           .new <- .getUdfInfo(n)
           if (any(names(.udfEnv$rxSEeqUsr) == n)) {
             .c <- .udfEnv$rxSEeqUsr[n]
             if (.c == .new[[1]]) {
               message("compiled with R user function '", n, "'; now there is a clashing C user function")
               .env$needRecompile <- TRUE
               message("triggered a recompile to use the C user function (they are always preferred)")
             } else {
               stop("there is both C and R user functions '", n, "' with a different number of arguments\n  since rxode2 prefers C, you will need to rename your R user function to use it")

             }
           }
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
  .env$needRecompile
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
#' @keywords internal
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
