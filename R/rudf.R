.udfEnv <- new.env(parent=emptyenv())
.udfEnv$fun <- list()

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
  return(list(nargs=.nargs,
              fun))
}

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
