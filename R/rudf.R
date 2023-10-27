.udfEnv <- new.env(parent=emptyenv())
.udfEnv$fun <- list()

.getUdfInfo <- function(fun) {
  .fun <- try(get(fun, mode="function"), silent=TRUE)
  if (inherits(.fun, "try-error")) {
    return(list(nargs=NA_integer_,
                sprintf("function '%s' is not supported; user not found",
                        fun)))
  }
  .formals <- formals(.fun)
  if (any(names(.formals) == "...")) {
    return(list(nargs=NA_integer_,
                "user defined R functions in rxode2 cannot have ... in part of the arguments"))
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
  .ret <- with(.envir, do.call(.fun, args))
  if (length(.ret) != 1L) return(NA_real_)
  .tmp <- try(as.double(.ret), silent=TRUE)
  if (inherits(.tmp, "try-error")) return(NA_real_)
  .ret
}
