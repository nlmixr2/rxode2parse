#' Get the internal breakdown of an evid
#'
#' @param i evid to breakdown
#' @return named evid integer vector
#' @export 
#' @author Matthew L. Fidler
#' @keywords internal
#' @examples
#' 
#' .getWh(1001)
#' .getWh(10401)
#' 
.getWh <- function(i) {
  checkmate::assertIntegerish(i,len=1, any.missing=FALSE)
  .Call(`_rxode2parse_getWh`, as.integer(i))
}

#' This converts NONMEM-style EVIDs to classic RxODE events
#'  
#' @param cmt compartment flag
#' @param amt dose amount
#' @param rate dose rate
#' @param dur dose duration
#' @param ii inter-dose interval
#' @param evid event id
#' @param ss steady state
#' @return classic evids, excluding evids that are added (you need to
#'   add them manually) or simply use etTran.  This is mostly for
#'   testing and really shouldn't be used directly.
#' @export
#' @author Matthew L. Fidler
#' @examples
#' .toClassicEvid(cmt=10, amt=3, evid=1)
#' .toClassicEvid(cmt=10, amt=3, rate=2, evid=1)
#' .toClassicEvid(cmt=10, amt=3, rate=-1, evid=1)
#' .toClassicEvid(cmt=10, amt=3, rate=-2, evid=1)
#' .toClassicEvid(cmt=10, amt=3, dur=2, evid=1)
#' .toClassicEvid(cmt=304, amt=3, dur=2, evid=1)
#' .toClassicEvid(cmt=7, amt=0, rate=2, evid=1, ss=1)
#' .toClassicEvid(cmt=-10, amt=3, evid=1)
#' .toClassicEvid(cmt=10, amt=3, evid=5)
#' .toClassicEvid(cmt=6, amt=3, evid=6)
#' .toClassicEvid(cmt=6, amt=3, evid=7)
#' .toClassicEvid(evid=2)
#' .toClassicEvid(evid=4)
.toClassicEvid <- function(cmt=1L, amt=0.0, rate=0.0, dur=0.0, ii=0.0, evid=0L, ss=0.0) {
  checkmate::assertIntegerish(cmt, any.missing=FALSE)
  checkmate::assertIntegerish(evid, any.missing=FALSE)
  checkmate::assertNumeric(amt)
  checkmate::assertNumeric(dur, any.missing=FALSE)
  checkmate::assertNumeric(ii)
  checkmate::assertNumeric(ss)
  .df <- data.frame(cmt=as.integer(cmt), evid=as.integer(evid), amt=as.double(amt),
                    rate=as.double(rate), dur=as.double(dur),
                    ii=as.double(ii),
                    ss=as.double(ss))
  .Call(`_rxode2parse_getClassicEvid`,
        .df$cmt, .df$amt, .df$rate, .df$dur,
        .df$ii, .df$evid, .df$ss)
}
#' This converts data to a format similar to what rxode2 will use
#'  
#' @param data input dataset
#' @return linCmt dataset
#' @noRd
#' @author Matthew L. Fidler
.toLinCmtData <- function(data) {
  .lowerNames <- tolower(names(data))
  if ("id" %in% .lowerNames) {
    stop("cannot solve multi-subject data, no 'id'", call.=FALSE)
  }
  if (!("evid" %in% .lowerNames)) {
    stop("need 'evid'", call.=FALSE)
  }
  if (!("amt" %in% .lowerNames)) {
    stop("need 'amt'", call.=FALSE)
  }
  if (!("time" %in% .lowerNames)) {
    stop("need 'time'", call.=FALSE)
  }
  .dat <- data
  names(.dat) <- .lowerNames
  if (!("cmt" %in% .lowerNames)) {
    .dat$cmt <- 1
  }
  .dat <- .dat[order(.dat$time),]
  for (.v in c("rate", "dur", "ii", "ss")) {
    if (!(.v %in% .lowerNames)) {
      .dat[[.v]] <- 0.0
    }
  }
  
  .evid <- .toClassicEvid(cmt=.dat$cmt, amt=.dat$amt,
                          rate=.dat$rate, dur=.dat$dur,
                          ii=.dat$ii, evid=.dat$evid,
                          ss=.dat$ss)
  .dat$nmEvid <- .dat$evid
  .dat$evid <- .evid
  .time <- .dat$time
  .evid <- .dat$evid
  .dat0 <- .dat[.dat$evid != 0,]
  .datL <- do.call("rbind",
                   lapply(seq_along(.dat0$evid),
                          function(i) {
                            .wh <- .getWh(.dat0$evid[i])
                            .tinf <- 0.0
                            if (.dat0$rate[i] > 0) {
                              # rate = amt/time; tinf= amt/rate
                              .tinf <- .dat0$amt[i]/.dat0$rate[i]
                            } else if (.dat0$dur[i] > 0) {
                              .tinf <- .dat0$dur[i]
                            }
                            .len <- length(.tinf)
                            data.frame(time=.dat0$time[i],
                                       dose=.dat0$amt[i],
                                       tinf=.tinf,
                                       ii=.dat0$ii[i],
                                       cmt=.dat0$cmt,
                                       evidF=as.integer(.wh["whI"]),
                                       evid0=as.integer(.wh["wh0"]),
                                       # these are saved/modified
                                       # because of time-varying
                                       # possibilities
                                       tlag=rep(0.0, .len),
                                       f=rep(1.0, .len),
                                       rate=rep(0.0, .len))
                          }))
  row.names(.datL) <- NULL
  list(time=as.double(.dat$time),
       evid=as.integer(.dat$evid),
       linDat=.datL,
       len=length(.datL$time))
}

.toLinCmtParam <- function(cmt, trans,
                           p1, v1, p2, p3, p4, p5,
                           tlag, F, rate1, dur1,
                           ka, tlag2, F2, rate2, dur2) {
  return(list(c(cmt=as.integer(cmt), trans=as.integer(trans)),
              c(p1=p1, v1=v1, p2=p2, p3=p3, p4=p4, p5=p5,
                tlag=tlag, F=F, rate1=rate1, dur1=dur1,
                ka=ka, tlag2=tlag2, F2=F2, rate2=rate2, dur2=dur2)))
}



rex::register_shortcuts("rxode2parse")
.rxDerivedReg <- rex::rex(
  start,
  or(
    group(or("V", "Q", "VP", "VT", "CLD"), number),
    "KA", "VP", "VT", "CLD", "V", "VC", "CL", "VSS", "K", "KE", "KEL",
    "Q", "VT", group("K", number, number), "AOB", "ALPHA", "BETA", "GAMMA",
    "A", "B", "C"
  ),
  end
)


#' Calculate derived parameters for the 1-, 2-, and 3- compartment
#' linear models.
#'
#' This calculates the derived parameters based on what is provided
#' in a data frame or arguments
#'
#' @param ... The input can be:
#'
#'
#'  * A data frame with PK parameters in it; This should ideally
#'  be a data frame with one pk parameter per row since it will
#'  output a data frame with one PK parameter per row.
#'
#'  * PK parameters as either a vector or a scalar
#'
#'
#' @param verbose boolean that when TRUE provides a message about the detected pk parameters
#'   and the detected compartmental model.  By default this is `FALSE`.
#'
#' @param digits represents the number of significant digits for the
#'   output; If the number is zero or below (default), do not round.
#'
#' @return Return a data.frame of derived PK parameters for a 1-, 2-,
#'   or 3-compartment linear model given provided clearances and
#'   volumes based on the inferred model type.
#'
#' The model parameters that will be provided in the data frame are:
#'
#' * `vc`: Central Volume (for 1-, 2- and 3-
#'   compartment models)
#'
#' * `kel`: First-order elimination rate (for 1-, 2-, and
#'   3-compartment models)
#'
#' * `k12`: First-order rate of transfer from central to
#'   first peripheral compartment; (for 2- and 3-compartment models)
#'
#' * `k21`: First-order rate of transfer from first
#'   peripheral to central compartment, (for 2- and 3-compartment
#'   models)
#'
#' * `k13`: First-order rate of transfer from central to
#'   second peripheral compartment; (3-compartment model)
#'
#' * `k31`: First-order rate of transfer from second
#'   peripheral to central compartment (3-compartment model)
#'
#' * `vp`: Peripheral Volume (for 2- and 3- compartment models)
#'
#' * `vp2`: Peripheral Volume for 3rd compartment (3- compartment model)
#'
#' * `vss`: Volume of distribution at steady state; (1-, 2-, and 3-compartment models)
#'
#' * `t12alpha`: \eqn{t_{1/2,\alpha}}; (1-, 2-, and 3-compartment models)
#'
#' * `t12beta`: \eqn{t_{1/2,\beta}}; (2- and 3-compartment models)
#'
#' * `t12gamma`: \eqn{t_{1/2,\gamma}}; (3-compartment model)
#'
#' * `alpha`: \eqn{\alpha}; (1-, 2-, and 3-compartment models)
#'
#' * `beta`: \eqn{\beta}; (2- and 3-compartment models)
#'
#' * `gamma`: \eqn{\beta}; (3-compartment model)
#'
#' * `A`: true `A`; (1-, 2-, and 3-compartment models)
#'
#' * `B`: true `B`; (2- and 3-compartment models)
#'
#' * `C`: true `C`; (3-compartment model)
#'
#' * `fracA`: fractional A; (1-, 2-, and 3-compartment models)
#'
#' * `fracB`: fractional B; (2- and 3-compartment models)
#'
#' * `fracC`: fractional C; (3-compartment model)
#'
#' @author Matthew Fidler and documentation from Justin Wilkins, \email{justin.wilkins@@occams.com}
#'
#' @references Shafer S. L. `CONVERT.XLS`
#'
#' @references Rowland M, Tozer TN. Clinical Pharmacokinetics and Pharmacodynamics: Concepts and Applications (4th). Clipping Williams & Wilkins, Philadelphia, 2010.
#'
#' @examples
#'
#' ## Note that rxode2 parses the names to figure out the best PK parameter
#'
#' params <- rxDerived(cl = 29.4, v = 23.4, Vp = 114, vp2 = 4614, q = 270, q2 = 73)
#'
#' ## That is why this gives the same results as the value before
#'
#' params <- rxDerived(CL = 29.4, V1 = 23.4, V2 = 114, V3 = 4614, Q2 = 270, Q3 = 73)
#'
#' ## You may also use micro-constants alpha/beta etc.
#'
#' params <- rxDerived(k12 = 0.1, k21 = 0.2, k13 = 0.3, k31 = 0.4, kel = 10, v = 10)
#'
#' ## or you can mix vectors and scalars
#'
#' params <- rxDerived(CL = 29.4, V = 1:3)
#'
#' ## If you want, you can round to a number of significant digits
#' ## with the `digits` argument:
#'
#' params <- rxDerived(CL = 29.4, V = 1:3, digits = 2)
#' @export
rxDerived <- function(..., verbose = FALSE, digits = 0) {
  .lst <- list(...)
  if (inherits(.lst[[1]], "data.frame")) {
    .lst <- .lst[[1]]
  }
  .namesU <- toupper(names(.lst))
  .w <- which(regexpr(.rxDerivedReg, .namesU) != -1)
  if (length(.w) > 1L) {
    if (verbose) {
      message("parameters: ", paste(names(.lst)[.w], collapse = ","))
    }
    .linCmt <- .Call(
      `_rxode2parse_linCmtParse`, names(.lst)[.w],
      c(
        "with(.lst,.Call(`_rxode2parse_calcDerived`, ", "list(", "0, 0, 0, 0, ",
        ", 0, 0, 0, 0),digits))"
      ),
      verbose
    )$str
    .env <- environment()
    return(eval(parse(text = .linCmt), envir = .env))
  } else {
    stop("cannot figure out PK parameters to convert", call. = FALSE)
  }
}

#' Linear compartmental concentrations based on input data  
#'  
#' @param data This dataest should have the following columns:
#' 
#'   - `evid` - `rxode2` style evid
#' 
#'   - `amt` - amount of dose
#' 
#'   - `time` - time of observation
#'
#' The following columns are optional:
#'
#' - `id` the identifier of the individual
#'
#'  - `cmt` the compartment of dose
#' 
#'  - `rate` the rate of dose
#'
#'  - `dur` the duration of the dose
#'
#'  - `ii` the inter-dose interval
#'
#'  - `ss` the steady state flag
#'
#' This cannot have the an `id` column
#' 
#' @param ... Parameters for linear compartment model
#' 
#' @inheritParams rxDerived
#' 
#' @return original dataset `data` with an additional column `Cp` with
#'   the
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' rxLinCmt(data=nlmixr2data::theo_sd, CL = 29.4, V = 3)
#' 
rxLinCmt <- function(data, ..., tlag=0, F=1, rate1=0, dur1=0,
                     tlag2=0, F2=1, rate2=0, dur2=0,
                     verbose = FALSE) {
  .lowNames <- tolower(names(data))
  .w <- which(.lowNames == "id")
  .wt <- which(.lowNames == "time")
  if (length(.wt) != 1L) stop("need one time column", call.=FALSE)
  if (length(.w) == 0L) {
    .ids <- 1
    .data <- data[order(data[, .wt]),]
    .datLin <- list(.toLinCmtData(.data))
  } else if (length(.w) == 1L) {
    .data <- data[order(data[, .w], data[, .wt]),]
    .ids <- unique(data[,.w])
    .datLin <- lapply(seq_along(.ids),
                      function(id){
                        .d <- .data[data[,.w] == id, ]
                        .d <- .d[, -.w]
                        .toLinCmtData(.d)
                      })
  } else {
    stop("can't determine 'id' column, there seems to be duplicates",
         call.=FALSE)
  }
  .lst <- list(...)
  .namesU <- toupper(names(.lst))
  .w <- which(regexpr(.rxDerivedReg, .namesU) != -1)
  if (length(.w) > 1L) {
    if (verbose) {
      message("parameters: ", paste(names(.lst)[.w], collapse = ","))
    }
    .linCmt <- .Call(
      `_rxode2parse_linCmtParse`, names(.lst)[.w],
      c(
        "with(.lst, .toLinCmtParam(", "", "tlag, F, rate1, dur1, ",
        ", tlag2, F2, rate2, dur2))"
      ),
      verbose
    )
    .param <- eval(parse(text=.linCmt$str))
    .res <- lapply(seq_along(.datLin),
                   function(i) {
                     .Call(`_rxode2parse_linCmtA`, .datLin[[i]],
                           .param)
                   })
    .data <- data.frame(.data, Cp=unlist(.res))
    return(.data)
  } else {
    stop("cannot figure out PK parameters to use for PK curve", call. = FALSE)
  }
}
