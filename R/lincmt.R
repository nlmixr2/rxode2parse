.EVIDF_INF_RATE <- 1L
.EVIDF_INF_DUR <- 2
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
                            .dose <- .dat0$amt[i]
                            if (.wh["whI"] == .EVIDF_INF_RATE ||
                                  .wh["whI"] == .EVIDF_INF_DUR) {
                              .dose <- .dat0$rate[i]
                            }
                            data.frame(time=.dat0$time[i],
                                       dose=.dose,
                                       tinf=.tinf,
                                       ii=.dat0$ii[i],
                                       cmt=.dat0$cmt[i]-1L,
                                       evidF=as.integer(.wh["whI"]),
                                       evid0=as.integer(.wh["wh0"]))
                          }))
  .len <- length(.datL$tinf)
  # these are saved/modified because of time-varying possibilities
  .datL <-  data.frame(.datL, tlag=rep(0.0, .len),
                       f=rep(1.0, .len),
                       rate=rep(0.0, .len))
  row.names(.datL) <- NULL
  .ret <-  list(time=as.double(.dat$time),
                evid=as.integer(.dat$evid),
                linDat=.datL,
                len=length(.datL$time))
  .ret
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
#' @param ... Parameters for linear compartment model, parsed automatically
#'
#' @param tlag This is the lag time of the first compartment.
#' @param F This is the bioavailability of the first compartment
#' @param rate1 This is the modeled rate for the first compartment
#' @param dur1 This is the modeled duration for the first compartment
#' @param tlag2 This is the lag time of the second compartment
#' @param F2 This is the modeled bioavailability of the second
#'   compartment
#' @param rate2 This is the modeled rate of the second compartment
#' @param dur2 This is the modeled duration of the second compartment
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
