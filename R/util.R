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
                            data.frame(time=.dat0$time[i],
                                       dose=.dat0$amt[i],
                                       tinf=.tinf,
                                       ii=.dat0$ii[i],
                                       evidF=.wh["whI"],
                                       evid0=.wh["wh0"])
                          }))
  row.names(.datL) <- NULL
  list(time=.dat$time, evid=.dat$evid, linDat=.datL,
       len=length(.datL$time))
}
