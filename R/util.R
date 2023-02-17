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
  checkmate::assertNumeric(amt, any.missing=FALSE)
  checkmate::assertNumeric(dur, any.missing=FALSE)
  checkmate::assertNumeric(ii, any.missing=FALSE)
  checkmate::assertNumeric(ss, any.missing=FALSE)
  .df <- data.frame(cmt=as.integer(cmt), evid=as.integer(evid), amt=amt,
                    rate=rate, dur=dur, ii=ii, ss=ss)
  .Call(`_rxode2parse_getClassicEvid`,
        .df$cmt, .df$amt, .df$rate, .df$dur,
        .df$ii, .df$evid, .df$ss)
}
