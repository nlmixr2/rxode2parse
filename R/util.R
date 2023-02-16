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
