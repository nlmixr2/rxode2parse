devtools::load_all()

d <- nlmixr2data::nmtest
names(d) <- sub("lagt", "lagcentral",
                sub("bioav", "fdepot",
                    sub("rat2", "ratecentral",
                        sub("dur2", "durcentral", names(d)))))

p <- TRUE
library(ggplot2)

solveEqual <- function(id) {
  d <- d[d$id == id,]
  s1 <- linCmt(d, cl=1.1, v=20/1000, ka=1.5) |>
    dplyr::filter(EVID == 0)
  print(ggplot(data=s1, aes(TIME, Cc)) +
          geom_line(col="red") +
          theme_bw() +
          geom_line(data=d, aes(time, cp), col="blue") )
}

solveEqual(1)

id <- unique(d$id)

p <- FALSE
lapply(id, function(i) {
  meths <- c("liblsoda", "lsoda", "dop853")
  modDat <- c("none", "rate", "dur")
  for (meth in meths) {
    for (modifyData in modDat) {
      for (addlDropSs in c(TRUE, FALSE)) {
        solveEqual(i, meth=meth, modifyData=modifyData, addlDropSs=addlDropSs)
      }
    }
  }
})

## invisible(lapply(unique(d$id), function(i){message(i);solveEqual(i);Sys.sleep(1)}))

## need to check id=25 type for infusions too (static, modeled) need
## to check for steady state when the infusion is still going at the
## time of steady-state release

## Need to check steady state infusion as well as infusions where ii=dur with lag times

## modeled equivalents of 425, 525
