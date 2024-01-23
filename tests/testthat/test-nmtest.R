system("rm -rv src/*.o src/*.so");devtools::load_all()

d <- nlmixr2data::nmtest
names(d) <- sub("lagt", "lagcentral",
                sub("bioav", "fcentral",
                    sub("rat2", "ratecentral",
                        sub("dur2", "durcentral", names(d)))))

p <- TRUE
library(ggplot2)

solveEqual <- function(id) {
  d2 <- d[d$id == id,]
  s1 <- linCmt(d2, cl=1.1, v=20, ka=1.5, sm=1000) |>
    dplyr::filter(EVID == 0) |>
    dplyr::mutate(Cc=Cc)
  assign("s1", s1, envir=globalenv())
  print(ggplot(data=s1, aes(TIME, Cc)) +
          geom_line(col="red", linewidth=1.2) +
          theme_bw() +
          geom_line(data=d2, aes(time, cp), col="blue", lty=2, linewidth=1.2) +
          ggtitle(paste0("id=", id, "; red=solved; blue: NONMEM")))
}

solveEqual(1)

solveEqual(2)

solveEqual(3)

solveEqual(4)

solveEqual(5)

solveEqual(6)

solveEqual(7)

solveEqual(8)

solveEqual(9)

solveEqual(13)

solveEqual(16)

solveEqual(17)

solveEqual(18)

# FIXME


solveEqual(10)

solveEqual(11)

solveEqual(12)

solveEqual(15)

solveEqual(19)

solveEqual(20)

solveEqual(21)

solveEqual(22)

solveEqual(23)

solveEqual(24)

solveEqual(25)

solveEqual(26)


## id <- unique(d$id)

## p <- FALSE
## lapply(id, function(i) {
##   meths <- c("liblsoda", "lsoda", "dop853")
##   modDat <- c("none", "rate", "dur")
##   for (meth in meths) {
##     for (modifyData in modDat) {
##       for (addlDropSs in c(TRUE, FALSE)) {
##         solveEqual(i, meth=meth, modifyData=modifyData, addlDropSs=addlDropSs)
##       }
##     }
##   }
## })

## invisible(lapply(unique(d$id), function(i){message(i);solveEqual(i);Sys.sleep(1)}))

## need to check id=25 type for infusions too (static, modeled) need
## to check for steady state when the infusion is still going at the
## time of steady-state release

## Need to check steady state infusion as well as infusions where ii=dur with lag times

## modeled equivalents of 425, 525
