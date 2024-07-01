system("rm -rv src/*.o src/*.so");devtools::load_all()

d <- nlmixr2data::nmtest
names(d) <- sub("lagt", "lagcentral",
                sub("bioav", "fcentral",
                    sub("rat2", "ratecentral",
                        sub("dur2", "durcentral", names(d)))))

p <- TRUE
library(ggplot2)

solveEqual <- function(id, modifyData = c("none", "dur", "rate"),
                       addlKeepsCov = TRUE, addlDropSs=TRUE, ss2cancelAllPending=FALSE) {
  d2 <- d[d$id == id,]
  s1 <- linCmt(d2, cl=1.1, v=20, ka=1.5, sm=1000) |>
    dplyr::filter(EVID == 0) |>
    dplyr::mutate(Cc=Cc)
  assign("s1", s1, envir=globalenv())
  ggplot(data=s1, aes(TIME, Cc)) +
          geom_line(col="red", linewidth=1.2) +
          theme_bw() +
          geom_line(data=d2, aes(time, cp), col="blue", lty=2, linewidth=1.2) +
          ggtitle(paste0("id=", id, "; red=solved; blue: NONMEM"))
}

solveEqual <- function(id, plot = p, modifyData = c("none", "dur", "rate"),
                       addlKeepsCov = TRUE, addlDropSs=TRUE, ss2cancelAllPending=FALSE, tol=0.01) {
  hasRate <- any(d[d$id == id & d$evid != 0,]$rate != 0)
  hasModeledRate <- any(d[d$id == id & d$evid != 0,]$mode == 1)
  hasModeledDur  <- any(d[d$id == id & d$evid != 0, ]$mode == 2)
  hasChangedF <- any(d[d$id == id & d$evid != 0, ]$bioav != 1)
  modifyData <- match.arg(modifyData)
  d <- d[d$id == id,]
  rate <- unlist(as.vector(d[d$evid != 0, "rate"]))
  ii0 <- all(d$ii == 0)
  oneRate <- (length(rate) == 1L)
  dose1 <- all(d[d$evid != 0, ]$cmt == 1)
  if (!addlDropSs) {
    if (any(d$ss == 2)) {
      return()
    }
  }
  if (modifyData == "rate" && hasRate && !hasModeledRate && !hasModeledDur && oneRate && !ii0 && !dose1) {
    if (p) message("modified rate to be modeled")
    rate <- d[d$evid != 0, "rate", drop=FALSE]
    rate <- rate$rate
    if (length(rate) == 1) {
      d$ratecentral <- rate
      d$rate <- ifelse(d$rate==0, 0, -1)
      d$mode <- 1
    }
    ## print(etTrans(d, f, addlDropSs=addlDropSs))
  } else if (modifyData == "dur" && hasRate && !hasModeledRate && !hasModeledDur && !hasChangedF && oneRate && !ii0 && !dose1) {
    if (p) message("modified dur to be modeled")
    rate <- as.numeric(d[d$evid != 0, "rate"])
    amt <- as.numeric(d[d$evid != 0, "amt"])
    if (length(rate) == 1) {
      d$durcentral <- amt/rate
      d$ratecentral <- ifelse(d$rate==0, 0, -2)
      d$mode <- 2
      assign(".d", d, envir=globalenv())
    }
  } else if (any(modifyData == c("dur", "rate"))) {
    if (p) {
      message("skipping because cannot be modified")
      print(list(modifyData=modifyData,
                 hasRate=hasRate,
                 hasModeledRate=hasModeledRate,
                 hasModeledDur= hasModeledDur,
                 hasChangedF=hasChangedF))
    }
    return(invisible())
  }
  ## if (plot) {
  ## print(etTrans(d, fl))
  if (p) {
    print(d[d$evid != 0, ])
    s1 <- linCmt(d, cl=1.1, v=20, ka=1.5, sm=1000,
                 addlKeepsCov = addlKeepsCov, addlDropSs=addlDropSs,
                 ss2cancelAllPending=ss2cancelAllPending) |>
      dplyr::filter(EVID == 0)
    print(ggplot(data=s1, aes(TIME, Cc)) +
            geom_line(col="blue") +
            geom_point(data=d, aes(x=time, y=cp), col="red") +
            theme_bw() +
            ggtitle(paste0("id=", id, " NONMEM: blue; rxode2: red; modify:",  modifyData, " addlDropSs: ", addlDropSs)))
  }
  ## } else {
    sub <- 0
    if (id %in% c(410, 411, 409, 415, 709, 510, 610)) {
      sub <- 24
    }
    if (id %in% c(809, 909, 1009)) {
      sub <- 48
    }
    if (id %in% 825) {
      sub <- 96
    }
    test_that(paste0("nmtest id:", id, " alag; modifyData:", modifyData,"; addlDropSs: ", addlDropSs),
    {
      s1 <- linCmt(d, cl=1.1, v=20, ka=1.5, sm=1000,
                   addlKeepsCov = addlKeepsCov, addlDropSs=addlDropSs,
                   ss2cancelAllPending=ss2cancelAllPending) |>
        dplyr::filter(EVID == 0)
      expect_equal(s1$Cc[s1$TIME >= sub],
                   d[d$id == id & d$evid == 0 & d$time >= sub,]$cp)
    })
  ## }
}

id <- unique(d$id)

#p <- FALSE
env <- new.env(parent=emptyenv())
env$err <- data.frame(i=integer(0), modifyData=character(0),
                      addlDropSs=logical(0))

lapply(id, function(i) {
  modDat <- c("none", "rate", "dur")
  for (modifyData in modDat) {
    for (addlDropSs in c(TRUE, FALSE)) {
      message("modDat: ", modifyData, " addlDropSs: ", addlDropSs, " i: ", i)
      t <- try(solveEqual(i, modifyData=modifyData, addlDropSs=addlDropSs))
      if (inherits(t, "try-error")) {
        env$err <- rbind(env$err,
                         data.frame(i=i, modifyData=modifyData,
                                    addlDropSs=addlDropSs))
      }
    }
  }
})


# FIXME

solveEqual(10)

solveEqual(12)

solveEqual(102)

solveEqual(109)

solveEqual(110)

solveEqual(210)

solveEqual(310)

solveEqual(410)

solveEqual(510)

solveEqual(610)

#?
solveEqual(415)

# Fixed

solveEqual(1)

solveEqual(2)

solveEqual(3)

solveEqual(4)

solveEqual(5)

solveEqual(6)

solveEqual(7)

solveEqual(8)

solveEqual(9)

solveEqual(11)

solveEqual(13)

solveEqual(15)

solveEqual(16)

solveEqual(17)

solveEqual(18)

solveEqual(19)

solveEqual(20)

solveEqual(21)

solveEqual(22)

solveEqual(23)

solveEqual(24)

solveEqual(25)

solveEqual(26)

solveEqual(111)

solveEqual(125)

solveEqual(211)

solveEqual(311)

solveEqual(325)

solveEqual(411)

solveEqual(425)

solveEqual(509)

solveEqual(525)


solveEqual(609)

solveEqual(625)

solveEqual(709)

solveEqual(725)

solveEqual(809)

solveEqual(909)

solveEqual(1009)

solveEqual(225)

# ?
solveEqual(409)

solveEqual(825)


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
