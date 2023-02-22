tol <- 1e-5
tolSS <- 5e-5

.qr <- function(x) {
  qs::qread(testthat::test_path(x))
}

.qs <- function(x, y) {
  qs::qsave(x,testthat::test_path(y))
}

etSsB <- .qr("test-lincmt-etSsB.qs")
etSsI <- .qr("test-lincmt-etSsI.qs")
etSsR <- .qr("test-lincmt-etSsR.qs")

.txt <- "non sens"

etSsB.ode.1c <- .qr("test-lincmt-etSsB.ode.1c.qs")

d <- etSsB.ode.1c[,c("evid", "cmt", "amt", "ii", "ss", "time")]

s1 <- rxLinCmt(d, V=20, CL=25)

test_that(sprintf("one compartment bolus steady state (%s)", .txt), {
  expect_equal(etSsB.ode.1c$C2, s1$Cp, tolerance=tol)
})

etSsI.ode.1c <- .qr("test-lincmt-etSsI.ode.1c.qs")

d <- etSsI.ode.1c[,c("evid", "cmt", "amt", "ii", "ss", "rate", "time")]

s1 <- rxLinCmt(d, V=20, CL=25)

test_that(sprintf("one compartment infusion tau steady state (%s)", .txt), {
  expect_equal(etSsI.ode.1c$C2, s1$Cp, tolerance = tolSS)
})


etSsR.ode.1c <- .qr("test-lincmt-etSsR.ode.1c.qs")

d <- etSsR.ode.1c[,c("evid", "cmt", "amt", "ii", "ss", "rate", "time")]

s1 <- rxLinCmt(d, V=20, CL=25)

test_that(sprintf("one compartment infusion tau steady state (%s)", .txt), {
  expect_equal(etSsR.ode.1c$C2, s1$Cp, tolerance = tolSS)
})


## library(ggplot2)

## ggplot(s1, aes(time, Cp)) + geom_line() +
##   geom_line(data=etSsR.ode.1c, aes(time, C2), col="red")
