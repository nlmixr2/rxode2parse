tolerance <- 1e-5
.qr <- function(x) {
  qs::qread(testthat::test_path(x))
}

etSsB <- .qr("test-lincmt-etSsB.qs")
etSsI <- .qr("test-lincmt-etSsI.qs")
etSsR <- .qr("test-lincmt-etSsR.qs")
tolerance <- 1e-5

.txt <- "non sens"

etSsB.ode.1c <- .qr("test-lincmt-etSsB.ode.1c.qs")

d <- etSsB.ode.1c[,c("evid", "cmt", "amt", "ii", "ss", "time")]

s1 <- rxLinCmt(d, V=20, CL=25)

test_that(sprintf("one compartment bolus steady state (%s)", .txt), {
  expect_equal(etSsB.ode.1c$C2, s1$Cp, tolerance=tolerance)
})

