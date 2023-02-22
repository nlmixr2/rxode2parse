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

etSsB.ode.2c <- .qr("test-lincmt-etSsB.ode.2c.qs")
d <- etSsB.ode.2c[,c("evid", "cmt", "amt", "ii", "ss", "time")]

s2 <- rxLinCmt(d, V = 40, CL = 18, V2 = 297, Q = 10)

test_that(sprintf("two compartment bolus steady state (%s)", .txt), {
  expect_equal(etSsB.ode.2c$C2, s2$Cp, tolerance=tol)
})


etSsI.ode.2c <- .qr("test-lincmt-etSsI.ode.2c.qs")
d <- etSsI.ode.2c[,c("evid", "cmt", "amt", "ii", "ss", "rate", "time")]

s2 <- rxLinCmt(d, V = 40, CL = 18, V2 = 297, Q = 10)

test_that(sprintf("two compartment infusion tau steady state (%s)", .txt), {
  expect_equal(etSsI.ode.2c$C2, s2$Cp, tolerance=tol)
})

etSsR.ode.2c <- .qr("test-lincmt-etSsR.ode.2c.qs")
d <- etSsR.ode.2c[,c("evid", "cmt", "amt", "ii", "ss", "rate", "time")]

s2 <- rxLinCmt(d, V = 40, CL = 18, V2 = 297, Q = 10)

test_that(sprintf("two compartment infusion steady state (%s)", .txt), {
  expect_equal(etSsR.ode.2c$C2, s2$Cp, tolerance=tol)
})
