digits <- 4
etSsB <- qs::qread("test-lincmt-etSsB.qs")
etSsI <- qs::qread("test-lincmt-etSsI.qs")
etSsR <- qs::qread("test-lincmt-etSsR.qs")
tolerance <- 1e-5

.txt <- "non sens"

etSsB.ode.1c <- qs::qread("test-lincmt-etSsB.ode.1c.qs")

d <- etSsB.ode.1c[,c("evid", "cmt", "amt", "ii", "ss", "time")]

s1 <- rxLinCmt(d, V=20, CL=25)

test_that(sprintf("one compartment bolus steady state (%s)", .txt), {
  expect_equal(etSsB.ode.1c$C2, s1$Cp, tolerance=1e-5)
})

