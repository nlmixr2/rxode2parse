test_that("single point 3-cmt linCmt bolus", {

  .v <- .solComp3(k10=0.1, k12=3, k21=1, k13=2, k31=0.5)

  pX <- c(10, 0, 0)
  Xo <- rep(0, 3)
  dT <- 1
  E <- exp(-.v$L*dT) # Exponentials

  Xo <- Xo + pX[1] * .v$C1 %*% E
  Xo <- Xo + pX[2] * .v$C2 %*% E
  Xo <- Xo + pX[3] * .v$C3 %*% E

  pX2 <- .solve1pt(pX, t=1, k10=0.1, k12=3, k21=1, k13=2, k31=0.5)

  pX3 <- pX2
  dim(pX3) <- c(3L, 1L)

  expect_equal(Xo, pX3)

  # solve next time too
  pX <- pX2
  Xo <- rep(0, 3)
  dT <- 1
  E <- exp(-.v$L*dT) # Exponentials

  Xo <- Xo + pX[1] * .v$C1 %*% E
  Xo <- Xo + pX[2] * .v$C2 %*% E
  Xo <- Xo + pX[3] * .v$C3 %*% E

  pX2 <- .solve1pt(pX, t=1, k10=0.1, k12=3, k21=1, k13=2, k31=0.5)
  pX3 <- pX2
  dim(pX3) <- c(3L, 1L)


  expect_equal(Xo, pX3)

})
