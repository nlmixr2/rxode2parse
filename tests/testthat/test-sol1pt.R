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


test_that("single point 3-cmt linCmt oral bolus", {

  .v <- .solComp3(k10=0.1, k12=3, k21=1, k13=2, k31=0.5)

  pX <- c(10, 0, 0, 0)
  Xo <- rep(0, 3)
  dT <- 1
  ka <- 1.0
  E <- exp(-.v$L*dT) # Exponentials

  Xo <- Xo + pX[2] * .v$C1 %*% E
  Xo <- Xo + pX[3] * .v$C2 %*% E
  Xo <- Xo + pX[4] * .v$C3 %*% E

  Ea <- exp(-ka*dT)

  Xo <- Xo + ka*pX[1]*(.v$C1 %*% ((E - Ea)/(ka - .v$L)))

  Xo <- c(pX[1] * Ea, Xo)

  pX2 <- .solve1pt(pX, ka=ka, t=1, k10=0.1, k12=3, k21=1, k13=2, k31=0.5, rate=c(0, 0))

  expect_equal(Xo, pX2)

  # next time point
  pX <- pX2

  Xo <- rep(0, 3)
  dT <- 1
  ka <- 1.0
  E <- exp(-.v$L*dT) # Exponentials

  Xo <- Xo + pX[2] * .v$C1 %*% E
  Xo <- Xo + pX[3] * .v$C2 %*% E
  Xo <- Xo + pX[4] * .v$C3 %*% E

  Ea <- exp(-ka*dT)

  Xo <- Xo + ka*pX[1]*(.v$C1 %*% ((E - Ea)/(ka - .v$L)))

  Xo <- c(pX[1] * Ea, Xo)

  pX2 <- .solve1pt(pX, ka=ka, t=1, k10=0.1, k12=3, k21=1, k13=2, k31=0.5, rate=c(0, 0))

  expect_equal(Xo, pX2)

})

test_that("single point 3-cmt linCmt rate", {

  .v <- .solComp3(k10=0.1, k12=3, k21=1, k13=2, k31=0.5)

  pX <- c(0, 0, 0)
  rate <- 1
  Xo <- rep(0, 3)
  dT <- 1
  E <- exp(-.v$L*dT) # Exponentials

  Xo <- Xo + pX[1] * .v$C1 %*% E
  Xo <- Xo + pX[2] * .v$C2 %*% E
  Xo <- Xo + pX[3] * .v$C3 %*% E
  Xo <- Xo + ((rate*.v$C1) %*% ((1 - E)/.v$L)) # Infusion

  pX2 <- .solve1pt(pX, t=1, k10=0.1, k12=3, k21=1, k13=2, k31=0.5, rate=1.0)

  pX3 <- pX2
  dim(pX3) <- c(3L, 1L)

  expect_equal(Xo, pX3)

  # now next time
  pX <- pX2
  rate <- 1
  Xo <- rep(0, 3)
  dT <- 1
  E <- exp(-.v$L*dT) # Exponentials

  Xo <- Xo + pX[1] * .v$C1 %*% E
  Xo <- Xo + pX[2] * .v$C2 %*% E
  Xo <- Xo + pX[3] * .v$C3 %*% E
  Xo <- Xo + ((rate*.v$C1) %*% ((1 - E)/.v$L)) # Infusion

  pX2 <- .solve1pt(pX, t=1, k10=0.1, k12=3, k21=1, k13=2, k31=0.5, rate=1.0)

  pX3 <- pX2
  dim(pX3) <- c(3L, 1L)

  expect_equal(Xo, pX3)

})


test_that("single point 3-cmt linCmt rate+ka", {

  pX <- c(10, 0, 0, 0)
  Xo <- rep(0, 3)
  dT <- 1
  ka <- 1.0
  rate <- 1

  E <- exp(-.v$L*dT) # Exponentials

  Xo <- Xo + pX[2] * .v$C1 %*% E
  Xo <- Xo + pX[3] * .v$C2 %*% E
  Xo <- Xo + pX[4] * .v$C3 %*% E
  Xo <- Xo + ((rate*.v$C1) %*% ((1 - E)/.v$L)) # Infusion

  Ea <- exp(-ka*dT)

  Xo <- Xo + ka*pX[1]*(.v$C1 %*% ((E - Ea)/(ka - .v$L)))

  Xo <- c(pX[1] * Ea, Xo)

  pX2 <- .solve1pt(pX, ka=ka, t=1, k10=0.1, k12=3, k21=1, k13=2, k31=0.5, rate=c(0, rate))

  expect_equal(Xo, pX2)

  # next time point
  pX <- pX2

  Xo <- rep(0, 3)
  dT <- 1
  ka <- 1.0
  E <- exp(-.v$L*dT) # Exponentials

  Xo <- Xo + pX[2] * .v$C1 %*% E
  Xo <- Xo + pX[3] * .v$C2 %*% E
  Xo <- Xo + pX[4] * .v$C3 %*% E
  Xo <- Xo + ((rate*.v$C1) %*% ((1 - E)/.v$L)) # Infusion

  Ea <- exp(-ka*dT)

  Xo <- Xo + ka*pX[1]*(.v$C1 %*% ((E - Ea)/(ka - .v$L)))

  Xo <- c(pX[1] * Ea, Xo)

  pX2 <- .solve1pt(pX, ka=ka, t=1, k10=0.1, k12=3, k21=1, k13=2, k31=0.5, rate=c(0, rate))

  expect_equal(Xo, pX2)
})

test_that("2 compartment linCmt() bolus", {

  .v <- .solComp2(k10=0.1, k12=3, k21=1)

  pX <- c(10, 0)
  Xo <- rep(0, 2)
  dT <- 1
  E <- exp(-.v$L*dT) # Exponentials

  Xo <- Xo + pX[1] * .v$C1 %*% E
  Xo <- Xo + pX[2] * .v$C2 %*% E

  pX2 <- .solve1pt(pX, t=1, k10=0.1, k12=3, k21=1)

  pX3 <- pX2
  dim(pX3) <- c(2L, 1L)

  expect_equal(Xo, pX3)

  # solve next time too
  pX <- pX2
  Xo <- rep(0, 2)
  dT <- 1
  E <- exp(-.v$L*dT) # Exponentials

  Xo <- Xo + pX[1] * .v$C1 %*% E
  Xo <- Xo + pX[2] * .v$C2 %*% E

  pX2 <- .solve1pt(pX, t=1, k10=0.1, k12=3, k21=1)
  pX3 <- pX2
  dim(pX3) <- c(2L, 1L)

  expect_equal(Xo, pX3)
})

test_that("single point 2-cmt linCmt oral bolus", {

  .v <- .solComp2(k10=0.1, k12=3, k21=1)

  pX <- c(10, 0, 0)
  Xo <- rep(0, 2)
  dT <- 1
  ka <- 1.0
  E <- exp(-.v$L*dT) # Exponentials

  Xo <- Xo + pX[2] * .v$C1 %*% E
  Xo <- Xo + pX[3] * .v$C2 %*% E

  Ea <- exp(-ka*dT)

  Xo <- Xo + ka*pX[1]*(.v$C1 %*% ((E - Ea)/(ka - .v$L)))

  Xo <- c(pX[1] * Ea, Xo)

  pX2 <- .solve1pt(pX, ka=ka, t=1, k10=0.1, k12=3, k21=1, rate=c(0, 0))

  expect_equal(Xo, pX2)

  # next time point
  pX <- pX2

  Xo <- rep(0, 2)
  dT <- 1
  ka <- 1.0
  E <- exp(-.v$L*dT) # Exponentials

  Xo <- Xo + pX[2] * .v$C1 %*% E
  Xo <- Xo + pX[3] * .v$C2 %*% E

  Ea <- exp(-ka*dT)

  Xo <- Xo + ka*pX[1]*(.v$C1 %*% ((E - Ea)/(ka - .v$L)))

  Xo <- c(pX[1] * Ea, Xo)

  pX2 <- .solve1pt(pX, ka=ka, t=1, k10=0.1, k12=3, k21=1, k13=2, k31=0.5, rate=c(0, 0))

  expect_equal(Xo, pX2)

})
