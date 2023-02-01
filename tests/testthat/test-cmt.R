test_that("cmt()", {
  ## from draft translation of ddmore model DDMODEL00000302
  
  expect_error(rxode2parse("
cmt(ddta1)
cmt(ddta2)
cmt(ddta3)
cmt(ddta4)
d/dt(rxddta1) <- K21 * rxddta2 - K12 * rxddta1 - K10 * rxddta1
d/dt(rxddta2) <-  - K21 * rxddta2 + K12 * rxddta1
CPT <- rxddta1 / scale1
d/dt(rxddta5) <- KTR * CPT - KTR * rxddta5
d/dt(rxddta3) <- KTR * rxddta5 - KTR * rxddta3
A5 <- rxddta5
A3 <- rxddta3
EDRUGT <- ALPHA * rxddta3
HAZT <- 0
if (t > 0) HAZT <- BETA * (EDRUGT^BETA) * (t^(BETA - 1)) * COVAR
d/dt(rxddta4) <- HAZT
d/dt(rxddta6) <- CPT
AUCT <- rxddta6
CAV <- AUCT / t
CP <- rxddta1 / scale1
EDRUG <- ALPHA * rxddta3
"), NA)
  
})
