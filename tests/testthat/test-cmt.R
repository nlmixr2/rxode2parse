test_that("cmt()", {
  ## from draft translation of ddmore model DDMODEL00000302
  
  expect_error(rxode2parse("
cmt(ddta1)
cmt(ddta2)
cmt(ddta3)
cmt(ddta4)
RTTE <- 0
AGEREF <- 65
ALBREF <- 4.0
BUNREF <- 16
BODYWEIGHTREF <- 85
if (SEXX == 1) SEX <- 1
if (SEXX == 2) {
SEX <- 0
BODYWEIGHTREF <- 68
}
BWT <- BW
if (BW ==  - 99) BWT <- BODYWEIGHTREF
THETA1 <- 0.312
THETA2 <- 1.21
THETA3 <- 0.957
THETA4 <-  - 1.02
THETA5 <- 1.48
THETA6 <- 1.02
THETA7 <- 0.476
THETA8 <- 0.527
THETA9 <- 0.484
THETA10 <- 0.303
THETA11 <- 0.149
THETA12 <- 0.223
THETA13 <-  - 0.212
LWT75 <- log(BWT / 75)
MUX1 <- THETA1 + THETA7 * LWT75 + THETA13 * log(NDOS / 2.4)
MUX2 <- THETA2 + THETA8 * LWT75 + THETA11 * SEX
MUX3 <- THETA3 + THETA9 * LWT75 + THETA12 * SEX
MUX4 <- THETA4 + THETA10 * LWT75
MUX5 <- THETA5
MUX6 <- THETA6
CLINF <- exp(MUX1 + ET1)
V1 <- exp(MUX2 + ET2)
V2 <- exp(MUX3 + ET3)
Q <- exp(MUX4 + ET4)
KDES <- exp(MUX5 + ET5)
CLT <- exp(MUX6 + ET6)
scale1 <- V1 / 1000
K12 <- Q / V1
K21 <- Q / V2
LOGKTR <- theta1 + eta1
KTR <- exp(LOGKTR)
ALPHA <- exp(theta2)
BETA <- exp(theta3)
COVAR <- theta4 * (BWT - 75) + theta5 * PPN
CL <- CLT * exp( - KDES * t) + CLINF
K10 <- CL / V1
d/dt(rxddta1) <- K21 * rxddta2 - K12 * rxddta1 - K10 * rxddta1
d/dt(rxddta2) <-  - K21 * rxddta2 + K12 * rxddta1
CPT <- rxddta1 / scale1
d/dt(rxddta5) <- KTR * CPT - KTR * rxddta5
d/dt(rxddta3) <- KTR * rxddta5 - KTR * rxddta3
A5 <- rxddta5
A3 <- rxddta3
EDRUGT <- ALPHA * rxddta3
HAZT <- 0
if (t > 0) HAZT <- BETA * (EDRUGT^BETA) * (t^(BETA - 1)) * exp(COVAR)
d/dt(rxddta4) <- HAZT
d/dt(rxddta6) <- CPT
AUCT <- rxddta6
CAV <- AUCT / t
CP <- rxddta1 / scale1
EDRUG <- ALPHA * rxddta3
HAZ <- 0
if (t > 0) HAZ <- BETA * (EDRUG^BETA) * (t^(BETA - 1)) * exp(COVAR)
SURV <- exp( - rxddta4)
PDF <- SURV * HAZ
if (DV == 0) {
Y <- SURV
CHLAST <- rxddta4
} else {
CHLAST <- CHLAST
}
if (DV == 1) {
Y <- PDF
}
if (icall == 4) {
if (newind != 2) R = rxunif()
DV <- 0
RTTE <- 0
if (CENS == 1) RTTE <- 1
if (R > SURV) {
DV <- 1
RTTE <- 1
}
Y <- DV
}
"), NA)
  
})
