.rxLinInfoKaReg <- rex::rex(
  start,
  or("fdepot","fcentral",
     "ratedepot","ratecentral",
     "durdepot","durcentral",
     "lagdepot","lagcentral"),
  end
)

#' R implementation of `rxode2parse`/`rxode2` `linCmt()` solutions
#'
#' This interfaces the C/C++ linCmt() interface for use in R without
#' having to compile any `rxode2` models
#'
#' @param data input dataset, may have some parameters in the dataset
#'   itself; It also uses the standard `rxode2` dataset format which
#'   is translated as usual with `etTrans()`.  This procedure will
#'   produce the underlying `rxode2` model code and variables to
#'   translate the compartmental model
#'
#' @param ... Parameters that can be passed/parsed and integrated into
#'   the `data` for a linear compartment solution.  You can either
#'   have the parameters in the dataset itself or in this parameter
#'   block.  The type of linear compartmental model is determined by
#'   the parameters you specify in either the `data` or this extra set
#'   of arguments.
#'
#' @param fdepot The bioavailibility of doses to the depot.  Also
#'   affects rate/dur. The default is `1` or there is no
#'   bioavailibilty effect.  In `rxode2` you can specify this in a
#'   `linCmt()` model by adding `f(depot)=`.
#'
#' @param fcentral The bioavailibility of doses to the central compartment.  Also
#'   affects rate/dur. The default is `1` or there is no
#'   bioavailibilty effect.  In `rxode2` you can specify this in a
#'   `linCmt()` model by adding `f(central)=`.
#'
#' @param ratedepot The rate of doses that are modeled with `rate=-1`
#'   when dosed to the `depot` compartment.  This can be added to a
#'   `rxode2` `linCmt()` model by adding `rate(depot)=` to the model.
#'
#' @param ratecentral The rate of doses that are modeled with
#'   `rate=-1` when dosed to the `central` compartment.  This can be
#'   added to a `rxode2` `linCmt()` model by adding `rate(central)=`.
#'
#' @param durdepot The duration of infusion in the depot compartment
#'   with doses that are modeled with `rate=-2`.  This can be added to
#'   a `rxode2` `linCmt()` model by adding `dur(depot)=` to the model.
#'
#' @param durcentral The duration of infusions in the central
#'   compartment with doses that are modeled with `rate=-2`. This can
#'   be added to a `rxode2` `linCmt()` model by adding `dur(central)=`
#'   to the model.
#'
#' @param lagdepot The lag time of the depot compartment.  Can be
#'   specified in an `rxode2` model by `alag(depot)=`
#'
#' @param lagcentral he lag time of the depot compartment.  Can be
#'   specified in an `rxode2` model by `alag(central)=`
#'
#' @param scale Scaling factor for the final `Cc` output
#'
#' @param gradient If `TRUE`, this will add the gradients for the
#'   linear compartment models to the output, otherwise it only shows
#'   the concentration in the central compartment
#'
#' @inheritParams etTransParse
#'
#' @return A dataframe containing the linear compartment solution of
#'   the model appended as a column to the input data as `Cc`
#'
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' d <- nlmixr2data::nmtest
#'
#' names(d) <- sub("lagt", "lagcentral",
#'             sub("bioav", "fdepot",
#'             sub("rat2", "ratecentral",
#'             sub("dur2", "durcentral", names(d)))))
#'
#' ret <- linCmt(d, cl=1.1, v=20, ka=1.5)
#'
linCmt <- function(data, ...,
                   fdepot=1, fcentral=1,
                   ratedepot=0, ratecentral=0,
                   durdepot=0, durcentral=0,
                   lagdepot=0, lagcentral=0,
                   addlKeepsCov=FALSE, addlDropSs=TRUE, ssAtDoseTime=TRUE, scale=1,
                   gradient=FALSE,
                   keep=NULL,
                   covsInterpolation = c("locf", "linear", "nocb", "midpoint")) {
  checkmate::assertNumeric(fdepot, len=1, lower=0, finite=TRUE, any.missing=FALSE)
  checkmate::assertNumeric(fcentral, len=1, lower=0, finite=TRUE, any.missing=FALSE)
  checkmate::assertNumeric(ratedepot, len=1, lower=0, finite=TRUE, any.missing=FALSE)
  checkmate::assertNumeric(ratecentral, len=1, lower=0, finite=TRUE, any.missing=FALSE)
  checkmate::assertNumeric(durdepot, len=1, lower=0, finite=TRUE, any.missing=FALSE)
  checkmate::assertNumeric(durcentral, len=1, lower=0, finite=TRUE, any.missing=FALSE)
  checkmate::assertNumeric(scale, len=1, lower=0, finite=TRUE, any.missing=FALSE)
  checkmate::assertNumeric(lagdepot, len=1, lower=0, finite=TRUE, any.missing=FALSE)
  checkmate::assertNumeric(lagcentral, len=1, lower=0, finite=TRUE, any.missing=FALSE)
  checkmate::assertLogical(addlKeepsCov, len=1, any.missing=FALSE)
  checkmate::assertLogical(addlDropSs, len=1, any.missing=FALSE)
  checkmate::assertLogical(ssAtDoseTime, len=1, any.missing=FALSE)
  checkmate::assertLogical(gradient, len=1, any.missing=FALSE)
  covsInterpolation <- match.arg(covsInterpolation)
  # First take care of the possible inputs to lag time
  .linNamesData <- names(data)
  .w <- grepl(.rxLinInfoKaReg, .linNamesData, ignore.case=TRUE)
  if (length(.w) > 0L) {
    names(data)[.w] <- tolower(names(data)[.w])
    .lagDataInfo <- names(data)[.w]
    .linNamesData <- names(data)
  } else {
    .lagDataInfo <- character(0)
  }
  .lagParInfo <- list(fdepot=fdepot, fcentral=fcentral,
                      ratedepot=ratedepot, ratecentral=ratecentral,
                      durdepot=durdepot, durcentral=durcentral,
                      lagdepot=lagdepot, lagcentral=lagcentral)
  if (length(.lagDataInfo) > 0) {
    .keepLagParInfo <- names(.lagParInfo)[!(names(.lagParInfo) %in% .lagDataInfo)]
    .lagParInfo <- lapply(.keepLagParInfo,
                          function(x) {
                            .lagParInfo[[x]]
                          })
  }
  if (length(.lagParInfo) > 0) {
    for (v in names(.lagParInfo)) {
      data[[v]] <- .lagParInfo[[v]]
    }
  }

  .w <- which(grepl(.rxDerivedReg, .linNamesData, ignore.case = TRUE))
  if (length(.w) > 0L) {
    .linNamesData <- .linNamesData[.w]
  } else {
    .linNamesData <- character(0)
  }
  .lst <- list(...)
  .linNamePar <- names(.lst)
  .w <- which(grepl(.rxDerivedReg, .linNamePar, ignore.case = TRUE))
  if (length(.w) > 0) {
    .linNamePar <- .linNamePar[.w]
    for (.i in .linNamePar) {
      data[[tolower(.i)]] <- .lst[[.i]]
    }
  } else {
    .linNamePar <- character(0)
  }
  .base <- paste(c(.linNamesData, .linNamePar), collapse=", ")
  .trans <- eval(str2lang(paste0(".rxTransInfo(", .base, ")")))
  if (is.na(.trans$str["ka"])) {
    .trans$str <- c(.trans$str[1:6],
                    c("lagcentral"="lagcentral","fcentral"="fcentral","ratecentral"="ratecentral",
                      "durcentral"="durcentral"),
                    .trans$str[7],
                    c("lagdepot"=NA_character_,"fdepot"=NA_character_,
                      "ratedepot"=NA_character_,"durdepot"=NA_character_))
    .w <- which(grepl(.rxLinInfoKaReg, .linNamesData, ignore.case=TRUE), .linNamePar)
    if (length(.w) > 0) {
      .linNameExtra <- .linNamesData[.w]
      for (.n in .linNameExtra) {
        data[[tolower(.n)]] <- .lst[[.n]]
      }
    }
  } else {
    .trans$str <- c(.trans$str[1:6],
                    c("lagdepot"="lagdepot","fdepot"="fdepot",
                      "ratedepot"="ratedepot","durdepot"="durdepot"),
                    .trans$str[7],
                    c("lagcentral"="lagcentral","fcentral"="fcentral",
                      "ratecentral"="ratecentral", "durcentral"="durcentral"))
  }

  if (is.na(.trans$str["ka"])) {
    .mv <- rxode2parse(paste0("Cc=linCmt(", .base, ")\n",
                              "f(central)=fcentral\n",
                              "rate(central)=ratecentral\n",
                              "dur(central)=durcentral\n",
                              "alag(central)=lagcentral\n"), linear=TRUE)
  } else {
    .mv <- rxode2parse(paste0("Cc=linCmt(", .base, ")\n",
                              "f(depot)=fdepot\n",
                              "f(central)=fcentral\n",
                              "rate(depot)=ratedepot\n",
                              "rate(central)=ratecentral\n",
                              "dur(depot)=durdepot\n",
                              "dur(central)=durcentral\n",
                              "alag(depot)=lagdepot\n",
                              "alag(central)=lagcentral\n"), linear=TRUE)
  }
  data$rxRowNum <- seq_along(data[,1])
  .data <- as.data.frame(etTransParse(data, .mv, dropUnits=TRUE, allTimeVar = TRUE,
                                      addlKeepsCov=addlKeepsCov, addlDropSs=addlDropSs,
                                      ssAtDoseTime=ssAtDoseTime,
                                      keep = c("rxRowNum", keep)))
  .ret <- lapply(unique(.data$ID), function(i) {
    .dati <- .data[.data$ID==i,]
    .l <- length(.dati$TIME)
    .extra <- as.data.frame(setNames(lapply(names(.trans$str), function(n) {
      .t <- .trans$str[n]
      if (!is.na(.t) && !any(names(.dati) == .t)) {
        .t <- NA_character_
      }
      if (is.na(.t)) {
        if (grepl("^f", n)) {
          return(rep(1.0, .l))
        } else {
          return(rep(0.0, .l))
        }
      }
      .i <- seq_along(.dati[,.t])
      .v <- .dati[,.t]
      .w <- which(!is.na(.v))
      .yleft <- .v[.w[1]]
      .yright <- .v[.w[length(.w)]]
      if (covsInterpolation == "locf") {
        fun <- approxfun(.i[.w], .v[.w], method="constant",
                         yleft = .yleft, yright = .yright)
      } else if (covsInterpolation == "nocb") {
        fun <- approxfun(.i[.w], .v[.w], method="constant",
                         yleft = .yleft, yright = .yright,  f = 1)
      } else if (covsInterpolation == "midpoint") {
        fun <- approxfun(.i[.w], .v[.w], method="constant",
                         yleft = .yleft, yright = .yright,  f = 0.5)
      } else if (covsInterpolation == "linear") {
        fun <- approxfun(.i[.w], .v[.w], yleft = .yleft, yright = .yright)
      } else {
        stop("unknown covsInterpolation", call.=FALSE)
      }
      fun(.i)
    }), names(.trans$str)))
    .ret <- .Call(`_rxode2parse_compC`,
                  list(.dati[, c("TIME", "EVID", "AMT", "II")], .extra, .trans$trans), .mv)
    data.frame(ID=i, .ret)
  })
  dplyr::as_tibble(do.call(`rbind`, .ret))

}
#' Calculate the lambdas and coefficients of the two compartment model
#'
#' @param k10 elimination rate
#' @param k12 rate from central to peripheral compartment
#' @param k21 rate from peripheral to central compartment
#' @return List with `L` vector and matrices `C1` and `C2`
#' @export
#' @keywords internal
#' @author Matthew L. Fidler based on `wnl` package/paper, implemented
#'   in C/C++
#' @examples
#' .solComp2(k10=0.1, k12=3, k21=1)
.solComp2 <- function(k10, k12, k21) {
  checkmate::assertNumeric(k10, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(k12, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(k21, lower=0, len=1, any.missing=FALSE)
  .ret <- .Call(`_rxode2parse_solComp2`, k10, k12, k21)
  if (is.null(.ret)) {
    stop("roots must be distinct real values", call.=FALSE)
  }
  .ret
}
#' Calculate the lambdas and coefficients of the three compartment model
#'
#' @inheritParams .solComp2
#' @param k13 rate from central to peripheral compartment #2
#' @param k31 rate from peripheral compartment #2 to central
#' @return List with `L` vector and matrices `C1`, `C2` and `C3`
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
#' @examples
#' .solComp3(k10=0.1, k12=3, k21=1, k13=2, k31=0.5)
.solComp3 <- function(k10, k12, k21, k13, k31) {
  checkmate::assertNumeric(k10, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(k12, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(k21, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(k13, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(k31, lower=0, len=1, any.missing=FALSE)
  .ret <- .Call(`_rxode2parse_solComp3`, k10, k12, k21, k13, k31)
  if (is.null(.ret)) {
    stop("roots must be distinct real values", call.=FALSE)
  }
  .ret
}
#' Solve 1 point from a compartmental model
#'
#' @param inp Input compartment amounts.  The first value represents
#'   the depot compartment (if ka > 0) followed by the compartmetns
#'   (1st, 2nd or 3rd) or the amount in the first compartment followed
#'   by the next compartments in order (2nd, 3rd).
#' @param t The time after this point to sove for
#' @inheritParams .solComp3
#' @param rate Rate that is active in the central compartment
#' @param ka The oral absorption rate
#' @return All compartment values solved to t after current value
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
.solve1pt <- function(inp, t, k10, k12=0.0, k21=0.0, k13=0.0, k31=0.0,
                      v=1.0, rate=0.0, ka=0.0) {
  checkmate::assertNumeric(k10, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(k12, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(k21, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(k13, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(k31, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(v, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(inp, any.missing=FALSE, lower=0, min.len=1, max.len=4)
  checkmate::assertNumeric(ka, any.missing=FALSE, lower=0, min.len=1)
  if (ka > 0) {
    checkmate::assertNumeric(rate, any.missing=FALSE, lower=0, len=2)
  } else {
    checkmate::assertNumeric(rate, any.missing=FALSE, lower=0, len=1)
  }
  .Call(`_rxode2parse_solve1ptLin`,
        inp, t, ka, k10, k12, k21, k13, k31, v, rate)
}
