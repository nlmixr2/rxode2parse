if (requireNamespace("devtools", quietly = TRUE)) {
  .var <- deparse(rxode2::rxSupportedFuns())
  .var[1] <- paste0(".parseEnv$.parseFuns <- ", .var[1])
  .pf <- devtools::package_file("R/parseFuns.R")
  unlink(.pf)
  parseFuns.R <- file(.pf, "wb")
  writeLines(.var, parseFuns.R)
  close(parseFuns.R)

  message("rebuild rxode2parse_control.h")
  .l <- readLines(file.path(system.file("include", package="rxode2"), "rxode2_control.h"))
  .l <- gsub("rxode2_control", "rxode2parse_control", .l)
  .pf <- devtools::package_file("inst/include/rxode2parse_control.h")
  unlink(.pf)
  rxode2parse_control.h <- file(.pf, "wb")
  writeLines(.l, parseFuns.R)
  close(rxode2parse_control.h)
  message("done")
}

