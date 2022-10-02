.makeGrammar <- function() {
  message("Update Parser h file");
  dparser::mkdparse(devtools::package_file("inst/tran.g"),
                    devtools::package_file("src/"),
                    grammar_ident="rxode2parse")
  file <- gsub("^([#]line [0-9]+ )\".*(src)/+(.*)\"","\\1\"\\2/\\3\"",
               readLines(devtools::package_file("src/tran.g.d_parser.c")))
  sink(devtools::package_file("src/tran.g.d_parser.h"))
  cat(paste(file,collapse="\n"))
  cat("\n")
  sink()
  unlink(devtools::package_file("src/tran.g.d_parser.c"))
}
