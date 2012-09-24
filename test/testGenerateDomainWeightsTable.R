library(RUnit)
library(tools)
# In R sourcing other files is not trivial, unfortunately.
# WARNING:
# This method ONLY works for project files in depth one sub dirs!
project.file.path <- function(...) {
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.dir <- dirname(file_path_as_absolute(script.name))
  project.dir <- sub(basename(script.dir),'',script.dir)
  normalizePath(file.path(project.dir,...))
}
src.project.file <- function(...) {
  source(project.file.path(...))
}
src.project.file('src','parseUniprotInterProMatch.R')

# Connect to redis
rc <- try(redisConnect(), silent=T)
if (identical(class(rc), 'try-error')) 
  stop("Could not connect to redis server. Did you start it?")
redisFlushAll()

# Generate test set:
iprDomAccessions("IPR000666")
iprDomAccessions("IPR000999")
iprDomAccessions("IPR010961")

iprDomCount("IPR000666")
iprDomCount("IPR000666")
iprDomCount("IPR000666")
iprDomCount("IPR000999")
iprDomCount("IPR000999")
iprDomCount("IPR010961")

iprDomVersatility("IPR000666", "IPR000999")
iprDomVersatility("IPR000666", "IPR010961")
iprDomVersatility("IPR000999", "IPR010961")

# Execute script:
scr.path <- project.file.path( "src", "generateDomainWeightsTable.R" )
system( paste("Rscript", scr.path, "50000000", "test_out.tbl") )

# Read result:
rslt <- read.table("test_out.tbl")
# print("RSLT")
# print( rslt[order(rownames(rslt)),])
expected <- as.data.frame(matrix(
  c(log2(50000000 / 3), 1/2, (0.5 * log2(50000000 / 3)),
    log2(50000000 / 2), 1/1, log2(50000000 / 2),
    log2(50000000 / 1), 1.0, log2(50000000 / 1)),
  nrow=3, byrow=T,
  dimnames=list(c("IPR000666", "IPR000999", "IPR010961"),
    c('IAF', 'IV', 'DW'))
  ))
# print( "EXPC" )
# print( expected[order(rownames(expected)),])
checkEquals(  rslt[order(rownames(rslt)),], expected[order(rownames(expected)),] )

# Clean up, boy:
redisFlushAll()
unlink("test_out.tbl")
