library(tools)
library(parallel)
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

# Print out Usage:
print( paste("Usage: Rscript",
    "computeInterProDomainWeightsForAccessions.R", 
    "one_uniprot_accession_per_line.txt",
    "[http://redis.server.url default 'localhost']",
    "[redis_port] default '6379'") )

# Read input
trailing.args <- commandArgs(trailingOnly=T)
if ( length(trailing.args) < 1 )
  stop("Missing arguments! See Usage for details!")

accessions <- as.character( read.table(trailing.args[[1]])$V1 )

redis.host <- if ( length(trailing.args) == 3 ) {
                trailing.args[[2]] 
              } else {
                'localhost'
              }

redis.port <- if ( length(trailing.args) == 3 ) {
                as.integer(trailing.args[[3]])
              } else {
                6379
              }

# Read data.
uniprot.xml.docs <- downloadUniprotDocsAndParse( accessions, noverbose=F )

# Start computation in parallel
rslt <- mclapply( uniprot.xml.docs, function(d) {
    # Connect to redis
    try.res <- try( redisConnect( host=redis.host, port=redis.port ) )
    if ( identical(class(try.res), 'try-error') )
      stop("Could not connect to redis server")
    
    # Parse and compute domain weights
    computeInterProDomainWeights( parseUniprotIprMatchDocument(d) )
    
    # Clean up, boy:
    redisClose()

    }, mc.preschedule=T, mc.cores=detectCores() )

# DONE
print("DONE")
