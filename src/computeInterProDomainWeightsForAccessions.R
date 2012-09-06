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

accessions <- read.table(trailing.args[[1]])$V1

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


# Start computation in parallel
mclapply( accessions, function(uniprot.acc) {
    # Connect to redis
    try.res <- try( redisConnect( host=redis.host, port=redis.port ) )
    if ( identical(class(try.res), 'try-error') )
    stop("Could not connect to redis server")

    # Read data
    u <- uniprotInterProMatchUrl( as.character(uniprot.acc) )
    d <- downloadXmlDoc( u )
    
    # Parse and compute domain weights
    if ( ! identical(class(d), 'try-error') )
      computeInterProDomainWeights( parseUniprotIprMatchDocument(d) )
    else
      print( paste("URL", u, "did return an error") )
    
    } )
