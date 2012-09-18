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
src.project.file('src','uniprotInterProMatchDataIO.R')

# Print out Usage:
print( paste("Usage: Rscript",
    "computeInterProDomainWeights.R", 
    "path/to/match_complete.xml",
    "start.line.number",
    "no.of.lines.to.read",
    "[http://redis.server.url ( default 'localhost' )",
    "redis_port ( default '6379' )]") )

# Read input
trailing.args <- commandArgs(trailingOnly=T)
if ( length(trailing.args) < 3 )
  stop("Missing arguments! See Usage for details!")

start.line.number <- as.integer( trailing.args[[2]] )
no.of.lines.to.read <- as.integer( trailing.args[[3]] )

redis.host <- if ( length(trailing.args) == 5 ) {
                trailing.args[[4]] 
              } else {
                'localhost'
              }

redis.port <- if ( length(trailing.args) == 5 ) {
                as.integer(trailing.args[[5]])
              } else {
                6379
              }

# Function to be invoked with every parsed protein tag:              
uniprot.xml.docs <- c()
epf <- function(d) {
  uniprot.xml.docs <<- append( uniprot.xml.docs, d )
}

# Read data.
extractSingleProteinTags( 
  do.call( 'paste', as.list( extractChunksWithCompleteProteinTags(
    trailing.args[[1]], start.line.number, no.of.lines.to.read))),
  each.prot.function=epf
)

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
  },
  mc.preschedule=T, mc.cores=detectCores()
)

# DONE
print("DONE")
