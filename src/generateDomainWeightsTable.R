library(tools)
library(rredis)

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

# Print out Usage:
print( paste("Usage: Rscript",
    "generateDomainWeightsTable.R", 
    "number_proteins_computed",
    "path/to/domain_weights_table.tbl",
    "[http://redis.server.url ( default 'localhost' )",
    "redis_port ( default '6379' )]") )

# Read input
trailing.args <- commandArgs(trailingOnly=T)
if ( length(trailing.args) < 2 )
  stop("Missing arguments! See Usage for details!")

# Cardinality of protein set:
cardinality.protein.set <- as.integer( trailing.args[[1]] )

redis.host <- if ( length(trailing.args) > 2 ) {
                trailing.args[[3]] 
              } else {
                'localhost'
              }

redis.port <- if ( length(trailing.args) == 4 ) {
                as.integer(trailing.args[[4]])
              } else {
                6379
              }

# We need to enable large recursions:
options( expressions=500000 )

# Function to compute the Domain Weight from the number of annotations and
# number of distinct neighbor domains in all annotated proteins:
dom.wght <- function( no.annotations, no.distinct.neighbors,
  no.all.proteins=cardinality.protein.set ) {
  #
  # Computes the product of inverse abundance frequency (IAF) and inverse
  # versatility (IV), where for a given domain d:
  #
  # IAF(d) = |set of proteins| / |proteins with domain d|
  # IV(d)  = 1 / |distinct neighbor domains of d|
  # 
  # |…| is the set's cardinality. Note, that '|distinct neighbor domains of d|'
  # will be interpreted as 1 to avoid division by zero.
  #
  # Args:
  #   no.annotations        : The number of proteins the domain d has been
  #                           annotated in.
  #   no.distinct.neighbors : The number of distinct neighbor domains found in
  #                           all examined proteins.
  #   no.all.proteins       : Number of all examined proteins.
  #
  # Returns: A list with the following key value pairs: 'IAF', 'IV' and the
  # domain weight 'DW'.
  #   
  
  iaf <- log2( no.all.proteins / no.annotations )
  # No division by zero:
  ns <- if( no.distinct.neighbors == 0 ) 1 else no.distinct.neighbors
  iv <- ( 1 / ns )
  list('IAF'=iaf, 'IV'=iv, 'DW'=(iaf * iv))
}

# Connect to redis
try.res <- try( redisConnect( host=redis.host, port=redis.port ) )
if ( identical(class(try.res), 'try-error') )
  stop("Could not connect to redis server")

# Table of domain weights:
dw <- matrix(numeric(), ncol=3, dimnames=list(c(),c('IAF','IV','DW')))
# Read out processed InterPro domains and compute their respective domain
# weights:
while ( ! is.null( ipr.id <- redisSPop("interpro_domain_ids") ) ) {
  no.annos <- as.integer( redisGet( paste(ipr.id, "_cnt", sep="")) )
  no.nghbrs <- as.integer( redisSCard( paste(ipr.id, "_nghbrs", sep="")) )
  wghts <- dom.wght(no.annos, no.nghbrs)
  dw <<- rbind( dw,
      matrix(wghts, nrow=1, 
        dimnames=list(ipr.id, names(wghts)))
    )
}

# Save output:
write.table( dw, file=trailing.args[[2]] )

print("DONE")
