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
src.project.file('src','uniprotInterProMatchDataIO.R')

# Print out Usage:
print( paste("Usage: Rscript",
    "--max-ppsize=500000",
    "computeInterProDomainWeights.R", 
    "path/to/split_of_match_complete.xml",
    "path/to/ordered_list_of_split_files.txt",
    "[http://redis.server.url ( default 'localhost' )",
    "redis_port ( default '6379' )]") )

# Read input
trailing.args <- commandArgs(trailingOnly=T)
if ( length(trailing.args) < 3 )
  stop("Missing arguments! See Usage for details!")

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

# Connect to redis
try.res <- try( redisConnect( host=redis.host, port=redis.port ) )
if ( identical(class(try.res), 'try-error') )
  stop("Could not connect to redis server")

# Function to be invoked with every parsed protein tag:              
uniprot.xml.docs <- c()
epf <- function(d) {
  uniprot.xml.docs <<- append( uniprot.xml.docs, d )
}

# Read data.
extractSingleProteinTags( 
  do.call( 'paste',
    as.list(
      extractChunksWithCompleteProteinTags(
        trailing.args[[1]], trailing.args[[2]], 
        readChunkFunc=readSplit,
        appendChunkFunc=appendSplitTillProteinEndTag(
          trailing.args[[1]], trailing.args[[2]]
        )
      )
    )
  ),
  each.prot.function=epf
)

# Log Message
print( paste("Parsing", length(uniprot.xml.docs), "Protein-Tags") )

# Start computation:
rslt <- unlist(
  lapply( uniprot.xml.docs,
    function(d) {
      # Parse and compute domain weights
      computeInterProDomainWeights( parseUniprotIprMatchDocument(d) )
    }
  ),
  use.names=F
)

# Debug:
# print(rslt)

# Feed result data into redis in a single transaction:
redisMulti()
not.printed.result <- lapply( rslt, eval )
redisExec()

# Clean up, boy:
redisClose()

# DONE
print("DONE")
