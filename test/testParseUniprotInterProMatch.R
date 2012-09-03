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

# Test uniprotInterProMatchUrl
print("Testing uniprotInterProMatchUrl(...)")
ipr.match.url <- uniprotInterProMatchUrl('A0A000')
checkEquals(ipr.match.url,
  'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=iprmc;id=A0A000;format=iprmcxml')

# Test parseEntry
print("Testing parseEntry(...)")
ipr.match.parsed <- try(parseEntry(ipr.match.url), silent=F)
checkTrue(!identical(class(ipr.match.parsed), 'try-error')



