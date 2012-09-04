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

# Test downloadXmlDoc
print("Testing downloadXmlDoc(...)")
rslt <- downloadXmlDoc(ipr.match.url)
checkTrue(identical(class(rslt),
    c("XMLInternalDocument", "XMLAbstractDocument", "oldClass" )))
rslt <- downloadXmlDoc("http://non.existing.url")
checkTrue(identical(class(rslt), 'try-error'))

# Downloaded test documented to perform following tests
exmpl.doc <- xmlInternalTreeParse("test/downloaded_iprmatch_A0A000.xml")

# Test getProtein
print("Testing getProtein(...)")
prt.rslt <- getProtein(exmpl.doc)
checkTrue(!is.null( prt.rslt ))
checkEquals(xmlGetAttr(prt.rslt, 'id'), 'A0A000')

# Test getIprScnMatches
print("Testing getIprScnMatches(...)")
rslt.ipr.matches <- getIprScnMatches(prt.rslt)
checkTrue(length(rslt.ipr.matches) == 5)
checkEquals(xmlGetAttr(rslt.ipr.matches[[1]], 'id'),
  'G3DSA:3.40.640.10')

# Test interProAnnotation
print("Testing interProAnnotation(...)")
checkEquals(interProAnnotation(rslt.ipr.matches[[1]]),
  list(id='IPR015421', start='55', end='271'))

# Test parseEntry
print("Testing parseEntry(...)")
ipr.match.parsed <- try(parseEntry(ipr.match.url), silent=F)
checkTrue(!identical(class(ipr.match.parsed), 'try-error'))
