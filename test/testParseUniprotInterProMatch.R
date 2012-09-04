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

# Downloaded test document to perform following tests
exmpl.doc <- xmlInternalTreeParse("test/downloaded_iprmatch_A0A000.xml")
inter.pro.accessions <- c('IPR015421', 'IPR015422',
  'IPR004839', 'IPR015424', 'IPR010961')

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
  matrix(c(start=55, end=271), nrow=1,
    dimnames=list(c('IPR015421'), c('start', 'end'))))

# Test iprAnnotationPositionsMatrix
print("Testing iprAnnotationPositionsMatrix(...)")
ipr.match.parsed <- try(iprAnnotationPositionsMatrix(exmpl.doc), silent=F)
checkTrue(!identical(class(ipr.match.parsed), 'try-error'))
checkEquals(class(ipr.match.parsed), 'matrix')
checkTrue(nrow(ipr.match.parsed) == 5)
checkTrue(ncol(ipr.match.parsed) == 2)

# Test neighbors
print("Testing neighbors(...)")
ngbs.rslt <- neighbors(ipr.match.parsed, 'IPR015421')
checkTrue( ! is.null(ngbs.rslt))
checkEquals( class(ngbs.rslt), 'list')
checkEquals( names(ngbs.rslt), c( 'start', 'end'))
checkEquals(ngbs.rslt$start, 'IPR004839')
checkEquals(ngbs.rslt$end, 'IPR015422')
ngbs.rslt <- neighbors(ipr.match.parsed, 'IPR015424')
checkTrue( ! is.null(ngbs.rslt))
checkEquals( class(ngbs.rslt), 'list')
checkEquals( names(ngbs.rslt), c( 'start', 'end'))
checkEquals(ngbs.rslt$start, NA)
checkEquals(ngbs.rslt$end, NA)

# Test parseUniprotIprMatchDocument
print("Testing parseUniprotIprMatchDocument(...)")
parse.rslt <- parseUniprotIprMatchDocument(exmpl.doc)
checkEquals( class(parse.rslt), 'matrix')
checkEquals( nrow(parse.rslt), 5)
checkEquals( ncol(parse.rslt), 2)
checkEquals( colnames(parse.rslt), c('start', 'end'))
checkEquals( rownames(parse.rslt), inter.pro.accessions )

# Test computeInterProDomainWeights
print("Testing computeInterProDomainWeights(...)")
computeInterProDomainWeights(parse.rslt)
checkEquals( redisSCard('ipr_accessions'), length(inter.pro.accessions) )
checkEquals( redisGet('IPR015422_cnt'), '1' )
checkEquals( redisSCard('IPR015422_nghbrs'), 2 )
