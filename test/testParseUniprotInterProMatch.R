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

# Downloaded test document to perform following tests
exmpl.doc <- xmlInternalTreeParse(
  project.file.path( "test","downloaded_iprmatch_A0A000.xml" ))
exmpl.doc.2 <- xmlInternalTreeParse(
  project.file.path( "test","downloaded_iprmatch_Q6GZX4.xml" ))
exmpl.doc.3 <- xmlInternalTreeParse(
  project.file.path( "test","downloaded_iprmatch_J006900.xml" ))
exmpl.doc.4 <- xmlInternalTreeParse(
  project.file.path( "test","downloaded_iprmatch_J000000.xml" ))
exmpl.doc.5 <- xmlInternalTreeParse(
  project.file.path( "test","downloaded_iprmatch_J000001.xml" ))

inter.pro.accessions <- c('IPR015421', 'IPR015422',
  'IPR004839', 'IPR015424', 'IPR010961')

# Test getProtein
print("Testing getProtein(...)")
prt.rslt <- getProtein(exmpl.doc)
checkTrue(!is.null( prt.rslt ))
checkEquals(xmlGetAttr(prt.rslt, 'id'), 'A0A000')
# NULL should be returned if document does not contain protein-tag:
fake <- getProtein(
  xmlInternalTreeParse("<tag id='my.tag'>Content</tag>"))
checkTrue(is.null(fake))
checkEquals( redisSPop("errors"),
  "Error: Document '<?xml version=\"1.0\"?>\n<tag id=\"my.tag\">Content</tag>\n ' did not contain a protein tag.")

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
ipr.match.parsed <- try(iprAnnotationPositionsMatrix(exmpl.doc), silent=T)
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
# IPR015424 contains IPR010961 
ngbs.rslt <- neighbors(ipr.match.parsed, 'IPR010961')
checkTrue( ! is.null(ngbs.rslt))
checkEquals( class(ngbs.rslt), 'list')
checkEquals( names(ngbs.rslt), c( 'start', 'end'))
checkEquals(ngbs.rslt$start, 'IPR015424' )
checkEquals(ngbs.rslt$end, 'IPR015424' )
# Test just two domains in protein
ngbs.rslt.2 <- neighbors( iprAnnotationPositionsMatrix(exmpl.doc.2),
  'IPR000666' )
checkTrue( ! is.null(ngbs.rslt.2) )
checkEquals( class(ngbs.rslt.2), 'list' )
checkEquals( names(ngbs.rslt.2), c('start', 'end') )
checkEquals( ngbs.rslt.2$start, NA )
checkEquals( ngbs.rslt.2$end, 'IPR015424' )
# Test two domains with identical positions should result in multiple 'start'
# neighbors
ngbs.rslt.3 <- neighbors( iprAnnotationPositionsMatrix(exmpl.doc.3),
  'IPR015424' )
checkTrue( ! is.null(ngbs.rslt.3) )
checkEquals( class(ngbs.rslt.3), 'list' )
checkEquals( names(ngbs.rslt.3), c('start', 'end') )
checkEquals( ngbs.rslt.3$start, c('IPR000666','IPR000999') )
checkEquals( ngbs.rslt.3$end, NA )

# Test parseUniprotIprMatchDocument
print( "Testing parseUniprotIprMatchDocument(...)" )
parse.rslt <- parseUniprotIprMatchDocument(exmpl.doc)
checkEquals( class(parse.rslt), 'matrix')
checkEquals( nrow(parse.rslt), 5)
checkEquals( ncol(parse.rslt), 2)
checkEquals( colnames(parse.rslt), c('start', 'end'))
checkEquals( rownames(parse.rslt), inter.pro.accessions )
# Assure NULL is returned and a try-error is logged, when a malformed
# XML document is attempted to be parsed:
parse.error <- parseUniprotIprMatchDocument( try(xmlInternalTreeParse("Foo Bar Baz"), silent=T) )
checkTrue( grepl('^Error', redisSPop("errors"), perl=T))
checkTrue( is.null(parse.error) )

# Test computeInterProDomainWeights
print("Testing computeInterProDomainWeights(...)")
computeInterProDomainWeights(parse.rslt)
checkEquals( redisSCard('interpro_domain_ids'), length(inter.pro.accessions) )
checkEquals( redisGet('IPR015422_cnt'), '1' )
checkEquals( redisSCard('IPR015422_nghbrs'), 2 )
computeInterProDomainWeights( parseUniprotIprMatchDocument(exmpl.doc.2) )
checkEquals( redisSCard('interpro_domain_ids'), length(inter.pro.accessions) + 1 )
checkEquals( redisGet('IPR015424_cnt'), '2' )
# And after computing another document, in which 'IPR015424' has TWO 'start'
# neighbors and no 'end' neighbors:
redisFlushAll()
computeInterProDomainWeights( parseUniprotIprMatchDocument(exmpl.doc.3) )
checkEquals( redisSCard('IPR015424_nghbrs'), 2 )
rks <- length( redisKeys() )
# Document with no valid InterPro Domains:
computeInterProDomainWeights( parseUniprotIprMatchDocument(exmpl.doc.4) )
checkEquals( length( redisKeys() ), rks )
# Document with only a single annotated InterPro Domain
computeInterProDomainWeights( parseUniprotIprMatchDocument(exmpl.doc.5) )
checkEquals( redisSCard('interpro_domain_ids'), 4 )
checkEquals( redisSCard('IPR000001_nghbrs'), 0 )

# Test logError
print("Testing logError(...)")
logError("Foo Bar Baz")
checkEquals(redisSPop("errors"), 'Foo Bar Baz')

# Clean up, girl:
redisFlushAll()
