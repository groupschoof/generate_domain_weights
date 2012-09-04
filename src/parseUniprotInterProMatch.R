library(XML)
library(parallel)
library(rredis)

uniprotInterProMatchUrl <- function(accession) {
  paste(
    'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=iprmc;id=',
    accession,
    ';format=iprmcxml',
    sep=''
    )
}

downloadXmlDoc <- function(uniprot.url, noverbose=T) {
  try(xmlInternalTreeParse(uniprot.url), silent=noverbose)
}

iprDomCount <- function(ipr.id) {
  redis.id <- paste(ipr.id, '_cnt', sep='')
  if (is.null(redisGet(redis.id)))
    redisSet(redis.id, charToRaw('1'))
  else
    redisIncr(redis.id)
}

iprDomVersatility <- function(ipr.id, neighbor.domain.id) {
  redis.id <- paste(ipr.id, '_nghbrs', sep='')
  redisSAdd(redis.id, neighbor.domain.id)
}

getProtein <- function(xml.doc) {
  getNodeSet(xml.doc, "//protein")[[1]]
}

getIprScnMatches <- function(prot.node) {
  sapply(getNodeSet(prot.node, "//ipr[@type='Domain']"),
    xmlParent)
}

interProAnnotation <- function(match.node) {
  ipa <- list(id=xmlGetAttr(
      getNodeSet(match.node, "ipr")[[1]],
      "id"))
  lcn <- getNodeSet(match.node, "lcn")[[1]]
  ipa['start'] <- xmlGetAttr(lcn, 'start')
  ipa['end'] <- xmlGetAttr(lcn, 'end')
  # return
  ipa
}

parseEntry <- function(xml.doc) {
  sapply(getIprScnMatches(getProtein(xml.doc)),
    interProAnnotation)
}


