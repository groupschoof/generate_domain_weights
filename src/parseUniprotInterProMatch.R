library(XML)
library(rredis)

# CONSTANTS:
locations <- c('start','end')

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

getProtein <- function(xml.doc) {
  getNodeSet(xml.doc, "//protein")[[1]]
}

getIprScnMatches <- function(prot.node) {
  sapply(getNodeSet(prot.node, "//ipr[@type='Domain']"),
    xmlParent)
}

interProAnnotation <- function(match.node) {
  ipr.id <- xmlGetAttr(
      getNodeSet(match.node, "ipr")[[1]],
      "id")
  lcn <- getNodeSet(match.node, "lcn")[[1]]
  st <- as.integer(xmlGetAttr(lcn, 'start'))
  en <- as.integer(xmlGetAttr(lcn, 'end'))
  matrix(c(st, en), nrow=1, dimnames=list(c(ipr.id),
      c('start', 'end')))
}

iprAnnotationPositionsMatrix <- function(xml.doc) {
  do.call('rbind', 
    lapply(getIprScnMatches(getProtein(xml.doc)),
      interProAnnotation))
}

neighbors <- function(position.matrix, row.name) {
  m <- position.matrix[rownames(position.matrix) != row.name, ]
  pos <- sort(c(m[,'start'], m[,'end']))
  setNames(
    lapply(locations, function(location) {
      diff.pos <- position.matrix[[ row.name, location ]] - pos
      if (identical(location, 'end'))
        diff.pos <- diff.pos * (-1)
      indices <- which(diff.pos >= 0)
      if ( length(indices) > 0 )
        names(which(diff.pos == min(diff.pos[ indices ])))
      else
        NA
    }),
    locations)
}

parseUniprotIprMatchDocument <- function(xml.doc) {
  iapm <- iprAnnotationPositionsMatrix(xml.doc)
  do.call('rbind', 
    lapply( rownames(iapm), function(row.name) {
        ngbs <- neighbors(iapm, row.name)
        matrix( ngbs, nrow=1, 
          dimnames=list(c(row.name), names(ngbs))
          )
      })
    )
}

iprDomAccessions <- function(ipr.id,
  accessions.set.key='ipr_accessions') {
  redisSAdd(accessions.set.key, charToRaw(ipr.id))
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
  redisSAdd(redis.id, charToRaw(neighbor.domain.id))
}

computeInterProDomainWeights <- function(ipr.domain.matrix) {
  lapply( rownames(ipr.domain.matrix), function(ipr.id) {
      iprDomAccessions(ipr.id)
      iprDomCount(ipr.id)
      lapply(locations, function(loc) {
          ngb.id <- ipr.domain.matrix[[ipr.id, loc]]
          if ( ! is.na(ngb.id) )
            iprDomVersatility(ipr.id, ngb.id)
          })
      })
  # return
  TRUE
}
