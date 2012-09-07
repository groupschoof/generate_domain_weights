library(XML)
library(rredis)
library(RCurl)

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

downloadUniprotDocsAndParse <- function(uniprot.accessions, noverbose=T) {
  uni.uris <- lapply( as.character( uniprot.accessions ),
    uniprotInterProMatchUrl )
  uni.docs <- getURL(uni.uris)
  uni.xmls <- lapply(uni.docs, function(d) {
      if( ! is.null(d) && ! is.na(d) &&
        ! identical(d,'') && ! grepl('^ERROR',d) ) 
      try( xmlInternalTreeParse(d), silent=noverbose )
    })
  # return only non null uni.xmls which includes ERRORs:
  ret.xmls <- uni.xmls[ ! as.logical( lapply(uni.xmls, is.null) ) ]
}

getProtein <- function(xml.doc) {
  p <- getNodeSet(xml.doc, "//protein")
  # return
  if ( length(p) == 0 ) {
    logError(paste("Error: Document '",
        toString.XMLNode( xml.doc ),
        "' did not contain a protein tag.",
        sep=''))
    NULL
  } else {
    p[[1]]
  }
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
  p <- getProtein(xml.doc)
  if ( ! is.null(p) )
    do.call('rbind', 
  lapply(getIprScnMatches( p ),
  interProAnnotation))
}

neighbors <- function(position.matrix, row.name) {
  # Exclude the row named <row.name> from position.matrix
  m <- position.matrix[rownames(position.matrix) != row.name, , drop=F ]
  # Single rowed matrix of all 'start' and 'end' positions where colnames are
  # the domain ids.
  pos <- matrix( c( m[,'start'], m[,'end'] ), nrow=1, 
    dimnames=list(c('pos'), c(rownames(m), rownames(m))) )
  # Result is a list with start and end neighbor for domain <row.name>
  setNames(
    lapply(locations, function(location) {
      # Current position to measure distances to:
      cp <- position.matrix[[ row.name, location ]]
      # Measure distance for each postion from left to start or from end to
      # right? The latter results in inversion of sign.
      dir.fct <- if (identical(location, 'start')) 1 else (-1)
      # Calculate distances:
      diff.pos <- dir.fct * (cp - pos)
      # Use only positive distances - look in the right direction, that is:
      dists <- diff.pos[ , diff.pos[,] >= 0, drop=F ]
      # Select those closest in the currently evaluated direction ( depends on
      # <location> ). No entries should be returned as 'NA', the matching
      # domain ids otherwise.
      if ( length(dists) > 0 ) colnames( dists[ , dists == min(dists), drop=F ] ) else NA
    }),
    locations)
}

parseUniprotIprMatchDocument <- function(xml.doc) {
  if ( ! is.null(xml.doc) && ! identical(class(xml.doc), 'try-error') ) {
    iapm <- iprAnnotationPositionsMatrix(xml.doc)
    if ( ! is.null(iapm) ) {
      if ( nrow(iapm) == 1 ) {
        matrix( c(NA, NA), nrow=1, dimnames=list(rownames(iapm), locations) )
      } else {
        do.call('rbind', 
          lapply( rownames(iapm), function(row.name) {
              ngbs <- neighbors(iapm, row.name)
              matrix( ngbs, nrow=1, 
                dimnames=list(c(row.name), names(ngbs))
                )
            })
          )
      }
    }
  } else {
    logError( xml.doc )
    # return
    NULL
  }
}

iprDomAccessions <- function(ipr.id,
  accessions.set.key='interpro_domain_ids') {
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

logError <- function(err.message) {
  redisSAdd("errors", err.message)
}

computeInterProDomainWeights <- function(ipr.domain.matrix) {
  if ( ! is.null( ipr.domain.matrix ) ) {
  # Compute domain weights for all domains annotated in the currently processed
  # document:
  lapply( rownames(ipr.domain.matrix), function(ipr.id) {
      # Remember all distinct domain accessions for convenience
      iprDomAccessions(ipr.id)
      # Count number of occurrences of current domain:
      iprDomCount(ipr.id)
      # Count neighbors at both positions 'start' and 'end':
      lapply(locations, function(loc) {
          # Count each distinct neighbor domain found (there can be more than
          # one at both the 'start' and 'end' positions):
          lapply( ipr.domain.matrix[[ipr.id, loc]], function(ngb.id) {
              if ( ! is.na(ngb.id) )
                iprDomVersatility(ipr.id, ngb.id)
              })
          })
      })
  }
  # No error, so return:
  TRUE
}

wasBusy <- function(d) {
  ns <- getNodeSet(d, "//h2")
  (length(ns) == 1 && identical(xmlValue(ns[[1]]), "Server Too Busy"))
}

findServerBusyResults <- function(xml.docs) {
  xml.docs[ as.logical( lapply(xml.docs, wasBusy) ) ]
}
