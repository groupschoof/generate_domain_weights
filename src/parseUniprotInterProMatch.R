library(XML)
library(rredis)
library(RCurl)

# CONSTANTS:
locations <- c('start','end')

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
  # Parses a Uniprot InterPro-Match XMl Document containing a single Protein
  # tag. Parsing extracts information of annotated InterPro domains and their
  # neighbors. 
  #
  # Arguments:
  #    xml.doc : The single protein tag as parsed by xmlInternalTreeParse (library 'XML')
  #
  # Returns: Matrix in which each row is named with an InterPro domain ID and
  # each column holds the 'start' and 'end' neighbors for the current row's
  # InterPro domain. This is of course in the context of the argument protein.
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
  # R-Code to add an InterPro-ID to the set of InterPro-IDs processed. NOTE:
  # Function returns quoted method call. Do not forget to 'eval(...)' it in a
  # redis transaction.
  #
  # Arguments:
  #    ipr.id : The InterPro-ID to store
  #    accessions.set.key : The key for the redis set to store all InterPro-IDs in.
  #
  # Returns: The quoted redis method call, so it can be executed later in a transaction.
  #
  redisSAdd(accessions.set.key, charToRaw(ipr.id))
}

iprDomCount <- function(ipr.id) {
  # R-Code to count number of Protein Annotations with the argument
  # InterPro-ID. NOTE: Function returns quoted method call. Do not forget to
  # 'eval(...)' it in a redis transaction.
  #
  # Arguments:
  #    ipr.id : The InterPro-ID to store
  #
  # Returns: The quoted redis method call, so it can be executed later in a transaction.
  #
  redis.id <- paste(ipr.id, '_cnt', sep='')
  if (is.null(redisGet(redis.id)))
    redisSet(redis.id, charToRaw('1'))
  else
    redisIncr(redis.id)
}

iprDomVersatility <- function(ipr.id, neighbor.domain.id) {
  # R-Code to count number of distinct direct neighbor domains in all Protein
  # Annotations with the argument InterPro-ID. NOTE: Function returns quoted
  # method call. Do not forget to 'eval(...)' it in a redis transaction.
  #
  # Arguments:
  #    ipr.id : The InterPro-ID to store found neighbor for.
  #    neighbor.domain.id : The accounted neighbor domain to store in the
  #                         corresponding reids set of ipr.id's neighbors.
  #
  # Returns: The quoted redis method call, so it can be executed later in a transaction.
  #
  redis.id <- paste(ipr.id, '_nghbrs', sep='')
  redisSAdd(redis.id, charToRaw(neighbor.domain.id))
}

logError <- function(err.message) {
  # Writes error message to redis error list.
  #
  # Arguments:
  #    err.message : The error message to store.
  #
  # Returns: Nothing
  redisLPush( "errors", err.message )
}

computeInterProDomainWeights <- function(ipr.domain.matrix) {
  # Generate R-code to be executed in a redis transaction. This R-code will
  # keep track of all encountered InterPro-IDs as well as count their number of
  # annotations and the number of distinct neighbor domains for each found
  # domain. NOTE: Remember to 'eval(...)' this function's returned list of
  # quoted R expressions.
  #
  # Arguments:
  #    ipr.domain.matrix : The matrix of all found InterPro domains, as
  #                        returned by 'parseUniprotIprMatchDocument'.
  #
  # Returns: Vector of quoted R-code to be executed in a single redis transaction.
  # 
  if ( ! is.null( ipr.domain.matrix ) ) {
    # Return:
    unlist(
      lapply( rownames(ipr.domain.matrix), function(ipr.id) {
          list(
            # Remember all distinct domain accessions for convenience
            bquote( iprDomAccessions( .(ipr.id) ) ),
            # Count number of occurrences of current domain:
            bquote( iprDomCount( .(ipr.id) ) ),
            # Count neighbors at both positions 'start' and 'end':
            lapply(locations, function(loc) {
                # Count each distinct neighbor domain found (there can be more than
                # one at both the 'start' and 'end' positions):
                lapply( ipr.domain.matrix[[ipr.id, loc]], function(ngb.id) {
                    if ( ! is.na(ngb.id) )
                      bquote( iprDomVersatility( .(ipr.id), .(ngb.id) ) )
                })
            })
          )
      })
    )
  }
}
