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

parseEntry <- function(ipr.match.url) {

}


