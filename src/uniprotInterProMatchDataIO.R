library(XML) 

readSegment <- function(path.to.file,
  start.line.no, read.lines=10000) {
  scan(path.to.file, skip=( start.line.no - 1),
    nlines=read.lines, what=character(),
    sep="\n")
}

appendSegmentTillProteinEndTag <- function(path.to.file,
  start.line.no, read.lines, prot.end.tag.regex="</protein") {

  lines.chunk <- readSegment( path.to.file, start.line.no, read.lines )
  end.ind <- which(grepl(prot.end.tag.regex, lines.chunk, fixed=T))

  if ( length(end.ind) > 0 ) {
    lines.chunk[ 1 : min(end.ind) ]
  } else {
    c( lines.chunk,
      appendSegmentTillProteinEndTag( path.to.file,
        ( start.line.no + read.lines + 1 ),
        read.lines, prot.end.tag.regex ) )
  }
}

extractCompleteProteinTags <- function(path.to.file,
  start.line.no, read.lines=10000,
  prot.start.tag.regex="<protein",
  prot.end.tag.regex="</protein") {

  lines.chunk <- readSegment(path.to.file, start.line.no, read.lines)
  beg.ind <- which(grepl(prot.start.tag.regex, lines.chunk, fixed=T))

  if( length(beg.ind) == 0) {
    # CASE: No protein start tag found, do nothing.
    NULL
  } else {
    if ( grepl(prot.end.tag.regex, lines.chunk[ length(lines.chunk) ], fixed=T) ) {
      # CASE: Protein end tag is the very last line of current chunk
      lines.chunk[ min(beg.ind) : length(lines.chunk) ]
    } else {
      # Read into next chunk up to the first protein end tag
      c( lines.chunk[ min(beg.ind) : length(lines.chunk) ], 
        appendSegmentTillProteinEndTag( path.to.file, ( start.line.no + read.lines ),
          read.lines, prot.end.tag.regex ) )
    }
  }
}

xmlUniprotInterProMatchProteinNodes <- function( xml.lines, protein.xpath="//protein" ) {
  xml.inp <- xmlInternalTreeParse( do.call( 'paste',
      as.list( c( "<dummyEncloseContentInSingleTag>",
          xml.lines, "</dummyEncloseContentInSingleTag>" ))))
  getNodeSet( xml.inp, protein.xpath )
}

