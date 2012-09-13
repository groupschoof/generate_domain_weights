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

extractChunksWithCompleteProteinTags <- function(path.to.file,
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

extractSingleProteinTags <- function( txt,
  prot.start.tag.regex="<protein",
  prot.end.tag.regex="</protein" ) {
  # Constructs a character vector with all protein tag entries found in
  # argument 'txt'. Each protein tag will be a single entry. Preceeding
  # and trailing 'junk' will be discarded. 
  #
  # Args:
  #  txt : The chunk of currently processed text to be parsed for protein tags.
  #
  #  prot.start.tag.regex : The regular expression used to find protein
  #                         start tags.
  #
  #  prot.end.tag.regex : The regular expression used to find protein
  #                       end tags.
  #
  # Returns: A character vector of all found protein tags.

  prot.beg.ind <- regexpr(prot.start.tag.regex, txt)[[1]]
  prot.end.ind <- regexpr(prot.end.tag.regex, txt)[[1]]

  # Return
  if ( prot.beg.ind < 0 || prot.end.ind < 0 ) {
    # CASE: Some rest text, that does not hold a complete protein tag.
    NULL
  } else {
    # CASE: At least a single protein tag encoded in current chunk of
    # text.
    prot.tag.txt <- substr( txt, prot.beg.ind, prot.end.ind )

    # More to parse?
    if ( prot.end.ind < nchar(txt) ) {
      # CASE: Trailing rest after the protein end tag in current chunk of text?
      all.results <- c( prot.tag.txt,
        extractSingleProteinTags(
          substr( txt, (prot.end.ind + 1), nchar(txt) ),
          prot.start.tag.regex, prot.end.tag.regex ))
      # return only non NULL entries:
      all.results[ ! is.null(all.results) ]
    } else {
      # CASE: Just the single protein tag found.
      prot.tag.txt
    }
  }
}

xmlUniprotInterProMatchProteinNodes <- function( xml.lines, protein.xpath="//protein" ) {
  xml.inp <- xmlInternalTreeParse( do.call( 'paste',
      as.list( c( "<dummyEncloseContentInSingleTag>",
          xml.lines, "</dummyEncloseContentInSingleTag>" ))))
  getNodeSet( xml.inp, protein.xpath )
}

