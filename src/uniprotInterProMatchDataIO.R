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
  # Returns: A character vector of all found protein tags or NULL.

  if ( is.null(txt) || is.na(txt) || length(txt) == 0 )
    return(NULL)

  prot.beg.ind <- regexpr(prot.start.tag.regex, txt)[[1]]
  print( prot.beg.ind )

  # RETURN:
  if ( prot.beg.ind < 0 ) {
    # CASE: No complete Protein Tag within this text chunk:
    NULL
  } else {
    # CASE: Found an opening Protein Tag, so cut off preceeding 'junk' until
    # first protein tag
    rest.txt <- substr( txt, prot.beg.ind, nchar(txt) )
    # Find index of first closing protein tag. Specifically its last char's
    # index.
    prot.end.match <- regexpr(prot.end.tag.regex, rest.txt)
    # No closing protein tag is an error, because function
    # 'extractChunksWithCompleteProteinTags' should not return a text chunk with
    # an opening protein tag but without a closing one.
    if (prot.end.match[[1]] < 0) {
      stop("Had a text chunk with an opening but no closing protein tag.")
    }
    prot.end.ind <- prot.end.match[[1]] + attr(prot.end.match, 'match.length')

    prot.xml.txt <- substr( rest.txt, 1, prot.end.ind )
    if ( prot.end.ind < nchar(rest.txt) ) {
      # CASE: More text to parse
      rslt <- c( prot.xml.txt, 
        extractSingleProteinTags(substr( rest.txt, (prot.end.ind + 1), nchar(rest.txt) )))
      # But only return non null entries:
      rslt[ ! is.null(rslt[]) ]
    } else {
      # CASE: Just the protein tag, already found
      prot.xml.txt
    }
  }
}

xmlUniprotInterProMatchProteinNodes <- function( xml.lines, protein.xpath="//protein" ) {
  xml.inp <- xmlInternalTreeParse( do.call( 'paste',
      as.list( c( "<dummyEncloseContentInSingleTag>",
          xml.lines, "</dummyEncloseContentInSingleTag>" ))))
  getNodeSet( xml.inp, protein.xpath )
}

