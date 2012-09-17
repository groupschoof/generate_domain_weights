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
  prot.end.tag.regex="</protein",
  append.lines=400) {

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
        appendSegmentTillProteinEndTag( path.to.file, ( start.line.no + read.lines + 1 ),
          append.lines, prot.end.tag.regex ) )
    }
  }
}

extractSingleProteinTags <- function( txt,
  prot.start.tag.regex="<protein",
  prot.end.tag.regex="</protein",
  noverbose=F) {
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
  #  prot.end.tag.regex : The regular expression used to find protein end tags.
  #                       (NOTE: Match will be one character longer than the
  #                       regular expression itself: 
  #                       So '</protein>' not '</protein' !)
  #
  # Returns: A character vector of all found protein tags or NULL.

  # Only process valid input:
  if ( is.null(txt) || is.na(txt) || length(txt) == 0 )
    return(NULL)

  # Look for a protein begin tag:
  prot.beg.ind <- regexpr(prot.start.tag.regex, txt)[[1]]

  # CASE: No complete Protein Tag within this text chunk:
  if ( prot.beg.ind < 0 ) {
    return(NULL)
  } 

  # CASE: Found an opening Protein Tag, so cut off preceeding 'junk' until
  # first protein tag
  rest.txt <- substr( txt, prot.beg.ind, nchar(txt) )
  # Find index of first closing protein tag. Specifically its last char's
  # index.
  prot.end.match <- regexpr(prot.end.tag.regex, rest.txt)
  # No closing protein tag means dealing with a chunk, where the protein end
  # tag is on the same line as the start tag. (The next chunk of text will be
  # dealing with this.)
  if (prot.end.match[[1]] < 0) {
    return(NULL)
  }
  prot.end.ind <- prot.end.match[[1]] + attr(prot.end.match, 'match.length')

  # Extract complete Protein Tag and parse this XML:
  parse.res <- try( xmlInternalTreeParse(
      substr( rest.txt, 1, prot.end.ind )),
    silent=noverbose )
  prot.xml.txt <- if ( identical( class(parse.res), 'try-error' ) ) {
    print( paste("XML parsing of the following TXT-Chunk caused an error:", txt, sep="\n") ) 
    NULL
  } else {
    parse.res
  }

  if ( prot.end.ind < nchar(rest.txt) ) {
    # CASE: More text to parse
    c( prot.xml.txt, 
      extractSingleProteinTags(substr( rest.txt, (prot.end.ind + 1), nchar(rest.txt) )))
  } else {
    # CASE: Just the protein tag, already found
    prot.xml.txt
  }
}
