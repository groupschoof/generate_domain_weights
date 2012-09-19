library(XML) 

readSegment <- function(path.to.file,
  start.line.no, read.lines=10000) {
  scan(path.to.file, skip=( start.line.no - 1),
    nlines=read.lines, what=character(),
    sep="\n")
}

readSplit <- function(path.to.split.file, ...) {
  scan(path.to.split.file, sep="\n", what=character())
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
        ( start.line.no + read.lines ),
        read.lines, prot.end.tag.regex ) )
  }
}

nextSplitFile <- function(path.to.curr.split,
  path.to.sorted.split.list) {
  fl.name <- basename(path.to.curr.split)
  sl <- as.character(read.table(path.to.sorted.split.list)$V1)
  normalizePath(file.path(
      dirname(path.to.curr.split),
      sl[ which(sl[] == as.character(fl.name)) + 1 ]
  ))
}

appendSplitTillProteinEndTag <- function(path.to.curr.split,
 path.to.sorted.split.list,
 prot.end.tag.regex="</protein") {

  path.to.file <- nextSplitFile( path.to.curr.split,
    path.to.sorted.split.list )
  lines.chunk <- readSplit( path.to.file )
  end.ind <- which(grepl(prot.end.tag.regex, lines.chunk, fixed=T))

  if ( length(end.ind) > 0 ) {
    lines.chunk[ 1 : min(end.ind) ]
  } else {
    c( lines.chunk,
      appendSplitTillProteinEndTag( path.to.curr.split,
        path.to.sorted.split.list,
        prot.end.tag.regex ) )
  }
}

extractChunksWithCompleteProteinTags <- function(path.to.file,
  start.line.no, read.lines=10000,
  prot.start.tag.regex="<protein",
  prot.end.tag.regex="</protein",
  append.lines=400,
  readChunkFunc=readSegment,
  appendChunkFunc=appendSegmentTillProteinEndTag(
      path.to.file,
      ( start.line.no + read.lines ),
      append.lines,
      prot.end.tag.regex 
    )
  ) {
  # First line number is always 1, never less:
  if( start.line.no < 1 )
    start.line.no <- 1

  lines.chunk <- readChunkFunc(path.to.file, start.line.no, read.lines)
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
        eval( appendChunkFunc ) )
    }
  }
}

extractSingleProteinTags <- function( txt,
  prot.start.tag.regex="<protein",
  prot.end.tag.regex="</protein",
  noverbose=T,
  each.prot.function) {
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
  #  noverbose : Defines level of logging warnings and errors. If true all try
  #              statements are executed in silent mode.
  #
  #  each.prot.function : A function to be invoked with each generated protein
  #                       tag as argument. This argument is the result of xmlInternalTreeParse.
  #
  # Returns: A character vector of all found protein tags or NULL.
  #

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

  # If an error occurred log the txt chunk that caused it.
  if ( identical( class(parse.res), 'try-error' ) ) {
    parse.res[1] <- paste(parse.res[1], "TXT-Chunk that caused the error:", txt, sep="\n") 
  }

  # Invoke arbitrary function with parsed protein tag:
  each.prot.function( parse.res )

  # Is there more to parse and process?
  if ( prot.end.ind < nchar(rest.txt) ) {
      extractSingleProteinTags(
        substr( rest.txt, (prot.end.ind + 1), nchar(rest.txt)),
        prot.start.tag.regex,
        prot.end.tag.regex,
        noverbose,
        each.prot.function
        )
  }

  # Return
  NULL
}
