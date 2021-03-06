library(RUnit)
library(tools)
# In R sourcing other files is not trivial, unfortunately.
# WARNING:
# This method ONLY works for project files in depth one sub dirs!
project.file.path <- function(...) {
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.dir <- dirname(file_path_as_absolute(script.name))
  project.dir <- sub(basename(script.dir),'',script.dir)
  normalizePath(file.path(project.dir,...))
}
src.project.file <- function(...) {
  source(project.file.path(...))
}
src.project.file('src','uniprotInterProMatchDataIO.R')

# Test extractChunksWithCompleteProteinTags
print("Testing extractChunksWithCompleteProteinTags(...)")
uni.ipr.scn.path <- project.file.path("test",
    "testUniprotInterProCompleteMatchFile.xml")
extr.lines <- extractChunksWithCompleteProteinTags( uni.ipr.scn.path, 1, 25 )
checkEquals( extr.lines,
  scan(uni.ipr.scn.path, sep="\n", nlines=41, what=character()) )

# Test no start and no end tag
extr.lines <- extractChunksWithCompleteProteinTags( uni.ipr.scn.path, 5, 20 )
checkTrue( is.null(extr.lines) ) 

# Test no start but an end tag
extr.lines <- extractChunksWithCompleteProteinTags( uni.ipr.scn.path, 20, 21 )
checkTrue( is.null(extr.lines) ) 

# Test a start but no end tag
extr.lines <- extractChunksWithCompleteProteinTags( uni.ipr.scn.path, 42, 20 )
checkEquals( extr.lines,
  scan(uni.ipr.scn.path, sep="\n", skip=41, nlines=31, what=character()) )

# Test a complete protein tag
extr.lines <- extractChunksWithCompleteProteinTags( uni.ipr.scn.path, 42, 30 )
checkEquals( extr.lines,
  scan(uni.ipr.scn.path, sep="\n", skip=41, nlines=31, what=character()) )

# Test multiple protein tags, last without end tag
extr.lines <- extractChunksWithCompleteProteinTags( uni.ipr.scn.path, 1, 65 )
checkEquals( extr.lines,
  scan(uni.ipr.scn.path, sep="\n", nlines=72, what=character()) )

# Test extractSingleProteinTags
print("Testing extractSingleProteinTags(...)")
extr.xml.txt.2 <- c()
epf <- function(p) {
  if ( ! identical( class(p), 'try-error' ) )
    extr.xml.txt.2 <<- append(extr.xml.txt.2, p)
}
uni.ipr.scn.path.2 <- project.file.path( "test", "testUniprotInterProCompleteMatchFile_2.xml" )
extr.lines.2 <- do.call( 'paste', as.list( extractChunksWithCompleteProteinTags( uni.ipr.scn.path.2, 1, 20 ) ) )
extractSingleProteinTags( extr.lines.2, each.prot.function=epf )
checkTrue( ! is.null( extr.xml.txt.2 ) )
checkEquals( length( extr.xml.txt.2 ), 1 ) 

# End Tag on Same Line as Start Tag and just read until that very line: Also
# one mal formed match tag in first protein:
extr.lines.2 <- do.call( 'paste', as.list( extractChunksWithCompleteProteinTags( uni.ipr.scn.path.2, 1, 70 ) ) )
extr.xml.txt.2 <- c()
extractSingleProteinTags( extr.lines.2, each.prot.function=epf )
checkTrue( ! is.null( extr.xml.txt.2 ) )
checkEquals( length( extr.xml.txt.2 ), 1 ) 
# # Second entry must be a try-error:
# checkEquals( class( extr.xml.txt.2[[2]] ), 'try-error' )

# Document excluding last entry has two valid proteins
extr.lines.2 <- do.call( 'paste', as.list( extractChunksWithCompleteProteinTags( uni.ipr.scn.path.2, 1, 71 ) ) )
extr.xml.txt.2 <- c()
extractSingleProteinTags( extr.lines.2, each.prot.function=epf )
checkTrue( ! is.null( extr.xml.txt.2 ) )
checkEquals( length( extr.xml.txt.2 ), 2 ) 

# Test using splits:
split.A.path <- project.file.path( "test", "splitA" )
split.list.path <- project.file.path( "test", "splitList" )
extr.lines.2 <- do.call( 'paste',
  as.list( extractChunksWithCompleteProteinTags(
    split.A.path, split.list.path,
    readChunkFunc=readSplit, 
    appendChunkFunc=appendSplitTillProteinEndTag(
      split.A.path, split.list.path)
    )
  )
)
extr.xml.txt.2 <- c()
extractSingleProteinTags( extr.lines.2,
  each.prot.function=epf )
checkTrue( ! is.null( extr.xml.txt.2 ) )
# Only a single valid protein in this chunk:
checkEquals( length( extr.xml.txt.2 ), 1 ) 

split.B.path <- project.file.path( "test", "splitB" )
extr.lines.2 <- do.call( 'paste',
  as.list( extractChunksWithCompleteProteinTags(
    split.B.path, split.list.path,
    readChunkFunc=readSplit, 
    appendChunkFunc=appendSplitTillProteinEndTag(
      split.B.path, split.list.path)
    )
  )
)
extr.xml.txt.2 <- c()
extractSingleProteinTags( extr.lines.2,
  each.prot.function=epf )
checkTrue( ! is.null( extr.xml.txt.2 ) )
# Two valid protein in this chunk:
checkEquals( length( extr.xml.txt.2 ), 2 ) 
