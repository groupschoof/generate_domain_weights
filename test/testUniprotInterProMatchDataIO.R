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

# Test extractCompleteProteinTags
print("Testing extractCompleteProteinTags(...)")
uni.ipr.scn.path <- project.file.path("test",
    "testUniprotInterProCompleteMatchFile.xml")
extr.lines <- extractCompleteProteinTags( uni.ipr.scn.path, 1, 25 )
print( extr.lines )
checkEquals( extr.lines, scan(uni.ipr.scn.path, sep="\n", nlines=41, what=character()) )

