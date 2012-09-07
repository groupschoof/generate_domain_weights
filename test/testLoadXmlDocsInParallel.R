library(RUnit)
library(tools)
library(parallel)
library(RCurl)
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
src.project.file('src','parseUniprotInterProMatch.R')

uni.accs <- as.character(read.table(project.file.path("test",
      "test_1000_sprot_accessions.txt"))$V1)

uris <- sapply(uni.accs, uniprotInterProMatchUrl)
print(head(uris))
res <- try(getURIAsynchronous(uris))

checkTrue( ! identical(class(res), 'try-error'))
