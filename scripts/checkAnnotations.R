#######################################################################
#                          Check annotations                          #
#---------------------------------------------------------------------#
# Author: Guillermo Peris                                             #
# Version: 08/01/2025                                                 #
#######################################################################

#-------------------------- Parameters--------------------------------#
#--- Organism: "Homo sapiens", "Mus musculus", "Danio rerio"
organism <- "Homo sapiens"
# Gene annotation release for ensembl
# Find here: https://www.ensembl.org/info/website/archives/assembly.html  
release <- "113"

#--------------------- Advanced parameters----------------------------#
# Number of lines to read.
maximalNumberOfLines <- 100000

#-------------------------- Libraries --------------------------------#
library(AnnotationHub)
library(UCSCRepeatMasker)
library(biomartr)
library(stringi)
#---------------------------------------------------------------------#

#---- Check available TE annotations.
ah <- AnnotationHub()
TE <- query(ah, c("RepeatMasker", organism))
#----

#---- Check available filter tags in gene annotation
gene_annot_file <- getGTF(
  db = "ensembl",
  organism = organism,
  remove_annotation_outliers = FALSE,
  path = file.path("/tmp"),
  release = release,
  mute_citation = TRUE
)

# Read file skipping #
lines <- readLines(gene_annot_file, n = maximalNumberOfLines)
lines <- lines[-grep('^#.*', lines)]

extractTags <- function(line) {
  m.cols <- unlist(strsplit(line, "\t"))[9]
  m.cols <- unlist(strsplit(m.cols, ";"))
  m.cols <- stri_replace_all_charclass(m.cols, "\\p{WHITE_SPACE}", "")
  m.cols <- gsub('\"', "", m.cols, fixed = TRUE)
  index.tags <- startsWith(m.cols, "tag")
  tags <- gsub("tag", "", m.cols[index.tags])
  tags
}

tags <- lapply(lines[1:10], function(x) extractTags(x))
uniqueTags <- unique(unlist(tags))

#------- Output
cat("RepeatMasker annotation: ")
TE
cat(paste0("Feature filters: NULL ", uniqueTags))
