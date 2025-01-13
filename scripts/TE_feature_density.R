#######################################################################
#                    TE density in gene features                      #
#---------------------------------------------------------------------#
# Author: Guillermo Peris                                             #
# Version: 13/01/2025                                                 #
# Function: Computes number or density of transposable elements in    #
#           gene features (5'UTR, exons, introns...)                  #
#######################################################################

#-------------------------- Libraries --------------------------------#
library(bedr)
library(reshape2)
library(stringr)
library(dplyr)
library(AnnotationHub)
library(UCSCRepeatMasker)
library(biomartr)
library(BRGenomics)
library(data.table)
library(tools)
source("scripts/auxiliaryFunctions.R")
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#-------------- Parameters (you should modify it! ) ------------------#
#---------------------------------------------------------------------#
#--- Organism: "Homo sapiens", "Mus musculus", "Danio rerio"
organism <- "Felis catus"  

#--- UCSC_TE_annot is record ID for RepeatMasker in Annotation Hub
#--- It depends on organism and assembly
#--- To search your interest ID, run script checkAnnotations.R
UCSC_TE_annot <- "AH98994"   

# Gene annotation release for ensembl
# Find here: https://www.ensembl.org/info/website/archives/assembly.html  
release <- "113" 

#--- interest_TEs: TE families to select. NULL for all TEs. 
interest_TEs <- c("LINE", "LTR", "DNA", "SINE")
# interest_TEs <- NULL

#--- interest_subF: filter TE subfamilies. NULL for all. Overrides interest_TEs.
interest_subF <- NULL
# interest_subF <- c("Alu", "L1")


#--- tag specifies filtering on gene annotation.
#--- Is important to filter spurious annotations
#--- values: NULL (no filter). For other values, run script checkAnnotations.R
tag <- "Ensembl_canonical"  
# tag <- NULL    

#--- minOverlap: default 1E-9. 
# 50 means only consider TEs that overlap at least 50 bp.
# In density analysis this applies to overlapping TE clusters, not individual TEs. 
minOverlap <- 50

#--- downstream: bp after transcription end site.
downstream <- 5000

#--- Relative path to output dir. "" for script folder.
OUTPUT_DIR <- "results/"  

#--- analysis: "density" (default) or "number" (number of overlapping TEs)
analysis <- "number"
#---------------------------------------------------------------------#


#---------------------------------------------------------------------#
#----------- Advanced parameters (you may modify it! ) ---------------#
#---------------------------------------------------------------------#

#--- fileGene: TRUE for file read, FALSE for file import from ensembl.
fileGene <- FALSE    
# if TRUE, gene_annot_file is absolute path to gff annotation
if(fileGene) {
  gene_annot_file  <- "data/Felis_catus.Felis_catus_9.0.113.gtf"  
} 
# Specify gene_annot_format. Use "gff" for ensembl import.
gene_annot_format <- "gff"

#--- fileTE: TRUE for file read, FALSE for file import from repeatMasker.
fileTE <- FALSE      
# if TRUE, TE_annot_file is absolute path to annotation.
if(fileTE) {
  TE_annot_file <- "data/cat_rmsk.txt"
  # Specify TE_annot_format
  TE_annot_format <- "rmsk"
}

#--- feature_types 
feature_types <-  c("gene", "five_prime_utr", "three_prime_utr", 
                 "exon", "intron", "downstream")
#---------------------------------------------------------------------#


#---------------------------------------------------------------------#
#---------------- Read and manage gene annotation --------------------#
#---------------------------------------------------------------------#

if(!fileGene) {
   gene_annot_file <- getGTF(
    db = "ensembl",
    organism = organism,
    remove_annotation_outliers = FALSE,
    path = file.path("/tmp"),
    release = release,
    mute_citation = TRUE
  )
}

#--- Read gene annotation and keep protein coding.
gene_annot <- importGtfTag(gene_annot_file, tag)
gene_annot <- gene_annot[gene_annot$gene_biotype == "protein_coding", ]

#--- Associate gene_id and gene_name for final table. 
gene_names <- data.frame(gene_id = gene_annot$gene_id,
                         gene_name = gene_annot$gene_name)  %>%
  group_by(gene_id) %>% 
  summarise(gene_name = unique(gene_name))
my.names <- gene_names$gene_name
names(my.names) <- gene_names$gene_id
#---------------------------------------------------------------------#


#---------------------------------------------------------------------#
#------------------ Read and manage TE annotation --------------------#
#---------------------------------------------------------------------#
if(!fileTE) {
  ah <- AnnotationHub()
  TE_annot <- ah[[UCSC_TE_annot]]
  TE_annot <- process_ensembl_rm(TE_annot)
} else {
  if(TE_annot_format == "rmsk") {
    TE_annot_file <- process_rmsk_file(TE_annot_file)
    TE_annot_format <- "bed"
  }
  TE_annot <- import(TE_annot_file, format=TE_annot_format)
  TE_annot <- process_file_rm(TE_annot)
}

#--- Filter interest TEs
if (!is.null(interest_subF)) {
  TE_annot <- TE_annot[TE_annot$subfamily %in% interest_subF, ]
} else if(!is.null(interest_TEs)) {
  TE_annot <- TE_annot[TE_annot$family %in% interest_TEs, ]
}  

#--- If analysis == "density" we need to merge all overlapping TEs in order to
#---    avoid counting several times the same bps in feature overlapped region
if(analysis == "density" | analysis != "number") {
  if(analysis != "density") {
    warning("analysis can only be density or number. density is used")
    analysis <- "density"
  }   
  TE_annot_reduced <- reduce(TE_annot, with.revmap = TRUE, ignore.strand = TRUE)
                         
} else {
   TE_annot_reduced <- TE_annot
} 

#-- Sort for downstream analysis
TE_annot_reduced <- sort(TE_annot_reduced, ignore.strand = TRUE)
  

# We save merged TE IDs.
if(analysis == "density") {
revmap <- mcols(TE_annot_reduced)$revmap
mcols(TE_annot_reduced)[, c("name", "TE_name", "subfamily", "family")] <- 
    DataFrame(name = paste(extractList(mcols(TE_annot)$name, revmap), collapse=";"),
              TE_name = paste(extractList(mcols(TE_annot)$TE_name, revmap), collapse=";"),
              subfamily = paste(extractList(mcols(TE_annot)$subfamily, revmap), collapse=";"),
              family = paste(extractList(mcols(TE_annot)$family, revmap), collapse=";"))
mcols(TE_annot_reduced)$revmap <- NULL
}

TE_annot_bed <- bedr.sort.region(granges_to_bed(TE_annot_reduced)) 
TE_annot_bed <- TE_annot_bed %>% 
  dplyr::rename(chrTE=chr, 
         startTE=start,
         endTE = end,
         strandTE = strand)
#---------------------------------------------------------------------#


#--- Begin loop for each feature.
system(paste0("mkdir -p ", OUTPUT_DIR))
for (feature_type in feature_types) {
  outputFile <- paste0(OUTPUT_DIR, "TE_", analysis, "_", feature_type, ".csv")
  
  #--- Introns must be extracted from gene annotation
  if(feature_type == "intron") {
    gene_annot_bed <- intronsFromGFF(gene_annot) 
    
  #--- Downstream regions are computed from 3'UTR.
  } else if (feature_type == "downstream") {
    three_prime <- annotationProcessing(gene_annot, "three_prime_utr")  
    gene_annot_bed <- downstreamProcessing(three_prime, downstream)
    
  #--- Other features.
  } else {
    gene_annot_reduced <- annotationProcessing(gene_annot, feature_type) 
    gene_annot_bed <- bedr.sort.region(granges_to_bed(gene_annot_reduced))
  }
  gene_annot_bed$gene_width = abs(gene_annot_bed$end - gene_annot_bed$start)
  #---------------------------------------------------------------------#
  
  
  #---------------------------------------------------------------------#
  #--------------------- Compute gene/TE overlap -----------------------#
  #---------------------------------------------------------------------#
  
  #--- Compute overlap.
  overlap <- bedr(
    input = list(a = gene_annot_bed, b = TE_annot_bed), 
    method = "intersect", 
    # -wao:  show amount of overlap (even for gene entries though no intersecting TEs)
    # -sorted: assumes inputs are sorted, so that algorithm is faster.
    params = "-wao -sorted")
  
  names(overlap)[ncol(overlap)] <- "width"
  overlap$gene_width <- as.numeric(overlap$gene_width)
  
  #--- Filter TE clusters with overlapping width under minOverlap bp 
  #--- We cannot remove rows because we need exons/introns width even if
  #--- there is no TE overlap or below minOverlap.
  no.overlap.index <- which(as.numeric(overlap$width) <= minOverlap)
  overlap[no.overlap.index, c("chr.bTE", "start.bTE", 
                                     "end.bTE", "strandTE",
                                     "name", "TE_name", 
                                     "subfamily", "family", "width")] <- 
                rep(c(".", "-1", "-1", ".", ".", ".", ".", ".","0"),
                    each = length(no.overlap.index))
  overlap$width <- as.numeric(overlap$width)
  
  
  #--- Compute TE density and TE_count for exons and introns.
  #--- First we group by feature_id (e.g. same exons or introns for a gene).
  #--- Second, we compute width of the overlap.
  #--- Third, we group by gene and sum up all widths.
  if(feature_type %in% c("exon", "intron") ){ 
    overlap_by_gene <- overlap %>% 
      # First we group by TEs in the same feature
      group_by(feature_id)  %>% 
      summarise(name = paste(name, collapse = ";"),
                gene_id =unique(gene_id),
                TE_name = paste(TE_name, collapse = ";"),
                subfamily = paste(subfamily, collapse = ";"),
                family = paste(family, collapse = ";"),
                width = sum(width),
                gene_width = unique(gene_width)) %>%
      group_by(gene_id) %>%   
      summarise(
        name = paste(name, collapse = ";"),
        TE_name = paste(TE_name, collapse = ";"),
        subfamily = paste(subfamily, collapse = ";"),
        family = paste(family, collapse = ";"),
        # We add TEs overlap widths
        width = sum(width),        
        # feature width where all TEs overlap
        feature_width = sum(gene_width)) 
    
    if(analysis == "density") {
      overlap_by_gene <- overlap_by_gene %>% 
        mutate(TE_density = abs(width/feature_width))
    } else {
      overlap_by_gene <- overlap_by_gene %>% 
        mutate(TE_count   = str_count(gsub("-", "_", TE_name), '\\w+'))
    }

    #--- Compute TE density and TE_count for other features.
  } else {
    overlap_by_gene <- overlap %>%              
      group_by(gene_id) %>%   
      summarise(
        name = paste(name, collapse = ";"),
        TE_name = paste(TE_name, collapse = ";"),
        subfamily = paste(subfamily, collapse = ";"),
        family = paste(family, collapse = ";"),
        width = sum(width),        
        feature_width = as.numeric(unique(gene_width))) 
    if(analysis == "density") {
      overlap_by_gene <- overlap_by_gene %>% 
        mutate(TE_density = abs(width/feature_width))
    } else {
      overlap_by_gene <- overlap_by_gene %>% 
        mutate(TE_count   = str_count(gsub("-", "_", TE_name), '\\w+'))
    }
    
  }
  
  
  #---------------------------------------------------------------------#
  #------------------------ Final processing ---------------------------#
  #---------------------------------------------------------------------#
  #--- Add gene names
  overlap_by_gene <- overlap_by_gene %>% 
    mutate(gene_name = my.names[gene_id]) %>% 
    relocate(gene_name, .after=gene_id)
  overlap_by_gene$gene_name <- my.names[overlap_by_gene$gene_id]
  
  #--- Check duplicates
  dups <- sum(duplicated(overlap_by_gene$gene_id))
  if (dups > 0) {
    warning("There are duplicated gene IDs")
    overlap_by_gene[duplicated(overlap_by_gene$gene_id), ]
  }
  
  #--- Check densities
  if(analysis == "density") {
    dens <- nrow(overlap_by_gene %>% filter(TE_density > 1)) 
    if (dens > 0) {
      warning("There are densities over one")
    }
  }
  
  #--- Save csv
  fwrite(overlap_by_gene, file=outputFile, sep="\t", quote = FALSE, row.names = FALSE , dec = ",")
}

if(!fileGene) {
   file.remove(gene_annot_file)
}

outputFile<- paste0(OUTPUT_DIR, "TE_", analysis, "_finalTable.csv")

#--- Construct table structure 
if(analysis == "density") {
  my.variable <- "TE_density"
} else {
  my.variable <- "TE_count"
}
feature_type <- feature_types[1]
file <- paste0(OUTPUT_DIR, "TE_", analysis, "_", feature_type, ".csv")
TE_table <- as.data.frame(fread(file, sep="\t", dec = ",", header = TRUE ))
TE_table_final <- TE_table[, c("gene_id", "gene_name", my.variable)]
colnames(TE_table_final) <- c("gene_id", "gene_name", paste0(feature_type, "_", my.variable))  
#---


for (feature_type in feature_types[-1]) {
  file <- paste0(OUTPUT_DIR, "TE_", analysis, "_", feature_type, ".csv")
  TE_table <- as.data.frame(fread(file, sep="\t", dec = ",", header = TRUE ))
  TE_table <- TE_table[, c("gene_id", my.variable)] 
  colnames(TE_table) <- c("gene_id", paste0(feature_type, "_", my.variable))  
  TE_table_final <- merge(TE_table_final, TE_table, by.x = "gene_id", all.x = TRUE)
}

#--- Change decimals and save csv
TE_table_final <- TE_table_final %>%  mutate_if(is.numeric, round, digits = 4)
fwrite(TE_table_final, file=outputFile, sep="\t", quote = FALSE, 
            row.names = FALSE , dec = ",", na = "NA")
