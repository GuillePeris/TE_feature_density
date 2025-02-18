#######################################################################
#            Auxiliary functions for TE_in_gene_density.R             #
#---------------------------------------------------------------------#
# Author: Guillermo Peris                                             #
# Version: 13/01/2025                                                 #
#######################################################################

#-----------------------------------------------#
#--- granges_to_bed: converts GRanges object ---# 
#--- to bed format keeping metacolumns.      ---#
#-----------------------------------------------#
granges_to_bed <- function(gr) {
  df <- data.frame(chr = as.character(seqnames(gr)),
                   start = start(gr),
                   end = end(gr), 
                   strand = strand(gr))
  my.cols <- colnames(mcols(gr))
  if(length(my.cols) > 0) {
    df <- cbind(df, mcols(gr))
  }
  return(df)
}


#-----------------------------------------------#
#--- collapseByGeneIDandFeatureID: reduces   ---# 
#--- GRanges by geneID and featureID and     ---#
#--- collapses all merge features.           ---#
#-----------------------------------------------#
collapseByGeneIDandFeatureID <- function(gr, feature) {
  # Split the GRanges object by gene_id
  gr_list <- split(gr, gr$gene_id)
  
  # Reduce the ranges and concatenate the ids
  reduced_list <- lapply(names(gr_list), function(gene_id) {
    x <- gr_list[[gene_id]]
    reduced <- reduce(x)
    overlaps <- findOverlaps(reduced, x)
    id_list <- vector("list", length(reduced))
    for (i in seq_along(reduced)) {
      hits <- subjectHits(overlaps[queryHits(overlaps) == i])
      id_list[[i]] <- paste(unique(mcols(x)[[feature]][hits]), collapse = ",")
    }
    mcols(reduced)$feature_id <- unlist(id_list)
    mcols(reduced)$gene_id <- gene_id
    return(reduced)
  })
  
  # Combine the reduced GRanges objects back into one GRanges object
  reduced_gr <- unlist(GRangesList(reduced_list))
  return(reduced_gr)
}


#-----------------------------------------------#
#--- downstreamProcessing: creates           ---# 
#--- downstream feature from 3'UTR.          ---#
#-----------------------------------------------#
downstreamProcessing <- function(three_prime, downstream) {
  start(three_prime[strand(three_prime) == "+"])   <- 
    end(three_prime[strand(three_prime) == "+"])   +  1
  end(three_prime[strand(three_prime) == "+"])   <- 
    end(three_prime[strand(three_prime) == "+"])   +  downstream
  end(three_prime[strand(three_prime) == "-"])   <- 
    start(three_prime[strand(three_prime) == "-"])   -  1
  start(three_prime[strand(three_prime) == "-"])   <- 
    start(three_prime[strand(three_prime) == "-"])   -  downstream
  gene_annot_bed <- granges_to_bed(three_prime)
  cols <- colnames(gene_annot_bed)
  gene_annot_bed <- bedtoolsr::bt.sort(gene_annot_bed) 
  colnames(gene_annot_bed) <- cols
  return(gene_annot_bed)
}

#-----------------------------------------------#
#--- annotationProcessing: process GRanges   ---# 
#--- depending on feature type.              ---#
#-----------------------------------------------#
annotationProcessing <- function(gene_annot, feature_type) {
  #--- Filter by gene type
  gene_annot_type <- gene_annot[gene_annot$type == feature_type, ]
  #--- Remove 1bp features
  gene_annot_type <- gene_annot_type[width(gene_annot_type) > 1, ]
  
  # For 5'UTR/3'UTR keep first/last exon
  if(feature_type %in% "five_prime_utr") {
    five_UTR_fwd <- gene_annot_type[which(strand(gene_annot_type)=="+")]
    five_UTR_rev <- gene_annot_type[which(strand(gene_annot_type)=="-")]
    five_UTR_last_fwd <- as.data.frame(five_UTR_fwd) %>% group_by(gene_id) %>% dplyr::slice(which.min(start))
    five_UTR_last_rev <- as.data.frame(five_UTR_rev) %>% group_by(gene_id) %>% dplyr::slice(which.max(end))
    gene_annot_type <- sort(makeGRangesFromDataFrame(rbind(five_UTR_last_fwd, five_UTR_last_rev), keep.extra.columns=TRUE))
  } else if (feature_type %in% "three_prime_utr") {
    three_UTR_fwd <- gene_annot_type[which(strand(gene_annot_type)=="+")]
    three_UTR_rev <- gene_annot_type[which(strand(gene_annot_type)=="-")]
    three_UTR_last_fwd <- as.data.frame(three_UTR_fwd) %>% group_by(gene_id) %>% dplyr::slice(which.max(end))
    three_UTR_last_rev <- as.data.frame(three_UTR_rev) %>% group_by(gene_id) %>% dplyr::slice(which.min(start))
    gene_annot_type <- sort(makeGRangesFromDataFrame(rbind(three_UTR_last_fwd, three_UTR_last_rev), keep.extra.columns=TRUE))
  }
  
  # Defining feature_id. In downstream pipeline features will be grouped
  # first by feature_id (e.g., exons in the same transcript) and then by gene_id.
  if(feature_type %in% c("five_prime_utr", "three_prime_utr", "exon")) {
    feature <- "exon_id"
  } else if (feature_type == "gene") {
    feature <- "gene_id"
  } 
  
  gene_annot_reduced <- collapseByGeneIDandFeatureID(gene_annot_type, feature)
  gene_annot_reduced <- sort(gene_annot_reduced, ignore.strand=TRUE)
  return(gene_annot_reduced)
}

#-----------------------------------------------#
#--- intronsFromGFF: annotates introns       ---#
#-----------------------------------------------#
intronsFromGFF <- function(gff) { 
  #--- Get genes
  genes <- gff[gff$type == "gene", ]
  genes <- sort(genes, ignore.strand=TRUE)
  genes_bed <- granges_to_bed(genes)
  my.colnames <- c("chr", "start", "end", "gene_id")
  genes_bed <- genes_bed[, my.colnames]
  genes_bed <- bedtoolsr::bt.sort(genes_bed)
  colnames(genes_bed) <- my.colnames
  
  
  #--- Get protein_coding exons and reduce by gene
  exons <- gff[gff$type == "exon", ] 
  exons_reduced <- reduce(exons)
  
  #--- Remove 1bp features
  exons_reduced <- exons_reduced[width(exons_reduced) > 1, ]
  exons_bed <- granges_to_bed(exons_reduced)
  exons.cols <- colnames(exons_bed)
  exons_bed <- bedtoolsr::bt.sort(exons_bed)
  colnames(exons_bed) <- exons.cols
  
  #--- Getting introns by substracting exons from genes
  complement <- bedtoolsr::bt.subtract(genes_bed, exons_bed)

  # complement <- bedr.sort.region(complement)
  complement <- bedtoolsr::bt.sort(complement) 
  colnames(complement) <- colnames(genes_bed)

  #--- We need to assign an arbitrary name to intron.
  complement$feature_id <- paste0("intron_", seq(1:nrow(complement)))
  return(complement) 
}

#-----------------------------------------------#
#--- process_file_rm: processes RepeatMasker ---#
#--- from user file.                         ---#
#-----------------------------------------------#
process_file_rm <- function(TE_annot) {
  #--- Remove non-canonical chromosomes and uninteresting columns.
  TE_annot <- tidyChromosomes(TE_annot)
  mcols(TE_annot)[, c("score", "itemRgb", "thick")] <- NULL
  
  #--- Split name in TE_name:subfamily:family.
  newColNames <- c("chr", "start", "end", "TE_name", "subfamily", "family", "milliDV", "TE_strand")
  newCols <-  reshape2::colsplit(str_replace_all(TE_annot$name, ":", "|"), "\\|", newColNames)
  mcols(TE_annot)[, c("TE_name", "subfamily", "family")] <- newCols[, c("TE_name", "subfamily", "family")]
  
  # Sort GRange by coordinates (ignore strand)
  TE_annot <- sort(TE_annot, ignore.strand=TRUE)
  
  return(TE_annot)
}
  
#-----------------------------------------------#
#--- process_file_rm: processes RepeatMasker ---#
#--- from ensembl download.                  ---#
#-----------------------------------------------#
process_ensembl_rm <- function(TE_annot) {
  
  # Remove non standard chromosomes.
  TE_annot <- tidyChromosomes(TE_annot)
  
  # Remove uninteresting repeats (same filter as SQuIRE).
  TE_annot <- TE_annot[!(TE_annot$repClass %in% 
                           c("Simple_repeat", "Satellite", "Unknown", "Low_complexity")), ]
  
  # Sort GRanges by coordinates (ignore strand)
  TE_annot <- sort(TE_annot, ignore.strand=TRUE)
  
  # Add metacols
  mcols(TE_annot)[, c("name", "TE_name", "subfamily", "family")] <- 
    data.frame(name = paste0(seqnames(TE_annot), "|",
                             start(TE_annot) + 1, "|",
                             end(TE_annot), "|",
                             TE_annot$repName, ":",
                             TE_annot$repFamily, ":",
                             TE_annot$repClass, "|",
                             TE_annot$milliDiv, "|",
                             strand(TE_annot)),
               TE_name = TE_annot$repName,
               subfamily = TE_annot$repFamily,
               family = TE_annot$repClass)
  mcols(TE_annot)[, 1:11] <- NULL
  
  return(TE_annot)
}

################################################
#-----------------------------------------------#
#--- importGtfTag: extracts lines from GTF   ---#
#--- file with arguments tag ("basic",       ---#
#--- "Ensembl_canonical", "MANE_select"...)  ---#
#--- and keeps all features "gene".          ---#
#--- This function tries to solve the        ---#
#--- problem of multiple tags in annotation. ---#
#-----------------------------------------------#
importGtfTag <- function(GTF_file, tag) {
  
  # Read file skipping #
  lines <- readLines(GTF_file)
  lines <- lines[-grep('^#.*', lines)]
  
  # Filter by feature
  lines.df <- fread(text = lines)
  lines.genes <- lines.df %>% filter(V3 == "gene")
  
  # Filter by tag
  lines.tag <- lines.df[grep(tag, lines.df$V9), ]
  
  # Join genes and other features
  lines.write <- rbind(lines.genes, lines.tag)
  
  # Convert to line and write to temp file
  tmpFile <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".gff")
  fwrite(lines.write, file=tmpFile, sep="\t", quote = FALSE, 
              col.names = FALSE, row.names = FALSE)
  
  # Read tmpFile as GTF
  tmp_annot <- import(tmpFile, format="gff")
  seqlevelsStyle(tmp_annot) <- "UCSC"
  tmp_annot <- tidyChromosomes(tmp_annot)
  tmp_annot <- sort(tmp_annot, ignore.strand=TRUE)
  
  # Clear tmpFile
  out <- file.remove(tmpFile)
  
  return(tmp_annot)
}

#-----------------------------------------------#
#--- process_rmsk_file: processes            ---#
#--- RepeatMasker file downloaded from UCSC  ---# 
#--- with rmsk format                        ---#
#-----------------------------------------------#
process_rmsk_file <- function(TE_annot_file){
  TE_table <- fread(TE_annot_file, sep="\t", dec = ",", header = TRUE )
  
  # Remove uninteresting repeats (same filter as SQuIRE).
  TE_table <- TE_table[!(TE_table$repClass %in% 
                           c("Simple_repeat", "Satellite", "Unknown", "Low_complexity")), ]
  
  
  TE_new_table <- data.frame(chr = TE_table$genoName,
                             start = TE_table$genoStart,
                             end = TE_table$genoEnd,
                             name = paste0(TE_table$genoName,  "|",
                                           TE_table$genoStart, "|",
                                           TE_table$genoEnd,   "|",
                                           TE_table$repName,   ":",
                                           TE_table$repFamily,  ":",
                                           TE_table$repClass, "|",
                                           TE_table$milliDiv,  "|",
                                           TE_table$strand),
                             score = TE_table$swScore,
                             strand = TE_table$strand   
)
  
  # outputFile
  outputFile <- paste0(file_path_sans_ext(TE_annot_file), ".bed")
  fwrite(TE_new_table, file=outputFile, sep="\t", quote = FALSE, 
              col.names = FALSE, row.names = FALSE , dec = ",")
  return(outputFile)
}
