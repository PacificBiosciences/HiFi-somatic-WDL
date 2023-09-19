# A set of copy number related helper functions
require(tidyverse)

#' @param segs A segment file with chromosome, Start and End columns.
#' @param cytoband_length Length of cytoband. Can be generated using load_cytoband() function.
#' @param special_arm_offset If supplied, offset the missing p arm at 13, 14, 15, 21 and 22.
#' @param sex If TRUE (default), remove X and Y chromosome from segments
#' @return Original with segments Start and End offset to absolute position for plotting
#' @export
set_absolute_pos <- function(segs, cytoband_length, sex=TRUE){
  # Sequenza uses start.pos and end.pos
  colnames(segs)[colnames(segs) == "start.pos"] <- "Start"
  colnames(segs)[colnames(segs) == "end.pos"] <- "End"
  # Disregards uppercase
  colnames(segs)[colnames(segs) == "start"] <- "Start"
  colnames(segs)[colnames(segs) == "end"] <- "End"
  colnames(segs)[colnames(segs) == "Chromosome"] <- "chromosome"
  # Add cumulative length and get absolute position for segments
  # Requires chromosome, Start and End (Case sensitive)
  if (sex){
    chrom <- c(seq(1,22), "X", "Y")
  } else {
    chrom <- seq(1,22)
  }
  cum_length <- 0
  for (chr in chrom){
    # Sometimes no segments at certain chromosome.
    if(nrow(segs[which(segs$chromosome == chr), ]) != 0){
      segs[which(segs$chromosome == chr), 'Start'] <- segs[which(segs$chromosome == chr), 'Start'] + cum_length
      segs[which(segs$chromosome == chr), 'End'] <- segs[which(segs$chromosome == chr), 'End'] + cum_length
      cum_length <- cum_length + cytoband_length[which(cytoband_length$chromosome == chr), 'length', drop=T]
    } else{
      cum_length <- cum_length + cytoband_length[which(cytoband_length$chromosome == chr), 'length', drop=T]
    }
  }

  if (sex){
    return(segs)
  } else{
    return(filter(segs, chromosome %in% chrom))
  }
}

#' @param cytoband_ifile Cytoband reference file (See default file 
#' at gis-scripts/references/cytoBand.txt)
#' @param offset_arm By default this will not remove p arm of acrocentric chromosome.
#' Set to yes to remove.
#' @param sex By default this loads the sex cytobands. If set to FALSE, remove sex chromosome
#' @return a cytoband list with various useful dataframe 
#' for plotting and annotation
#' @export
load_cytoband <- function(cytoband_ifile="/home/chuakp/gis-scripts/references/cytoBand.txt",
                          sex=TRUE){
			    
  # Read in cytoband for use with segments annotation etc
  cytoband <- read_tsv(cytoband_ifile, col_names = F, col_types = cols())

  colnames(cytoband) = c("chromosome", "Start", "End", "arm", "Misc")
  cytoband <- cytoband %>% mutate(chromosome = gsub("chr(.*)", "\\1", chromosome))
  
   # Remove sex chromosome if specified
  if(sex==FALSE){
    cytoband <- cytoband %>% 
      filter(!chromosome %in% c("X", "Y"))
  }


  cytoband_length <- cytoband %>% group_by(chromosome) %>% summarise(length = sum(End-Start))
  # Get arm position
  cytoband_arm <- cytoband %>% filter(Misc=="acen")
  cytoband_arm <- set_absolute_pos(cytoband_arm, cytoband_length)

  # Find chromosome position
  cytoband_pos <- cytoband_length
  cytoband_pos <- cytoband_pos %>% arrange(factor(cytoband_pos$chromosome, levels = c(seq(1,22), "X", "Y")))
  # Get absolute position of chrom end
  cytoband_pos <- cytoband_pos %>% mutate(Start=cumsum(as.numeric(length)) - length,
                                          End=cumsum(as.numeric(length)))
  
  # Set alternate color
  cytoband_pos <- cytoband_pos %>% mutate(
    color=if_else(
      row_number()%%2 == 0, "grey", "white"
    )
  )
  
  # Set label for p and q
  cytoband_arm <- cytoband_arm %>% group_by(chromosome) %>% summarise(Start=min(Start), End=max(End))
  cytoband_arm <- cytoband_arm %>% mutate(start_label = paste0(chromosome, "p"),
        				  end_label = paste0(chromosome, "q"))

  # Set label for mid point of p and q, set color for use with mixed labelling later
  cytoband_arm <- cytoband_arm %>% mutate(midpoint=(Start + End)/2)
  cytoband_arm$color <- "black"
  
  cytoband_arm <- cytoband_arm %>% 
    inner_join(cytoband_pos %>% 
                 dplyr::rename(Start_p="Start", End_q="End") %>% 
                 select(-color), by="chromosome")
  
  cytoband_arm$End_p <- cytoband_arm$Start
  cytoband_arm$Start_q <- cytoband_arm$End

  return(list(cytoband_coord=cytoband,
	cytoband_length=cytoband_length, 
	cytoband_arm=cytoband_arm,
	cytoband_pos=cytoband_pos))
}


#' @param segs A segment file with Start, End, and CNt (Total copy number)
#' @param region_toRemove A dataframe with chromosome, Start, End and length
#' containing the region that you want to remove from plots. Normally this is the 
#' gvar (variable length) and stalk (p arm for 13, 14, 15, 21 and 22) region in 
#' cytoband file. Example of region_toRemove can be generated using cytoband list from
#' load_cytoband():
#' cytoband <- load_cytoband()
#' cytoband_toRemove <- cytoband$cytoband_coord %>% 
#'   filter(Misc %in% c("gvar", "stalk")) %>% 
#'   mutate(length = End - Start) %>% 
#'   arrange(factor(chromosome, levels = c(1:22,"X","Y")))
#' cytoband_toRemove <- set_absolute_pos(cytoband_toRemove,cytoband$cytoband_length,
#'                                       sex = FALSE)
#' @param chromPos This is used to offset chromosome label. It's necessary because
#' chromosome position spans over whole the region to offset, hence the condition
#' where start of segments is after the end of region to remove is FALSE.
#' @return # Segments with unwanted region removed. Region to remove should be in 
#' absolute position already. The code works by going throw
#' each row of the regions to remove and minus the length of the region
#' for all segments with "Start" more than or equal to the end of the region.
#' CAUTION: This code assumes there's no segments overlapping the region you want to remove!
#' You may want to use bedtools intersect first before using this if there are regions
#' overlapping.
#' @export
remove_region <- function(segs, region_toRemove){
  # Sort region_toRemove by start point
  region_toRemove <- arrange(region_toRemove, Start)
  # print(paste0("Total number of segments: ", nrow(segs)))
  for (i in 1:nrow(region_toRemove)){
    regStart <-  region_toRemove[i, "Start", drop=TRUE]
    regEnd <-  region_toRemove[i, "End", drop=TRUE]
    regLength <- region_toRemove[i, "length", drop=TRUE]
    # print(paste0("Removing region:", regStart, "-", regEnd))
    # Remove segments within the region
    segs <- segs %>% filter(!(Start >= regStart & End <= regEnd))
    # Handles different types of segments 
    overlapEnd <- segs[which(segs$Start <= regStart & segs$End < regEnd & segs$End > regStart),]
    overlapStart <- segs[which(segs$Start > regStart & segs$Start < regEnd & segs$End >= regEnd),]
    coverRegion <- segs[which(segs$Start <= regStart & segs$End >= regEnd),]
    beforeRegion <- segs[which(segs$Start <= regStart & segs$End <= regStart),]
    afterRegion <- segs[which(segs$Start >= regEnd & segs$End >= regEnd),]
    # Sanity check
    totalRow <- nrow(overlapEnd) + nrow(overlapStart) + 
      nrow(coverRegion) + nrow(beforeRegion) + nrow(afterRegion)
    if (nrow(segs) != totalRow){
      print(paste0("Number of rows processed: ", totalRow,
                   " is different from original number of rows: ", nrow(segs)))
    }
    
    if(nrow(overlapEnd) != 0){
      overlapEnd$End <- regStart
      overlapEnd[, c("Start", "End")] <- overlapEnd[, c("Start", "End")] - regLength
    }
    
    if(nrow(overlapStart) != 0){
      overlapStart$Start <- regEnd
      overlapStart[, c("Start", "End")] <- overlapStart[, c("Start", "End")] - regLength
    }
    
    if(nrow(coverRegion) != 0){
      coverRegion$End <- coverRegion$End - regLength
    }
    
    if(nrow(afterRegion) != 0){
      afterRegion[, c("Start", "End")] <- afterRegion[, c("Start", "End")] - regLength
    }
    
    # Merge back the processed segments
    segs <- bind_rows(overlapStart, overlapEnd, coverRegion, beforeRegion, afterRegion)
    
    # Need to offset the region to remove itself
    region_toRemove[region_toRemove$Start >= regEnd, c("Start", "End")] <-
      region_toRemove[region_toRemove$Start >= regEnd, c("Start", "End")] - regLength
    
    # For handling midpoint labelling
    if (any(colnames(segs) == "midpoint")){
      segs[segs$midpoint >= regEnd, c("midpoint")] <-
        segs[segs$midpoint >= regEnd, c("midpoint")] - regLength
    }
  }
  return(segs)
}

#' @param segs A segment file with Start, End, and CNt (Total copy number)
#' @param toMed The column that you want to calculate median of and center based on. Don't
#' need to quote the argument. By default this is CNt. You can use logR
#' @return segment file with weighted mean ploidy and median ploidy calculated
#' as well as an adjustedCN column (Total CN minus median ploidy)
#' @export
calc_medianPloidy <- function(segs, toMed = CNt){
  # Get median ploidy adjusted copy number. Will return the segments with 
  # median ploidy, adjustedCN (total CN minus median ploidy) and weighted mean
  # ploidy (This is the same as the mean ploidy in Sequenza confint output)
  
  # Calculate median ploidy on per base level. Basically how many
  # bases have X copy number. Then look for the median copy number
  CNt <- enquo(toMed)
  colname_toUse <- as.character(enexpr(toMed))
  median_colName <- quo_name(paste0("median", colname_toUse))
  median_adjustName <- quo_name(paste0("medianCentered_", colname_toUse))
  mode_colName <- quo_name(paste0("mode", colname_toUse))
  
  medPloid <- segs %>% mutate(seg_length = End - Start) %>% 
    group_by(Tumor_Sample_Barcode, !!CNt) %>%
    summarise(totalNum = sum(as.numeric(seg_length))) %>% 
    arrange(!!CNt) %>% 
    mutate(cumsum_length = cumsum(as.numeric(totalNum)), 
           larger_thanMed = (cumsum_length > (sum(as.numeric(totalNum)) + 1)/2)) %>%
    filter(larger_thanMed) %>% filter(row_number() == 1) %>% select(Tumor_Sample_Barcode, !!CNt) %>% 
    ungroup %>%  rename(`:=`(!!median_colName, !!CNt))
  
  segs <- inner_join(segs, medPloid, by = "Tumor_Sample_Barcode")
  colnames(segs)[colnames(segs) == "start.pos"] <- "Start"
  colnames(segs)[colnames(segs) == "end.pos"] <- "End"
  segs <- segs %>% group_by(Tumor_Sample_Barcode) %>%
    mutate(seg_length = End - Start, 
           !!median_adjustName := !!CNt - !!sym(median_colName), 
           weighted_meanPloidy = sum(as.numeric(seg_length) * !!CNt)/sum(as.numeric(seg_length))) %>%
  ungroup()
  
  return(segs)
}

#' @param segs A segment file with Start, End, and CNt (Total copy number)
#' @param cytoband cytoband coordinates generated using load_cytoband
#' @param highlight_seg A vector of Chromosome, Start and End to indicate the region
#' to plot. By default whole genomes are plotted
#' @param highlight_pos A vector of 3, chromosome, Start and End to highlight as
#' a line on the plot
#' @return segment file with weighted mean ploidy and median ploidy calculated
#' as well as an adjustedCN column (Total CN minus median ploidy)
#' @export
plot_all_segments <- function(segs, cytoband, 
			      highlight_seg = NULL, highlight_pos = NULL){
  # Takes in a segments tibble with median adjusted CN and plot segments for all
  # Tumor_Sample_Barcode and cytoband info (Generated with helper function
  # load_cytoband above)
  # If you want to plot just specific place, supply a list of vector of
  # chromosome, start and end (in this order), otherwise all segments are
  # plotted
  
  cmap <- rev(c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac'))
  breaklist <- seq(-4,4,1)
  
  # Only plot the region you want
  if(!is.null(highlight_seg)){
    highlight <- tibble(chromosome=highlight_seg[1], Start=highlight_seg[2], End=highlight_seg[3])
    highlight <- set_absolute_pos(highlight, cytoband$cytoband_length)
    xscale_lim <- c(highlight$Start, highlight$End)
  } else{
    xscale_lim <- c(0, NA)
  }
  
  # Plot segments
  p <- ggplot(data=segs) +
    geom_segment(aes(x=Start, xend=End, 
                     y=Tumor_Sample_Barcode, yend=Tumor_Sample_Barcode, 
                     color=factor(adjustedCN)), size=3) +
    geom_vline(xintercept=cytoband$cytoband_pos$pos) +
    geom_vline(xintercept = cytoband$cytoband_arm$midpoint, alpha=0.5, linetype=2) +
    scale_x_continuous(breaks=cytoband$cytoband_arm$midpoint, 
                       labels = cytoband$cytoband_arm$chromosome,
                       expand=c(0, 0), limits = xscale_lim) +
    scale_color_manual(breaks=breaklist, values = cmap, limits = breaklist) +
    theme_bw() +
    labs(x = "Chromosome", y="Tumor Sample Barcode") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
  if (!is.null(highlight_pos)){
    p <- p + 
    geom_vline(xintercept = highlight_pos, alpha=0.4, linetype="dotted")
  }
  return(p)
}

# Helper function for get_gene_in_segments
splitColumnByOverlap <- function(query, subject, column="ENTREZID", ...){
  require(GenomicRanges)
  olaps <- findOverlaps(query, subject, ...)
  # Keep only if overlap more than 80% of the genes
  overlaps <- pintersect(query[queryHits(olaps)], subject[subjectHits(olaps)])
  percentageOverlaps <- width(overlaps) / width(query[queryHits(olaps)])
  olaps <- olaps[percentageOverlaps >= 0.8]
  f1 <- factor(subjectHits(olaps),
               levels=seq_len(subjectLength(olaps)))
  splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
}

#' This function takes a segments and annotate each segments with coding genes from
#' RefSeq (default). If given a list of drivers, will only subset to drivers within
#' the RefSeq bed file.
#' @param query_df A tibble with chromosome, start and end as the first three columns.
#' Extra columns will be kept.
#' @param filterDrivers If given a list, will filter the refseq bed to only the list of genes.
#' Only the genes in the filtered bed will be annotated on the tibble.
#' @param ref_seq_bed Any bed file with first three columns as chr, start and end. Last 
#' column should be gene names.
#' @return Original dataframe with genes within the region annotated.
#' @export
get_gene_in_segments <- function(query_df, 
                                 filterDrivers=NULL, 
                                 ref_seq_bed="~/gis-scripts/references/ref_seq_hg19_clean_uniqGene_sorted.bed"){
  require(tidyverse)
  # Rename chromosome
  colnames(query_df)[colnames(query_df) == "Chromosome"] <- "chromosome"
  query_grange <- GenomicRanges::makeGRangesFromDataFrame(query_df,
                                                          keep.extra.columns = TRUE)
  # Load ref_seq 
  refseq <- read.delim(ref_seq_bed, header = FALSE)
  colnames(refseq) <- c("chromosome", "Start", "End", "RefSeq_Genes")
  if (!is.null(filterDrivers)){
    refseq <- refseq[refseq$RefSeq_Genes %in% filterDrivers, ]
  }
  refseq <- GenomicRanges::makeGRangesFromDataFrame(refseq, keep.extra.columns = TRUE)
  annotated_genes <- splitColumnByOverlap(refseq, query_grange, "RefSeq_Genes")
  # Merge annotated gene back to original df
  query_df <- query_df %>% mutate(rnum=row_number()) %>% 
    rowwise() %>% mutate(gene=paste(unique(annotated_genes[[rnum]]), collapse =",")) %>% 
    ungroup() %>% select(-rnum)
  return(query_df)
}

# Below are a bunch of functions used to parse GISTIC file and plots GISTIC q-value plots
# Function to parse all_lesions file
parse_lesion_file <- function(lesion_file){
  
  all_lesions <- read_tsv(lesion_file, col_types=cols())
  
  # Remove the useless row of CN values
  all_lesions <- all_lesions %>% filter(!grepl("CN values", `Unique Name`))
  
  # Remove probes info on peal limit
  all_lesions <- all_lesions %>% 
    mutate(`Peak Limits` = gsub("\\(.*\\)", "", `Peak Limits`))
  
  # Split peak limits
  all_lesions <- all_lesions %>%
    separate(`Peak Limits`, sep=":|-", into=c("chromosome", "Start", "End"))
  
  # Remove chr
  all_lesions <- all_lesions %>% 
    mutate(chromosome = gsub("chr", "", chromosome))
  
  # Save only the useful column
  all_lesions <- all_lesions %>%
    select(chromosome, Start, End,
           `Unique Name`, `Descriptor`,
           `q values`) %>% 
    mutate_at(c("chromosome", "Start", "End"), as.numeric)
}

set_zero_centromere <- function(scores_file, forPlot, cytoband_arm){
  # browser()
  # Tricky, if there's a broad peak across centromere, we cannot just set 0 for centromere start and end.
  # Basically, we need to break the peak that goes across centromere into three parts and set the centromere
  # to 0.
  for (chr in unique(scores_file$chromosome)){
    cytoband_chr <- cytoband_arm %>% filter(chromosome == chr)
    segments_to_break <- scores_file %>% 
      filter(chromosome == chr & Start <= cytoband_chr$Start & End >= cytoband_chr$End)
    if(nrow(segments_to_break) != 0){
      seg_append <- tibble(toplot=c(rep(segments_to_break$toplot, 2), rep(0, 2)), 
                           position=c(cytoband_chr$Start, cytoband_chr$End, cytoband_chr$Start+1, cytoband_chr$End-1))
    } else{
      seg_append <- tibble(toplot=rep(0, 2),
                           position=c(cytoband_chr$Start, cytoband_chr$End))
    }
    forPlot <- forPlot %>% bind_rows(seg_append)
  }
  return(forPlot)
}

parse_peak_genes <- function(sig_gene, drivers_list){
  # Parse annoying gistic peak gene file
  peak_info <- t(sig_gene[1:3, ])
  colnames(peak_info) <- peak_info[1, 1:3]
  # Unselect first and last row
  peak_info <- peak_info[-c(1,nrow(peak_info)), ,drop=FALSE]
  rowlabels_cytoband <- rownames(peak_info)
  peak_info <- as_tibble(peak_info) %>% mutate(cytoband = rowlabels_cytoband)
  
  collapse_driver_gene <- function(cytoband){
    genes <- sig_gene[4:(nrow(sig_gene) - 1), cytoband]
    # Remove NA
    genes <- genes[!is.na(genes),,drop=TRUE]
    # Match CGC
    genes <- genes[match(drivers_list, genes)]
    genes <- genes[!is.na(genes)]
    return(paste(genes, collapse = ","))
  }
  
  collapse_gene <- function(cytoband){
    genes <- sig_gene[4:(nrow(sig_gene) - 1), cytoband]
    # Remove NA
    genes <- genes[!is.na(genes),,drop=TRUE]
    # Match CGC
    # genes <- genes[match(drivers_list, genes)]
    # genes <- genes[!is.na(genes)]
    return(paste(genes, collapse = ","))
  }
  
  peak_info <- peak_info %>% rowwise() %>% 
    mutate(driver_genes = collapse_driver_gene(cytoband),
           all_genes_gistic = collapse_gene(cytoband)) %>% ungroup
  
  peak_info <- peak_info %>% separate(`wide peak boundaries`, 
                                      into=c("chromosome", "Start", "End"), 
                                      sep = ":|-")
  peak_info$chromosome <- gsub("chr", "", peak_info$chromosome)
  # Convert columns to numeric
  peak_info <- peak_info %>% mutate_at(c("q value", "residual q value", 
                                         "chromosome", "Start", "End"), as.numeric)
  # Get -log10 of q-value
  peak_info <- peak_info %>% mutate(neg_log10_qval = -log10(`q value`))
  
  return(peak_info)
}
