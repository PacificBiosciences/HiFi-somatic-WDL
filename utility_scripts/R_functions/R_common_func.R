## A set of common functions I use across other scripts

# Colormap from stackexchange, good for class palette
c25 <- c("dodgerblue2","#E31A1C", # red
         "green4",
         "#6A3D9A", # purple
         "#FF7F00", # orange
         "black","gold1",
         "skyblue2","#FB9A99", # lt pink
         "palegreen2",
         "#CAB2D6", # lt purple
         "#FDBF6F", # lt orange
         "gray70", "khaki2",
         "maroon","orchid1","deeppink1","blue1","steelblue4",
         "darkturquoise","green1","yellow4","yellow3",
         "darkorange4","brown")

# Coding definition (non-syn including silent) for TCGA MAF format
coding_def <- c("De_novo_Start_InFrame", "De_novo_Start_OutOfFrame",
                "Start_Codon_Del", "Start_Codon_SNP", "Frame_Shift_Del",
                "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
                "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
                "Silent", "Splice_Site", "Start_Codon_Ins",
                "Translation_Start_Site", "Stop_Codon_Del")

#' Function to quickly visualize color palette in bar chart
#' Input is a vector of color.
bar_color <- function(color_vec){
  barplot(rep(1, length(color_vec)), 
  col=color_vec, 
  names.arg=1:length(color_vec))
}


calc_MATH <- function(vaf_vec){
  # Calculate MATH Score
  # vaf_vec is a vector of the VAF
  # Get MAF for all sample in a maf file by doing:
  # math_score <- maf %>% group_by(Tumor_Sample_Barcode) %>% summarise(MATH = calc_MATH(VAF))
  vaf_mad <- mad(vaf_vec)
  math_score <- 100 * vaf_mad / median(vaf_vec)
}


get_subset <- function(maf_file, output_prefix, samples){
  # This function takes in a MAF file, output prefix, and a list of samples, and output
  # the VAF matrix of the mutations shared across the samples and those that are missing
  # in at least one sectors.
  # Example usage: get_subset("tmp.maf", "noPON_filter", c("A126", "A552"))
  
  all_mutations <- read_tsv(maf_file, 
                            col_types = cols(".default" = "c"))
  all_mutations <- all_mutations %>% separate(allelic_depth, into=c("ref_freq", "alt_freq"),
                                              sep=",")
  all_mutations <- all_mutations %>% mutate_at(c("ref_freq", "alt_freq"), as.numeric)
  
  all_mutations <- all_mutations %>% mutate(VAF=alt_freq/(alt_freq+ref_freq),
                                            read_depth=(alt_freq+ref_freq))
  
  # Filter low read depth
  all_mutations <- all_mutations %>% filter(read_depth > 10)
  
  subset <- all_mutations %>% filter(Tumor_Sample_Barcode %in% samples) %>% 
    select(Tumor_Sample_Barcode, Hugo_Symbol, Chromosome, Start_position, End_position, 
           Reference_Allele, Tumor_Seq_Allele2, Protein_Change, VAF, `CGC_Cancer Somatic Mut`) %>% 
    unite(mutation_id, Hugo_Symbol, Protein_Change, sep=":")
  
	# Spread into sample VAF matrix
  mut_mat <- subset %>% spread(key = Tumor_Sample_Barcode, value = VAF,
                               fill = 0)
  
	# Write output files
  mut_mat %>% write_tsv(paste0(output_prefix, "_all_coding_", paste(samples, collapse = "_"), ".tsv"))
  mut_mat %>% filter_at(vars(samples), any_vars(. == 0)) %>% 
    write_tsv(paste0(output_prefix, "_unique_coding_", paste(samples, collapse = "_"), ".tsv"))
  mut_mat %>% filter_at(vars(samples), all_vars(. != 0)) %>% 
    write_tsv(paste0(output_prefix, "_shared_coding_", paste(samples, collapse = "_"), ".tsv"))
  
	# Binarize into absent/present matrix
  bin_mat <- mut_mat
  bin_mat[, c("A139-PDX", "TS10", "TS32")][bin_mat[, c("A139-PDX", "TS10", "TS32")] > 0] <- 1
  
  bin_mat <- bin_mat %>% mutate(rsum =rowSums(.[c("A139-PDX", "TS10", "TS32")]))
  
  total_overlap <- bin_mat %>% filter(rsum == 3)
  unique_mut <- bin_mat %>% filter(rsum == 1)
	print(paste0("Shared mutations: ", total_overlap))
	print(paste0("Unique mutations: ", unique_mut))
}

# Given an early late TSV file, find out the gene interested in early or late situation
get_earlylate_gene <- function(earlylate, gene_list){
  # Use this function to get a tibble that spreads into early and late for the gene list input, e.g.
  # get_earlylate_gene(earlylate_dataframe, c("EGFR", "TP53"))
  earlylate_gene_df <- tibble(Tumor_Sample_Barcode=unique(earlylate$Tumor_Sample_Barcode))
  for (gene_interested in gene_list){
    var_name <- paste0(gene_interested, "_Protein_Change")
    var_name_early <- paste0("Early_", gene_interested)
    var_name_late <- paste0("Late_", gene_interested)
    
    earlylate_gene_df <- earlylate %>% 
      filter(Variant_Classification %in% coding_def, Hugo_Symbol==gene_interested) %>% 
      group_by(Tumor_Sample_Barcode, comb.timing) %>% 
      summarise(!! var_name := paste(Protein_Change, collapse = ";")) %>% 
      spread(key = comb.timing, value = !!var_name) %>% 
      rename(!!var_name_early := "Early", !!var_name_late := "Late") %>% 
      full_join(earlylate_gene_df, by="Tumor_Sample_Barcode")
  }
  return(earlylate_gene_df)
}

load_drivers <- function(){
  # Use known drivers list
  luaddrivers <- read_tsv("~/gis-scripts/references/96_drivers_LUAD.txt", col_types = cols()) 
  tcgadrivers <- read_tsv("~/gis-scripts/references/tcga_2018_driver_list.txt", col_types = cols()) %>% 
    filter(Cancer %in% c("LUAD", "PANCAN"))
  cgcdrivers <- read_csv("~/gis-scripts/references/cancer_gene_census.csv", col_types=cols()) 
  # Exclude fusion, promoter and translocation (only) drivers. Select somatic only
  cgcdrivers <- cgcdrivers %>% filter(str_detect(`Mutation Types`, "Mis|N|F|D|S|O|A") & Somatic == "yes")
  # signalling pathway
  signalling <- read_tsv("~/gis-scripts/references/TCGA_2018_signalling_pathways.txt", col_types=cols())
  # Cancermine. Gives priority to NSCLC
  priority_cancermine <- c("non-small cell lung carcinoma","lung cancer",
			   "lung adenocarcinoma","lung small cell carcinoma","lung squamous cell carcinoma",
			   "lung carcinoma","lung giant cell carcinoma","lung large cell carcinoma",
			   "lung papillary adenocarcinoma","lung oat cell carcinoma","mucinous lung adenocarcinoma")
  cancermine <- read_tsv("~/gis-scripts/references/cancermine_collated_2018-10-17.tsv", col_types=cols()) %>%
    filter(grepl("lung", cancer_normalized)) 
  cancermine$cancer_normalized <- factor(cancermine$cancer_normalized, levels=priority_cancermine)
  cancermine <- cancermine %>% arrange(cancer_normalized) %>%
    group_by(gene_normalized) %>% mutate(rnum = row_number()) %>% 
    filter(rnum==1) %>% select(-rnum)

  # Drivers to subset
  all_drivers <- c(luaddrivers$`Tier1 LUAD updated`, 
		   tcgadrivers$Gene, cgcdrivers$`Gene Symbol`,
		   signalling$Hugo_Symbol,
		   # For cancermine, only load those with more than 3 citations
		   cancermine$gene_normalized[cancermine$citation_count>=3]) %>% unique

  return(list(all_drivers = all_drivers,
	      luaddrivers = luaddrivers,
	      tcgadrivers = tcgadrivers,
	      cgcdrivers = cgcdrivers,
	      signalling = signalling,
	      cancermine = cancermine))
}

# Oncoprint sort function taken from github Armish:
# https://gist.github.com/armish/564a65ab874a770e2c26
# If binarize = TRUE, automatically change everything that's not empty string to 0
memoSort <- function(M, binarize=TRUE) {
  # browser()
  # Handle those with only 1 driver
  # Binarize matrix for sorting
  if(binarize){
    M[M!=""] <- 1
    M[M==""] <- 0
    M <- apply(M, 2, as.numeric)
  }
  if (nrow(M) == 1){
    return(M)
  }
  # run as.data.frame to convert spaces to dot in rownames
  M <- as.data.frame(M)
  o_M <- M
  M <- M
  # Sort with binarized matrix
  M[M!=0] <- 1
  geneOrdering <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE);
  geneOrder <- geneOrdering$ix
  geneOrderNames <- names(geneOrdering$x)
  scoreCol <- function(x) {
    score <- 0;
    for(i in 1:length(x)) {
      if(x[i]) {
        score <- score + 2^(length(x)-i);
      }
    }
    return(score);
  }
  scores <- apply(M[geneOrder, ], 2, scoreCol);
  sampleOrdering <- sort(scores, decreasing=TRUE, index.return=TRUE);
  sampleOrder <- sampleOrdering$ix
  sampleNameOrder <- names(sampleOrdering$x)
  # Group the private mutations together
  o_M <- o_M[geneOrder, sampleOrder]
  privateMut <- names(rowSums(M)[rowSums(M) == 1])
  privateMut_mat <- M[privateMut,,drop=FALSE]
  # Set aside non private mutations
  otherMut <- geneOrderNames[!geneOrderNames %in% privateMut]
  # Sort by sample order if needed
  if(is.data.frame(privateMut_mat) | is.matrix(privateMut_mat)){
    privateMut_mat <- as.data.frame(privateMut_mat)
    privateMut_mat <- arrange(rownames_to_column(privateMut_mat, var = "genes"),
                              !!!rlang::syms(sampleNameOrder))
    privateOrder <- privateMut_mat$genes
  } else{
    # Sometimes there's no private mutations
    privateOrder <- NULL
  }

  # Return correct order
  o_M <- o_M[c(otherMut, privateOrder), ]

  return(o_M);
  # return(list(geneOrder, sampleOrder));
}

convert_lm_to_formula <- function(lm_model){
  ## Use this function to convert a linear formula to a string for annotation
  ## purpose. Takes in a summary(lm) object.
  ## Should work for multiple regression as well as it loops through all
  ## coefficients.
  lm_summary <- summary(lm_model)
  lm_confint <- confint(lm_model)
  formula <- paste0(rownames(attr(terms(lm_summary), "factors"))[1], " = ")
  for (i in 2:nrow(lm_summary$coefficients)){
    formula <- paste0(formula, formatC(lm_summary$coefficients[i, 'Estimate'], format="e", digit=1), 
                      rownames(lm_summary$coefficients[i,,drop=FALSE]),
                      "(p = ", formatC(lm_summary$coefficients[i, 'Pr(>|t|)'], format="e", digit=1), 
                      "; 95CI = ", formatC(lm_confint[i, '2.5 %'], format="e", digit=1), ' - ',
                      formatC(lm_confint[i, '97.5 %'], format="e", digit=1),
                      ")")
  }
  formula <- paste0(formula, " + ", formatC(lm_summary$coefficients[1, 'Estimate'], format="e", digit=1))
  formula <- paste0(formula, "\nResidual =  ", formatC(lm_summary$sigma, format="e", digit=1))
  formula <- paste0(formula, ", Adj R^2 =  ", formatC(lm_summary$adj.r.squared, format="e", digit=1))
}

#' This function output a ggplot of timeline visualization.
#' The input is a dataframe of:
#' content, start, end, and group
#' where the y-axis will be split into groups (e.g. category of event),
#' and content is the event name to annotate. 
timevis_gg <- function(timelineDF, title=""){
  require(lubridate)
  # For clinical record coloring consistent with copy number plots
  c5 <- c("#FB9A99", "maroon", "palegreen2",
	  "red", "steelblue4", "darkturquoise")
  # Plot treatments duration and other date separate in dots and segments respectively
  dots <- timelineDF %>% filter(start==end)
  non_dots <- timelineDF %>% filter(start!=end) %>% 
    mutate(interval=interval(start, end), midpoint = as.Date(start + (as.duration(interval)/2)))
  
  clin_vis <- ggplot() + 
    geom_linerange(data=non_dots, aes(x=group, ymin=start, ymax=end, group=content, color=group), 
                   size=3, alpha=1, position = position_dodge(width=0.5), show.legend = FALSE) +
    geom_quasirandom(data=dots, aes(y=start, x=group, color=group), 
                     show.legend = FALSE, groupOnX = TRUE, size=3) +
    geom_text_repel(data=dots, aes(y=start, x=group, label=content), 
                    min.segment.length = 0, alpha=1,
                    position=position_quasirandom(groupOnX = TRUE)) +
    geom_text_repel(data=non_dots, aes(x=group, y=midpoint, label=content, group=content),
                    position = position_dodge(width=0.5), alpha=1) +
    geom_vline(aes(xintercept = seq(1.5, length(unique(tmp$group)) -0.5, 1)), linetype=2) +
    scale_y_date(date_breaks = "6 months", date_minor_breaks = "3 months") +
    scale_color_manual(values=c5) +
    coord_flip() +
    labs(x="Types", y="Date", title = title) +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, vjust=0.5))
  
  return(clin_vis)
}

#' Small function to get max, median and max of a vector
minMedmax <- function(numVec, removeNA=TRUE){
  res <- c(
    min=min(numVec, na.rm=removeNA),
    max=max(numVec, na.rm=removeNA),
    median=median(numVec, na.rm=removeNA)
  )
  return(res)
}

#' choose color from c25 (excluding black as I don't like black...)
#' priority from 8 to 15 as I like the pastel color
choose_color <- function(vec, c25, c25orNot=TRUE){
  c25 <- c25[c25!="black"]
  colList <- c()
  # Remove NA (ComplexHeatmap will color with grey)
  vec <- vec[!is.na(vec)]
  # Get unique factors
  fac <- unique(factor(vec))
  # Convert to numeric
  numvec <- as.numeric(fac)
  fac <- as.character(fac)
  # If using c25 colors
  if(c25orNot){
    if(length(fac) < 8){
      colorPal <- c25[7:14]
      for (i in 1:length(fac)){
        colList[fac[i]] = colorPal[numvec[i]]
      }
    } else if (length(fac) >= 8 & length(fac) < 24){
      colorPal <- c25
      for (i in 1:length(fac)){
        colList[fac[i]] = colorPal[numvec[i]]
      }
    } else{
      print("Too many unique categories to map!")
      return(0)
    }
  } else{
    # Directly map if not using c25
    for (i in 1:length(fac)){
      colorPal <- c25
      colList[fac[i]] = colorPal[numvec[i]]
    }
  }
  return(colList)
}

#' A small function to remove duplicate genes on an
#' expression matrix by sorting to MAD and keeping the higher
#' MAD gene. exprMat is a gene matrix with gene names as the 
#' row names and samples as column names.
dedup_RNA_byMAD_toMatrix <- function(expr, samples, geneCol="rsem_geneNames"){
  exprMat <- as.matrix(expr[, samples])
  rownames(exprMat) <- expr[, geneCol, drop=TRUE]
  mad <- apply(exprMat, 1, mad)
  exprMat <- exprMat[order(mad, decreasing = TRUE),,drop=FALSE]
  exprMat <- exprMat[!duplicated(rownames(exprMat)),,drop=FALSE]
  return(exprMat)
}

# A set of colors for complex heatmap
# Generate color palette for bottom annotation
anno_pal <- RColorBrewer::brewer.pal(4, name="Dark2")
anno_pal10 <- viridis::viridis(10)
anno_pal2 <- c('#f2f0f7','#6a51a3')
anno_pal4 <- c('#f2f0f7','#cbc9e2','#9e9ac8','#6a51a3')
better_redBlue <- c('#ef8a62','#FFFFFF','#67a9cf')

#' N-1 chi-square test recommended by Frank Harrell:
#' https://stats.stackexchange.com/questions/14226/given-the-power-of-computers-these-days-is-there-ever-a-reason-to-do-a-chi-squa
chisq_nMinusOne <- function(table){
  n = sum(table)
  stats <- chisq.test(table, correct=FALSE)
  # Don't do yates correction (Over-correction!)
  pval <- 1 - pchisq(stats$statistic * ((n-1)/n), df = stats$parameter)
  return(pval)
}

#' Function to do upper-quartile normalization on RNA-seq expression matrix
upper_quartile_norm <- function(expr_mat, method="all"){
  # Remove row with all zero
  expr_mat <- expr_mat[rowSums(expr_mat > 0) > 0, ]
  if(method=="all"){
    quantile <- apply(expr_mat, 2,
                      function(x){quantile(x, 0.75)})
  } else{
    quantile <- apply(expr_mat, 2,
                      function(x){quantile(x[x > 0], 0.75)})
  }
  expr_mat <- t(t(expr_mat) / quantile)
  # Multiple by mean of upper quantile for higher values, nicer for analysis
  expr_mat <- expr_mat * mean(quantile)
  return(expr_mat)
}

#' Function to regress gene. Purities should be a dataframe with samples as rownames,
#' a column named purity, and any extra column specified in color, shape and wrap
regress_gene <- function(expr_mat, purity_df, gene="ROS1",
                         color=NULL, shape=NULL, wrap=NULL){
  gene_expr <- t(expr_mat[gene,, drop=FALSE])
  # Remove samples with no purity in purity_df
  gene_expr <- gene_expr[intersect(rownames(purity_df), rownames(gene_expr)),,drop=FALSE]
  
  gene_expr <- as_tibble(gene_expr, rownames = "subject_id")
  purity_df <- as_tibble(purity_df, rownames = "subject_id")
  
  # Merge purities
  gene_expr <- gene_expr %>% 
    inner_join(purity_df, by="subject_id")

  # Handle special character such as NKX2-1 for aes_string with backtick
  gene = paste0("`", gene, "`")
  
  ggplot(gene_expr, aes_string(x="purity", y=gene)) +
    geom_point(aes_string(color=color, shape=shape), show.legend = TRUE) +
    facet_wrap(wrap) +
    # scale_shape_manual(values=c(4,2)) +
    # lims(x=c(0, 1)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
    geom_smooth(method="lm", fullrange=TRUE) +
    labs(title = paste0(gene, " Expression against Purity"), x="Tumor Purity") +
    theme(axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5),
          panel.spacing = unit(2, "lines"))
}

#' Function to plot PCA with ggplot and color according to dataframe
#' gene_mat is a matrix with genes as rownames and samples as colnames,
#' annotateDF is a dataframe with samples as rownames and any columns
#' to be used for color and shape
plotPCA_colorShape <- function(gene_mat, annotateDF, color=NULL, shape=NULL, size=NULL){
  # Get only samples with information
  # gene_mat <- gene_mat[, intersect(rownames(gene_mat), rownames(annotateDF)), drop=FALSE]
  pca_choice <- prcomp(t(gene_mat), scale=TRUE)
  # Get PCA variance for ggplot annotation
  var_pc <- pca_choice$sdev * pca_choice$sdev
  var_pc <- var_pc / sum(var_pc)
  
  sample_pc <- tibble(subject_id = rownames(pca_choice$x),
                            PC1 = pca_choice$x[, "PC1"],
                            PC2 = pca_choice$x[, "PC2"])
  annotateDF <- as_tibble(annotateDF, rownames="subject_id")
  
  toPlot <- annotateDF %>% 
    inner_join(sample_pc, by="subject_id")
  
  ggplot(toPlot, aes(x=PC1, y=PC2)) +
    geom_point(aes_string(color=color, shape=shape, size=size)) +
    labs(x = paste0("PC1:" , round(var_pc[1] * 100, 2), "%"), 
         y = paste0("PC2: ", round(var_pc[2] * 100, 2), "%")) +
    coord_fixed()
}

# Function to plot submap p-values matrix using viridis colormap
# row_labels and col_labels should be the name of the subtypes used
# for dataset A and dataset B respectively when running submap
plot_submap <- function(submap_res, row_labels, col_labels){
  submap_res <- read_tsv(submap_res,
                         skip = 2)
  submap_res_mat <- as.matrix(submap_res[3:ncol(submap_res)])
  colnames(submap_res_mat) <- col_labels
  rownames(submap_res_mat) <- row_labels
  # Annotate p-value with customized label
  dim_mat <- dim(submap_res_mat)
  label_pval <- matrix(paste(matrix("p = ", nrow=dim_mat[1], ncol=dim_mat[2]), 
                             signif(submap_res_mat, 2), sep=""), 
                       nrow=dim_mat[1], ncol=dim_mat[2])
  pheatmap(submap_res_mat,
           cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = label_pval, color = viridis(20, direction = -1),
           cellwidth = 50, cellheight = 50, legend = FALSE)
}

# Given a vector, generate all combination of 2 for stat_compare_means
generate_combination <- function(vec){
  combination <- combn(unique(as.character(vec)), 2)
  list_com = list()
  for (i in 1:ncol(combination)){
    list_com[[i]] <- c(combination[1, i], combination[2, i])
  }
  return(list_com)
}

## This code is used to obtain the annotation for ENID using HGNC database and
## mygene query service (Different results!). A generated annotation data frame
## can be found under references/annotation_RSEM_ENID.tsv
# annotation <- data.frame(row.names = merged_count$gene_id, gene_name=merged_count$gene_name,
#                          gene_type=merged_count$gene_type)
# hgnc_annot <- read.delim("~/mount/aquila_home/references/hgnc_complete_set_2018-10-1.txt")
# annotation$ENID <- gsub("(ENSG[0-9]+)\\..*", "\\1", rownames(annotation))
# 
# annotation$HGNC_entrezID <- hgnc_annot[match(annotation$ENID, hgnc_annot$ensembl_gene_id), 'entrez_id']
# annotation$HGNC_symbol <- hgnc_annot[match(annotation$ENID, hgnc_annot$ensembl_gene_id), 'symbol']
# # Sanity check on HGNC annotation vs RSEM pipeline annotation, HGNC is updated to latest version
# # so there could be quite a lot of gene names change!
# annotation[which(as.character(annotation$gene_name) != as.character(annotation$HGNC_symbol)),]
# 
# # Check for genes without HGNC symbol or duplicated entrezID, here all of them are actually NA
# # 21842 cannot be mapped to HGNC database
# annotation[which(duplicated(annotation$entrezID) & !is.na(annotation$entrezID)),]
# 
# # Test with mygene
# mygene_annot <- queryMany(annotation$ENID, scopes="ensembl.gene",
#                           fields=c("entrezgene", "symbol"), species="human", returnall=TRUE)
# 
# annotation$mygene_entrez <- mygene_annot$response[match(annotation$ENID, mygene_annot$response$query), "entrezgene"]
# annotation$mygene_symbol <- mygene_annot$response[match(annotation$ENID, mygene_annot$response$query), "symbol"]
# 
# annotation[which(duplicated(annotation$mygene_entrez)),]
# annotation[which(as.character(annotation$mygene) != as.character(annotation$HGNC_symbol)),]
# 
# write.table(annotation, "annotation_mygene_hgnc.tsv", sep="\t",
#             quote=FALSE)
# Use biomart
# library(biomaRt)
# grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
#                  path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
# 
# genes.with.id=getBM(attributes=c("entrezgene", "external_gene_name"),
#                     values=genes$ENID, mart= grch37, filters="ensembl_gene_id") # fuction to get  gene id's and gene name from data base

