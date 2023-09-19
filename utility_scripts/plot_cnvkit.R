# This script is used to plot the CNVkit output. It is based on the
# custom function scripts written by Khi Pin

# Run this script from the command line with the following command:
# Rscript plot_cnvkit.R <BEDNAME>
# Bed must have the following columns: chromosome, start, end, cn, log2
# Requires tidyverse and data.table R packages to be installed

# Load libraries
library(tidyverse)
library(data.table)
source("/Users/khipinchua/pb_bitbucket/wdl-hifisomatic/utility_scripts/R_functions/copy_number_functions.R")
theme_set(theme_bw(base_size = 14))

cytoband <- load_cytoband("/Users/khipinchua/pb_bitbucket/wdl-hifisomatic/utility_scripts/R_functions/cytoBand_hg38.txt.gz")

# Load in bed from command line argument
args <- commandArgs(trailingOnly = TRUE)
BEDGRAPHNAME <- args[1]

# Read all bedgraph files into a dataframe. For each
# bedgraph, set colnames to be chromosome, Start, End, CN
bedgraph <- fread(
    BEDGRAPHNAME,
    select = c("chromosome", "start", "end", "cn", "log2"),
    # Start and end are integer
    colClasses = list(
        character = c("chromosome"),
        numeric = c("start", "end"),
        numeric = c("cn", "log2")
    )
) %>%
    mutate(chromosome = gsub("chr(.*)", "\\1", chromosome))

bedgraph <- set_absolute_pos(bedgraph, cytoband_length = cytoband$cytoband_length)

# If cn is more than 5, set to 5 and add a label column for plotting
bedgraph <- bedgraph %>%
    mutate(
        cn = ifelse(cn > 5, 5, cn),
        label = ifelse(cn > 5, cn, "")
    )

plot_cn <- ggplot(bedgraph, aes(x = Start, y = cn)) +
    geom_segment(aes(x = Start, xend = End, y = cn, yend = cn),
        show.legend = FALSE, linewidth = 1
    ) +
    # Annotate with the label column
    annotate(
        geom = "text",
        x = bedgraph$Start,
        y = bedgraph$cn,
        label = bedgraph$label,
        size = 3
    ) +
    scale_x_continuous(
        breaks = cytoband$cytoband_arm$midpoint,
        labels = cytoband$cytoband_arm$chromosome,
        minor_breaks = NULL,
        limits = c(min(cytoband$cytoband_arm$Start_p), max(cytoband$cytoband_arm$End_q)),
        # No left and right padding with expand
        expand = c(0, 0)
    ) +
    geom_vline(
        xintercept = unique(
            cytoband$cytoband_arm$Start_p,
            cytoband$cytoband_arm$End_q
        ),
        linetype = 2, alpha = 0.5
    ) +
    labs(x = "", y = "Copy Number") +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, angle = -40))

plot_depth <- ggplot(bedgraph, aes(x = Start, y = log2)) +
    geom_segment(aes(x = Start, xend = End, y = log2, yend = log2),
        show.legend = FALSE, linewidth = 1
    ) +
    scale_x_continuous(
        breaks = cytoband$cytoband_arm$midpoint,
        labels = cytoband$cytoband_arm$chromosome,
        minor_breaks = NULL,
        limits = c(min(cytoband$cytoband_arm$Start_p), max(cytoband$cytoband_arm$End_q)),
        # No left and right padding with expand
        expand = c(0, 0)
    ) +
    geom_vline(
        xintercept = unique(
            cytoband$cytoband_arm$Start_p,
            cytoband$cytoband_arm$End_q
        ),
        linetype = 2, alpha = 0.5
    ) +
    scale_y_continuous(
        limits = c(-1, 1)) +
    labs(x = "", y = "Log2 Ratio") +
    theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, angle = -40))

# Remove extension from BEDGRAPHNAME.
# Extension can be anything
ofilename <- str_extract(
        BEDGRAPHNAME, "(.*)\\..*", group = 1
    )

ggsave(
    filename = paste0(ofilename, ".cn.pdf"),
    plot = plot_cn,
    width = 16, height = 8, useDingbats = FALSE
)

ggsave(
    filename = paste0(ofilename, ".log2_depth.pdf"),
    plot = plot_depth,
    width = 16, height = 8, useDingbats = FALSE
)

