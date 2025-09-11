# %%
import pycirclize
from pycirclize import Circos
from pycirclize.utils import load_eukaryote_example_dataset
import pandas as pd
import re
import argparse
import json

def expand_info_field(df):
    # Split INFO field into individual components
    info_fields = df['INFO'].str.split(';')

    # Dictionary to store all key-value pairs
    info_dict = {}

    # Process each row
    for idx, fields in enumerate(info_fields):
        row_dict = {}
        if fields is not None:  # Handle potential NaN values
            for field in fields:
                if '=' in field:  # Only process fields with key=value format
                    key, value = field.split('=', 1)  # Split on first '=' only
                    row_dict[key] = value
        info_dict[idx] = row_dict

    # Create DataFrame from dictionary
    info_df = pd.DataFrame.from_dict(info_dict, orient='index')

    # Merge expanded columns with original DataFrame
    result_df = pd.concat([df, info_df], axis=1)

    return result_df

def extract_chr2_pos2(df):
    # Filter for BND entries
    bnd_df = df[df['SVTYPE'] == 'BND'].copy()  # Use copy() to avoid warnings
    # If empty, return empty DataFrame
    if bnd_df.empty:
        return bnd_df

    # Extract chr2 and pos2 using regex
    # ALT field format example: ]chr2:123456]N
    pattern = r'[\[\]]?(chr\w+):(\d+)[\[\]]'

    # Create new columns using .loc to avoid warnings
    bnd_df.loc[:, 'CHR2'] = bnd_df['ALT'].str.extract(pattern)[0]
    bnd_df.loc[:, 'POS2'] = bnd_df['ALT'].str.extract(pattern)[1].astype(int)

    # If BCSQ column exists, extract gene name (second field split by |)
    if 'BCSQ' in bnd_df.columns:
        bnd_df.loc[:, 'GENE'] = bnd_df['BCSQ'].str.split('|').str[1]
        # If NA, replace with empty string
        bnd_df.loc[:, 'GENE'] = bnd_df['GENE'].fillna('')

    return bnd_df

def create_fusion_annotation(df):
    # Create empty FUSION column using .loc to avoid warnings
    df = df.copy()  # Make a copy to avoid modifying original
    df.loc[:, 'FUSION'] = ''

    # Filter for BND records only
    bnd_records = df[df['SVTYPE'] == 'BND'].copy()

    for idx, row in bnd_records.iterrows():
        # Find mate record using MATEID
        mate_records = df[df['ID'] == row['MATE_ID']]
        mate = mate_records.iloc[0] if not mate_records.empty else None

        # Check if mate exists and both have GENE information
        if mate is not None and row['GENE'] != '' and mate['GENE'] != '':
            # If both genes are the same, skip
            if row['GENE'] == mate['GENE']:
                continue
            # Sort genes alphabetically and join with ::
            genes = sorted([row['GENE'], mate['GENE']])
            fusion = '::'.join(genes)

            # Update FUSION column for both records using .loc
            df.loc[idx, 'FUSION'] = fusion
            mate_idx = mate_records.index[0]
            df.loc[mate_idx, 'FUSION'] = fusion

    # Fill empty FUSION with empty string
    df.loc[:, 'FUSION'] = df['FUSION'].fillna('')

    return df

# %%

def save_empty_circos(output_file):
    circos_empty = Circos.initialize_from_bed("/app/hg38.bed", space = 3)
    circos_empty.text("No SVs found", size=15)
    fig = circos_empty.plotfig(dpi = 300)
    fig.savefig(output_file, dpi = 300, bbox_inches = 'tight')

def canonicalize_breakpoints(row) -> str:
    """
    Return a canonical CHR:POS-CHR2:POS2 string for a mate pair so that
    A-B and B-A are represented identically. Prefer parsing from MATE_PAIR
    if available; otherwise, fall back to CHROM/POS and CHR2/POS2.
    """
    endpoints = []
    try:
        # Try parsing from MATE_PAIR like "chr1:123-chr2:456"
        if 'MATE_PAIR' in row and isinstance(row['MATE_PAIR'], str) and row['MATE_PAIR']:
            matches = re.findall(r'(chr[\w]+):(\d+)', row['MATE_PAIR'])
            if len(matches) >= 2:
                endpoints = [(matches[0][0], int(matches[0][1])), (matches[1][0], int(matches[1][1]))]
    except Exception:
        # Fall back to CHROM/POS pairs on any parsing issue
        endpoints = []

    if not endpoints:
        try:
            endpoints = [
                (str(row['CHROM']), int(row['POS'])),
                (str(row['CHR2']), int(row['POS2']))
            ]
        except Exception:
            # As a last resort, build with raw strings without sorting
            return f"{row.get('CHROM', '')}:{row.get('POS', '')}-{row.get('CHR2', '')}:{row.get('POS2', '')}"

    # Sort by chromosome string then position to canonicalize order
    endpoints.sort(key=lambda x: (x[0], x[1]))
    (c1, p1), (c2, p2) = endpoints
    return f"{c1}:{p1}-{c2}:{p2}"

def main():
    parser = argparse.ArgumentParser(description='Generate circos plot from VCF file')
    parser.add_argument('vcf_file', help='Input VCF file (can be .vcf or .vcf.gz)')
    parser.add_argument(
        'sample_name', help='Sample name to be used in the plot',
        default='SAMPLE'
        )
    parser.add_argument(
        'mitelman_mcgene', help='Mitelman fusion database MCGENE file'
    )
    parser.add_argument(
        '--max-gene-labels', type=int, default=60,
        help='Maximum number of gene labels to display (default: 60). Set to 0 for no limit.'
    )
    parser.add_argument(
        '--chr-label-size', type=int, default=10,
        help='Chromosome label font size (default: 10)'
    )
    parser.add_argument(
        '--sector-space', type=int, default=5,
        help='Space between chromosomes/sectors (default: 5)'
    )
    parser.add_argument(
        '--link-alpha', type=float, default=0.8,
        help='Transparency for translocation lines (0.0-1.0, default: 0.8)'
    )
    args = parser.parse_args()

    # Generate output filename by replacing .vcf or .vcf.gz with .pdf
    output_file = args.sample_name + '.pdf'
    output_png = args.sample_name + '.png'
    # Output table for Mitelman fusion pairs
    mitelman_out = args.sample_name + '.mitelman_fusions.tsv'
    # Output table for all sample pairs that involve any Mitelman-known gene
    # Note: do not use 'mitelman' in filename since pairs may not be Mitelman fusions
    mitelman_gene_pairs_out = args.sample_name + '.known_gene_pairs.tsv'

    # Initialize output TSVs as header-only to ensure files exist
    # even if we return early (no SVs or no BNDs found).
    pd.DataFrame(columns=["Fusion", "Breakpoints"]).to_csv(mitelman_out, sep='\t', index=False)
    pd.DataFrame(columns=["Gene Pairs", "Breakpoints"]).to_csv(mitelman_gene_pairs_out, sep='\t', index=False)

    # Read VCF file using pandas. Use "##" to skip the header. columns are
    # #CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  NPM30000DD  NPM30000DE
    try:
        df = pd.read_csv(args.vcf_file, sep="\t",
                        header = 0, comment = "#")
    except pd.errors.EmptyDataError:
        # Initialise empty df if no SV is found in the dataframe
        df = pd.DataFrame()

    # If df is empty, save a PDF that says "No SVs found"
    if df.empty:
        save_empty_circos(output_file)
        return

    # Set column names
    df.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
    df = expand_info_field(df)
    bnd_only = extract_chr2_pos2(df)
    # If bnd_only is empty, save empty fig
    if bnd_only.empty:
        save_empty_circos(output_file)
        return

    bnd_only = create_fusion_annotation(bnd_only)
    # # Remove all intrachromosomal stuff, keep only interchromosomal
    # bnd_only = bnd_only[bnd_only['CHR2'] != bnd_only['CHROM']]
    # For each BND, find the mate with MATEID column (match ID with MATEID),
    # and check if GENE of the mate is not empty. If not empty, add
    # a column called "FUSION" by concatenating GENE and GENE of the mate (sort and sep with ::)


    # Require POS2 to be at least 100 kbp away from POS
    bnd_only = bnd_only[abs(bnd_only['POS2'] - bnd_only['POS']) > 100000]
    # If empty, save empty fig
    if bnd_only.empty:
        save_empty_circos(output_file)
        return

    # Only keep a single event for Saumya
    # bnd_only = bnd_only[bnd_only['POS'] == 86936276]

    # Remove BND from original df
    df = df[df['SVTYPE'] != 'BND']
# %%

# %%
    # Load hg38 dataset (https://github.com/moshi4/pycirclize-data/tree/main/eukaryote/hg38)
    circos = Circos.initialize_from_bed("/app/hg38.bed", space=args.sector_space)

    circos.text("All translocations (GRCh38)", size=15)

    # Add cytoband tracks from cytoband file
    circos.add_cytoband_tracks((95, 100), "/app/hg38_cytoband.tsv")

    # Read cytoband TSV
    cytoband = pd.read_csv("/app/hg38_cytoband.tsv",
                        sep = "\t")

    # Read mitelman gene list
    # # Firstly, read the json schema file
    # schema = pd.read_json("/Users/khipinchua/softwares/references/mitelman_db/MCGENE.JSON.SCHEMA")

    # column names from schema "name" column. Read the mitelman DB gene list using the JSON schema
    mitelman_db = pd.read_csv(args.mitelman_mcgene, sep = "\t")
    # Find only rows where Gene contains "::"
    mitelman_db = mitelman_db[mitelman_db['Gene'].str.contains("::")]
    # Split the gene column by "::" and sort, then merge back
    mitelman_db['Fusion'] = mitelman_db['Gene'].str.split("::").apply(lambda x: "::".join(sorted(x)))
    # Group by fusion, and keep only if N>=3
    mitelman_db = mitelman_db.groupby('Fusion').filter(lambda x: len(x) >= 3)
    unique_mitelman_fusion = mitelman_db['Fusion'].unique()
    # Build set of known genes from Mitelman fusions
    mitelman_known_genes = set()
    for fusion in unique_mitelman_fusion:
        parts = fusion.split('::')
        for g in parts:
            if g:
                mitelman_known_genes.add(g)

    # Subset df to include only chromosomes that's in cytoband #chrom column
    bnd_only = bnd_only[bnd_only['CHROM'].isin(cytoband['#chrom'])]
    bnd_only = bnd_only[bnd_only['CHR2'].isin(cytoband['#chrom'])]

    # If empty, save empty fig
    if bnd_only.empty:
        save_empty_circos(output_file)
        return

    # Build canonical breakpoint strings per fusion, ensuring A-B and B-A collapse
    bnd_with_fusion = bnd_only[bnd_only['FUSION'].fillna('') != '']
    fusion_to_breakpoints = {}
    if not bnd_with_fusion.empty:
        for _, r in bnd_with_fusion.iterrows():
            # Always canonicalize to avoid duplicates from reciprocal records
            bp = canonicalize_breakpoints(r)
            fusion_to_breakpoints.setdefault(r['FUSION'], set()).add(bp)

    # Now that bnd_only has been subset to valid chromosomes, write matched fusions table
    sample_fusions = bnd_only['FUSION'].fillna('')
    sample_fusions = sample_fusions[sample_fusions != '']
    if not sample_fusions.empty:
        matched = sorted(set(sample_fusions) & set(unique_mitelman_fusion))
        matched_rows = []
        for f in matched:
            bps = sorted(fusion_to_breakpoints.get(f, []))
            matched_rows.append({"Fusion": f, "Breakpoints": ";".join(bps)})
        pd.DataFrame(matched_rows).to_csv(mitelman_out, sep='\t', index=False)

        # Also save all sample pairs where at least one gene is known in Mitelman
        def involves_known_gene(fusion: str) -> bool:
            genes = fusion.split('::') if isinstance(fusion, str) and fusion else []
            return any(g in mitelman_known_genes for g in genes)

        involved_pairs = sorted([f for f in set(sample_fusions) if involves_known_gene(f)])
        if involved_pairs:
            involved_rows = []
            for f in involved_pairs:
                bps = sorted(fusion_to_breakpoints.get(f, []))
                involved_rows.append({"Gene Pairs": f, "Breakpoints": ";".join(bps)})
            pd.DataFrame(involved_rows).to_csv(
                mitelman_gene_pairs_out, sep='\t', index=False
            )

    # Count total number of genes to label across all chromosomes
    total_genes = 0
    for sector in circos.sectors:
        bnd_chr = bnd_only[bnd_only['CHROM'] == sector.name]
        bnd_chr = bnd_chr[bnd_chr['FUSION'] != '']
        total_genes += len(bnd_chr)

    # Set maximum number of gene labels to display (0 means no limit)
    max_gene_labels = args.max_gene_labels
    should_label_genes = (max_gene_labels == 0) or (total_genes <= max_gene_labels)

    # Function to get clean chromosome name
    def get_chr_display_name(sector_name):
        # Remove "chr" prefix and return clean chromosome name
        return sector_name.replace('chr', '')

    # Plot chromosome labels on segments
    for i, sector in enumerate(circos.sectors):
        display_name = get_chr_display_name(sector.name)

        # Place label at the leftmost edge and lower than the chromosome segment
        # Position at the very start of the segment and below the cytoband track
        sector.text(
            display_name,
            x=sector.start + (sector.size * 0.02),  # Just 2% into the segment (leftmost edge)
            r=90,  # Below the cytoband track (95-100) for better separation
            size=args.chr_label_size,
            weight='bold',
            ha='left',  # Left-align text to keep it at the edge
            va='center'  # Vertically center
        )
        # Annotate gene name only if total count is reasonable
        if should_label_genes:
            bnd_chr = bnd_only[bnd_only['CHROM'] == sector.name]
            # Only label if fusion is not empty
            bnd_chr = bnd_chr[bnd_chr['FUSION'] != '']

            # Priority system: if we're approaching the limit, prioritize known fusions
            if not bnd_chr.empty and total_genes > max_gene_labels * 0.8 and max_gene_labels > 0:
                # Prioritize known fusions when approaching limit
                known_fusions = bnd_chr[bnd_chr['FUSION'].isin(unique_mitelman_fusion)]
                unknown_fusions = bnd_chr[~bnd_chr['FUSION'].isin(unique_mitelman_fusion)]

                # Take all known fusions plus some unknowns to fill quota
                remaining_quota = max(1, max_gene_labels // len(circos.sectors))
                if len(unknown_fusions) > remaining_quota:
                    unknown_fusions = unknown_fusions.head(remaining_quota)

                bnd_chr = pd.concat([known_fusions, unknown_fusions]).sort_values('POS')

            if not bnd_chr.empty:
                bnd_track = sector.get_track('cytoband')

                # Collect all gene positions and labels for this chromosome
                gene_positions = bnd_chr['POS'].tolist()
                gene_labels = bnd_chr['GENE'].tolist()

                # Sort genes by position to help with overlap resolution
                sorted_data = sorted(zip(gene_positions, gene_labels))
                gene_positions, gene_labels = zip(*sorted_data) if sorted_data else ([], [])

                # Group nearby genes to prevent overlapping labels (increased distance)
                grouped_annotations = []
                current_group_pos = []
                current_group_labels = []
                min_distance = 8000000  # Increased from 5MB to 8MB minimum distance between labels

                for i, (pos, label) in enumerate(zip(gene_positions, gene_labels)):
                    if not current_group_pos or (pos - current_group_pos[-1] >= min_distance):
                        # Start a new group or add to current if far enough
                        if current_group_pos:
                            # Finalize previous group
                            avg_pos = sum(current_group_pos) // len(current_group_pos)
                            combined_label = '/'.join(current_group_labels)
                            grouped_annotations.append((avg_pos, combined_label))

                        current_group_pos = [pos]
                        current_group_labels = [label]
                    else:
                        # Add to current group (genes too close together)
                        current_group_pos.append(pos)
                        if label not in current_group_labels:  # Avoid duplicates
                            current_group_labels.append(label)

                # Don't forget the last group
                if current_group_pos:
                    avg_pos = sum(current_group_pos) // len(current_group_pos)
                    combined_label = '/'.join(current_group_labels)
                    grouped_annotations.append((avg_pos, combined_label))

                # Use annotate method for better overlap resolution
                for pos, label in grouped_annotations:
                    # Determine if this gene/fusion is in Mitelman database
                    # Check if any gene in the combined label is part of a known fusion
                    is_known_fusion = False
                    genes_in_label = label.split('/')

                    for gene in genes_in_label:
                        # Find all fusions in this chromosome that contain this gene
                        gene_fusions = bnd_chr[bnd_chr['GENE'] == gene]['FUSION'].tolist()
                        for fusion in gene_fusions:
                            if fusion in unique_mitelman_fusion:
                                is_known_fusion = True
                                break
                        if is_known_fusion:
                            break

                    # Set text color based on whether it's a known fusion
                    text_color = "red" if is_known_fusion else "black"

                    bnd_track.annotate(
                        x=pos,
                        label=label,
                        min_r=None,  # defaults to track outer radius
                        max_r=None,  # defaults to min_r + 5
                        label_size=10,
                        line_kws=dict(color="grey", lw=0.5),
                        text_kws=dict(color=text_color),
                    )

    # Set color, grey if fusion is empty, red if fusion is in mitelman fusion,
    # light grey if other
    # Define colors for different fusion types
    def get_link_color(fusion):
        if fusion == '':
            return 'lightgrey'
        elif fusion in unique_mitelman_fusion:
            return 'red'
        else:
            return 'black'

    # Update the link coloring using the color function with transparency
    bnd_only.apply(lambda row: circos.link(
        (row['CHROM'], row['POS'], row['POS']),
        (row['CHR2'], row['POS2'], row['POS2']),
        color=get_link_color(row['FUSION']),
        alpha=args.link_alpha
    ), axis=1)
    # %%

    # %%
    fig = circos.plotfig(dpi = 300)
    fig.savefig(output_file, dpi = 300, bbox_inches = 'tight')
    fig.savefig(output_png, dpi = 300, bbox_inches = 'tight')
# %%
if __name__ == '__main__':
    main()
