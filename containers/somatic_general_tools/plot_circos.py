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
    bnd_df = df[df['SVTYPE'] == 'BND']
    # If empty, return empty DataFrame
    if bnd_df.empty:
        return bnd_df
    
    # Extract chr2 and pos2 using regex
    # ALT field format example: ]chr2:123456]N
    pattern = r'[[\]]?(chr\w+):(\d+)[[\]]'
    
    # Create new columns
    bnd_df['CHR2'] = bnd_df['ALT'].str.extract(pattern)[0]
    bnd_df['POS2'] = bnd_df['ALT'].str.extract(pattern)[1].astype(int)

    # If BCSQ column exists, extract gene name (second field split by |)
    if 'BCSQ' in bnd_df.columns:
        bnd_df['GENE'] = bnd_df['BCSQ'].str.split('|').str[1]
        # If NA, replace with empty string
        bnd_df['GENE'] = bnd_df['GENE'].fillna('')
    
    return bnd_df

def create_fusion_annotation(df):
    # Create empty FUSION column
    df['FUSION'] = ''
    
    # Filter for BND records only
    bnd_records = df[df['SVTYPE'] == 'BND'].copy()
    
    for idx, row in bnd_records.iterrows():
        # Find mate record using MATEID
        mate = df[df['ID'] == row['MATE_ID']].iloc[0] if not df[df['ID'] == row['MATE_ID']].empty else None
        
        # Check if mate exists and both have GENE information
        if mate is not None and row['GENE'] != '' and mate['GENE'] != '':
            # If both genes are the same, skip
            if row['GENE'] == mate['GENE']:
                continue
            # Sort genes alphabetically and join with ::
            genes = sorted([row['GENE'], mate['GENE']])
            fusion = '::'.join(genes)
            
            # Update FUSION column for both records
            df.loc[idx, 'FUSION'] = fusion
            df.loc[df['ID'] == row['MATE_ID'], 'FUSION'] = fusion

    # Fill empty FUSION with empty string
    df['FUSION'] = df['FUSION'].fillna('')
    
    return df

# %%

def save_empty_circos(output_file):
    circos_empty = Circos.initialize_from_bed("/app/hg38.bed", space = 3)
    circos_empty.text("No SVs found", size=15)
    fig = circos_empty.plotfig(dpi = 300)
    fig.savefig(output_file, dpi = 300, bbox_inches = 'tight')

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
    args = parser.parse_args()

    # Generate output filename by replacing .vcf or .vcf.gz with .pdf
    output_file = args.sample_name + '.pdf'
    output_png = args.sample_name + '.png'

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
    circos = Circos.initialize_from_bed("/app/hg38.bed", space = 3)

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

    # Subset df to include only chromosomes that's in cytoband #chrom column
    bnd_only = bnd_only[bnd_only['CHROM'].isin(cytoband['#chrom'])]
    bnd_only = bnd_only[bnd_only['CHR2'].isin(cytoband['#chrom'])]

    # If empty, save empty fig
    if bnd_only.empty:
        save_empty_circos(output_file)
        return

    # Plot chromosome name
    for sector in circos.sectors:
        sector.text(sector.name, size=10)
        # Annotate gene name
        bnd_chr = bnd_only[bnd_only['CHROM'] == sector.name]
        # Only label if fusion is not empty
        bnd_chr = bnd_chr[bnd_chr['FUSION'] != '']
        bnd_track = sector.get_track('cytoband')
        label_pos = bnd_chr['POS']
        labels = bnd_chr['GENE']
        bnd_track.xticks(
            label_pos,
            labels,
            label_orientation="vertical",
            outer = False,
            show_bottom_line=True,
            label_size=10,
            line_kws=dict(ec="grey"),
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

    # Update the link coloring using the color function
    bnd_only.apply(lambda row: circos.link(
        (row['CHROM'], row['POS'], row['POS']),
        (row['CHR2'], row['POS2'], row['POS2']),
        color=get_link_color(row['FUSION'])
    ), axis=1)
    # %%

    # %% 
    fig = circos.plotfig(dpi = 300)
    fig.savefig(output_file, dpi = 300, bbox_inches = 'tight')
    fig.savefig(output_png, dpi = 300, bbox_inches = 'tight')
# %%
if __name__ == '__main__':
    main()
