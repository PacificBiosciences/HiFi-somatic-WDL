#!/bin/bash
set -euo pipefail

# This script takes in a CSV file containing
# tumor and normal pair for each patient,
# then look for the BAM files in data_dir and 
# generate a json for WDL workflow

# Usage: ./generate_wdl_json.sh <input_csv> <data_dir>

# Example input csv:
# sample_group,Normal,Tumor
# patient1,normal1,tumor1
# This expects normal1 and tumor1 to have their own folder in data_dir, e.g.
# The directory structure:
# data_dir
# ├── normal1
# │   ├── m54002_210708_210708.bam
# |── tumor1
# |   ├── m54002_210708_222312.bam
# Please change the parameters specified in the block of json text below!

# Once the jsons are generated, an example of batch slurm submission:
# for i in sample_json/*.json; do fname=$(basename ${i}); sample=$(echo ${fname} | perl -pe 's/input\.(.*)\.json/$1/g'); quick_sub miniwdl_hifisomatic_${sample} 4 "conda activate miniwdl && miniwdl run ~/pb_bitbucket/wdl-hifisomatic/hifisomatic.wdl -i ${i} -d hifisomatic_${sample} --cfg miniwdl.cfg"; done

# Author: Khi Pin, Chua
# Date: 2023-8-15

csv=$1
data_dir=$(readlink -f $2)

mkdir -p sample_json

# Split CSV. First column is patient, second column is normal, third column is tumor
for line in $(cat $csv | tail -n+2)
do
split_line=$(echo $line | csvtk csv2tab)
patient=$(echo $split_line | cut -f1 -d' ')
normal=$(echo $split_line | cut -f2 -d ' ')
tumor=$(echo $split_line | cut -f3 -d ' ')
echo "Generating json for $patient"
# Generate block for tumor
tumor_bam=$(find $data_dir/$tumor/ -regextype posix-extended -regex ".*m[568].*_[0-9]{6}_[0-9]{6}.*.bam" | sed 's/^/\ \ \ \ \ \ \ \ \ \ "/g' | sed 's/$/",/g')
if [[ -z $tumor_bam ]]
then
echo "No tumor bam found for $patient"
exit 1
fi
normal_bam=$(find $data_dir/$normal/ -regextype posix-extended -regex ".*m[568].*_[0-9]{6}_[0-9]{6}.*.bam" | sed 's/^/\ \ \ \ \ \ \ \ \ \ "/g' | sed 's/$/",/g')
if [[ -z $normal_bam ]]
then
echo "No normal bam found for $patient"
exit 1
fi

rm -f sample_json/input.${patient}.json

echo -ne "{
  \"hifisomatic.cohort\": {
    \"patients\": [\n" >> sample_json/input.${patient}.json

echo -ne "      {
        \"patient_names\": \"$patient\",\n" >> sample_json/input.${patient}.json

echo -ne "        \"tumor_bams\": [
$tumor_bam
        ],
        \"normal_bams\": [
$normal_bam
        ]
      }
    ]
  }," >> sample_json/input.${patient}.json

echo -ne "
  \"hifisomatic.ref_fasta\": \"/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta\",
  \"hifisomatic.ref_fasta_index\": \"/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.fai\",
  \"hifisomatic.ref_bed\": \"/path/to/chr.bed\",
  \"hifisomatic.hprc_sniffles_control_vcf\": \"/path/to/HPRC_sniffles207.vcf.gz\",
  \"hifisomatic.hprc_sniffles_control_vcf_index\": \"/path/to/HPRC_sniffles207.vcf.gz.tbi\",
  \"hifisomatic.skip_align\": false,
  \"hifisomatic.cnvkit_refflat\": \"/path/to/refFlat.hg38.txt\",
  \"hifisomatic.cnvkit_threads\": 16,
  \"hifisomatic.pbmm2_threads\": 24,
  \"hifisomatic.merge_bam_threads\": 8,
  \"hifisomatic.samtools_threads\": 8,
  \"hifisomatic.clairs_threads\": 16,
  \"hifisomatic.clairs_platform\": \"hifi_revio\",
  \"hifisomatic.sniffles_threads\": 16,
  \"hifisomatic.sniffles_trf_bed\": \"/path/to/human_GRCh38_no_alt_analysis_set.trf.bed\",
  \"hifisomatic.call_small_variants\": true,
  \"hifisomatic.strip_kinetics\": false,
  \"hifisomatic.normal_pileup_mincov\": 10,
  \"hifisomatic.vep_cache\": \"/path/to/homo_sapiens_refseq_vep_110_GRCh38.tar.gz\",
  \"hifisomatic.annotsv_cache\": \"/path/to/annotsv.tar.gz\"
}" >> sample_json/input.${patient}.json

rm -f tmp.tumor.json tmp.normal.json
done
