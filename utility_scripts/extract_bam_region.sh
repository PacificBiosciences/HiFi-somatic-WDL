#!/bin/bash

set -euo pipefail

# This script is used to extract region from all bam files in the WDL
# workflow directory. The workflow directory is passed as the first argument.
# The second argument is the region to extract. The third argument is the
# output directory. The fourth argument is the output file name prefix.
# Requires samtools, fd (fd-find, install via Cargo or bioconda), bcftools
# and bedtools

# Run by (Region can have comma or no comma):
# ./extract_bam_region.sh WORKFLOW_DIR chr1:1,000-2,000 OUTPUT_DIR OUTPUT_PREFIX THREADS
# Example command to extract EGFR:
# bash ~/pb_bitbucket/wdl-hifisomatic/utility_scripts/extract_bam_region.sh ./ chr7:55,019,017-55,211,628 ./test_vis EGFR 16

workdir=$1
region=$2
outdir=$3
mkdir -p "$outdir"
prefix=$4
threads=$5

# Split region into chr start and end. Remove comma from start and end if present
chr=$(echo "$region" | cut -d':' -f1)
start=$(echo "$region" | cut -d':' -f2 | cut -d'-' -f1 | sed 's/,//g')
end=$(echo "$region" | cut -d':' -f2 | cut -d'-' -f2 | sed 's/,//g')

echo -e "${chr}\t${start}\t${end}" > ${prefix}_${chr}_${start}-${end}.bed

# In workflow directory, look for */_LAST/out/ for bam files
tumor_bam=$(fd -p -j "${threads}" --type symlink 'tumor_bams_phased/.*\.bam$' "${workdir}"/*/_LAST/out/ -x readlink -f {})
normal_bam=$(fd -p -j "${threads}" --type symlink 'normal_bams_phased/.*\.bam$' "${workdir}"/*/_LAST/out/ -x readlink -f {})
vcfs=$(fd -p -j "${threads}" --type symlink 'small_variant_vcf/.*\.vcf.gz$' "${workdir}"/*/_LAST/out/ -x readlink -f {})
sniffles_vcfs=$(fd -p -j "${threads}" --type symlink 'sniffles_somatic_vcf_filterHPRC/.*\.vcf.gz$' "${workdir}"/*/_LAST/out/ -x readlink -f {})
severus_vcfs=$(fd -p -j "${threads}" --type symlink 'Severus_filterHPRC_vcf/.*\.vcf.gz$' "${workdir}"/*/_LAST/out/ -x readlink -f {})
combined_cpg_pileup_normal=$(fd -p -j "${threads}" --type symlink 'pileup_normal_bed/.*\.normal\.cpg\.combined\.bed$' "${workdir}"/*/_LAST/out/ -x readlink -f {})
combined_cpg_pileup_tumor=$(fd -p -j "${threads}" --type symlink 'pileup_tumor_bed/.*\.tumor\.cpg\.combined\.bed$' "${workdir}"/*/_LAST/out/ -x readlink -f {})
hap1_cpg_pileup_normal=$(fd -p -j "${threads}" --type symlink 'pileup_normal_bed/.*\.normal\.cpg\.hap1\.bed$' "${workdir}"/*/_LAST/out/ -x readlink -f {})
hap1_cpg_pileup_tumor=$(fd -p -j "${threads}" --type symlink 'pileup_tumor_bed/.*\.tumor\.cpg\.hap1\.bed$' "${workdir}"/*/_LAST/out/ -x readlink -f {})
hap2_cpg_pileup_normal=$(fd -p -j "${threads}" --type symlink 'pileup_normal_bed/.*\.normal\.cpg\.hap2\.bed$' "${workdir}"/*/_LAST/out/ -x readlink -f {})
hap2_cpg_pileup_tumor=$(fd -p -j "${threads}" --type symlink 'pileup_tumor_bed/.*\.tumor\.cpg\.hap2\.bed$' "${workdir}"/*/_LAST/out/ -x readlink -f {})

# Extract region from each bam file and output to outdir
# Convert tumor_bam to array first
tumor_bam=($tumor_bam)
normal_bam=($normal_bam)
vcfs=($vcfs)
sniffles_vcfs=($sniffles_vcfs)
severus_vcfs=($severus_vcfs)
combined_cpg_pileup_normal=($combined_cpg_pileup_normal)
combined_cpg_pileup_tumor=($combined_cpg_pileup_tumor)
hap1_cpg_pileup_normal=($hap1_cpg_pileup_normal)
hap1_cpg_pileup_tumor=($hap1_cpg_pileup_tumor)
hap2_cpg_pileup_normal=($hap2_cpg_pileup_normal)
hap2_cpg_pileup_tumor=($hap2_cpg_pileup_tumor)
for ((i=0;i<${#tumor_bam[@]};++i)); do
    sample=$(basename "${tumor_bam[$i]}" | sed 's/.tumor.aligned.hiphase.bam//g')
    # Extract BAM region
    mkdir -p "${outdir}/${sample}"/BAM/
    mkdir -p "${outdir}/${sample}"/VCF/
    mkdir -p "${outdir}/${sample}"/pileup/
    echo -e "Extracting tumor BAM for ${sample}"
    samtools view -@"${threads}" -bh "${tumor_bam[$i]}" "$region" > "${outdir}/${sample}/BAM/${sample}.tumor_${prefix}.bam"
    samtools index -@"${threads}" "${outdir}/${sample}/BAM/${sample}.tumor_${prefix}.bam"
    # Find normal from list of normal_bam
    echo "${normal_bam[@]}"
    normal=$(echo "${normal_bam[@]}" | tr ' ' '\n' | grep "${sample}\.")
    echo -e "Extracting normal BAM for ${sample}"
    samtools view -@"${threads}" -bh "${normal}" "$region" > "${outdir}/${sample}/BAM/${sample}.normal_${prefix}.bam"
    samtools index -@"${threads}" "${outdir}/${sample}/BAM/${sample}.normal_${prefix}.bam"

    # Find all SNV/INDEL VCF for the sample and 
    # extract according to region desired
    echo "${vcfs[@]}" | tr ' ' '\n' | grep "${sample}\." > "${outdir}/${sample}/VCF/vcf.list"
    echo -e "Extracting VCF for ${sample}"
    while IFS= read -r vcf
    do
        vcf_name=$(basename "${vcf}")
        # Index if tbi does not exist
        if [[ ! -f "${vcf}.tbi" ]]
        then
            tabix -p vcf "${vcf}"
        fi
        bcftools view ${vcf} \
            -R ${prefix}_${chr}_${start}-${end}.bed \
            -Oz -o "${outdir}/${sample}/VCF/${vcf_name%.vcf.gz}_${prefix}.vcf.gz"
        tabix -p vcf "${outdir}/${sample}/VCF/${vcf_name%.vcf.gz}_${prefix}.vcf.gz"
    done < <(cat "${outdir}"/"${sample}"/VCF/vcf.list)

    # Do the same for SV (Severus and Sniffles)
    sniffles_vcf=$(echo "${sniffles_vcfs[@]}" | tr ' ' '\n' | grep "${sample}\.sniffles")
    severus_vcf=$(echo "${severus_vcfs[@]}" | tr ' ' '\n' | grep "severus_somatic_${sample}\.")
    sniffles_vcf_name=$(basename "${sniffles_vcf}")
    severus_vcf_name=$(basename "${severus_vcf}")
    echo -e "Extracting Sniffles VCF for ${sample}"
    bcftools view ${sniffles_vcf} \
        -R ${prefix}_${chr}_${start}-${end}.bed \
        -Oz -o "${outdir}/${sample}/VCF/${sniffles_vcf_name%.vcf.gz}_${prefix}.vcf.gz"
    tabix -p vcf "${outdir}/${sample}/VCF/${sniffles_vcf_name%.vcf.gz}_${prefix}.vcf.gz"
    echo -e "Extracting Severus VCF for ${sample}"
    bcftools view ${severus_vcf}\
        -R ${prefix}_${chr}_${start}-${end}.bed \
        -Oz -o "${outdir}/${sample}/VCF/${severus_vcf_name%.vcf.gz}_${prefix}.vcf.gz"
    tabix -p vcf "${outdir}/${sample}/VCF/${severus_vcf_name%.vcf.gz}_${prefix}.vcf.gz"

    # Extract pileup
    comb_normal=$(echo "${combined_cpg_pileup_normal[@]}" | tr ' ' '\n' | grep "${sample}\.")
    comb_tumor=$(echo "${combined_cpg_pileup_tumor[@]}" | tr ' ' '\n' | grep "${sample}\.")
    hap1_normal=$(echo "${hap1_cpg_pileup_normal[@]}" | tr ' ' '\n' | grep "${sample}\.")
    hap1_tumor=$(echo "${hap1_cpg_pileup_tumor[@]}" | tr ' ' '\n' | grep "${sample}\.")
    hap2_normal=$(echo "${hap2_cpg_pileup_normal[@]}" | tr ' ' '\n' | grep "${sample}\.")
    hap2_tumor=$(echo "${hap2_cpg_pileup_tumor[@]}" | tr ' ' '\n' | grep "${sample}\.")
    
    bedtools intersect -a "${comb_tumor}" -b ${prefix}_${chr}_${start}-${end}.bed > "${outdir}/${sample}/pileup/${sample}.tumor.cpg.combined_${prefix}.bed"
    bedtools intersect -a "${comb_normal}" -b ${prefix}_${chr}_${start}-${end}.bed > "${outdir}/${sample}/pileup/${sample}.normal.cpg.combined_${prefix}.bed"
    # Extract hap1 hap2 etc only if the variable is not empty
    if [[ -n "${hap1_normal}" ]]
    then
        echo -e "Extracting pileup tumor hap 1 for ${sample}"
        # echo -e "bedtools intersect -a ${hap1_tumor} -b ${prefix}_${chr}_${start}-${end}.bed > ${outdir}/${sample}/pileup/${sample}.tumor.cpg.hap1_${prefix}.bed"
        bedtools intersect -a "${hap1_tumor}" -b ${prefix}_${chr}_${start}-${end}.bed > "${outdir}"/"${sample}"/pileup/${sample}.tumor.cpg.hap1_${prefix}.bed
        echo -e "Extracting pileup normal hap 1 for ${sample}"
        # echo -e "bedtools intersect -a ${hap1_normal} -b ${prefix}_${chr}_${start}-${end}.bed > ${outdir}/${sample}/pileup/${sample}.normal.cpg.hap1_${prefix}.bed"
        bedtools intersect -a "${hap1_normal}" -b ${prefix}_${chr}_${start}-${end}.bed > "${outdir}"/"${sample}"/pileup/${sample}.normal.cpg.hap1_${prefix}.bed
        echo -e "Extracting pileup tumor hap 2 for ${sample}"
        # echo -e "bedtools intersect -a ${hap2_tumor} -b ${prefix}_${chr}_${start}-${end}.bed > ${outdir}/${sample}/pileup/${sample}.tumor.cpg.hap2_${prefix}.bed"
        bedtools intersect -a "${hap2_tumor}" -b ${prefix}_${chr}_${start}-${end}.bed > "${outdir}"/"${sample}"/pileup/${sample}.tumor.cpg.hap2_${prefix}.bed
        echo -e "Extracting pileup normal hap 2 for ${sample}"
        # echo -e "bedtools intersect -a ${hap2_normal} -b ${prefix}_${chr}_${start}-${end}.bed > ${outdir}/${sample}/pileup/${sample}.normal.cpg.hap2_${prefix}.bed"
        bedtools intersect -a "${hap2_normal}" -b ${prefix}_${chr}_${start}-${end}.bed > "${outdir}"/"${sample}"/pileup/${sample}.normal.cpg.hap2_${prefix}.bed
    fi
done