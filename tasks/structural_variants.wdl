version 1.0

# Use Severus to call SV
task Severus_sv {
  input {
    File tumor_bam
    File tumor_bam_index
    File? normal_bam
    File? normal_bam_index
    File trf_bed
    File? phased_vcf
    File? PON_tsv
    String pname
    Int threads
    Int min_supp_reads
  }

  Float file_size = ceil(size(tumor_bam, "GB") + size(normal_bam, "GB") + size(phased_vcf, "GB") + 10)

  command <<<
    set -euxo pipefail
    
    echo "Running Severus for ~{pname}"

    severus --version

    severus \
      --target-bam ~{tumor_bam} \
      ~{"--control-bam " + normal_bam} \
      ~{"--phasing-vcf " + phased_vcf} \
      ~{"--PON " + PON_tsv} \
      --out-dir ~{pname + "_severus"} \
      -t ~{threads} \
      --vntr-bed ~{trf_bed} \
      --min-support ~{min_supp_reads} \
      --resolve-overlaps \
      --between-junction-ins \
      --single-bp

    # Compress SVs plots HTML inside somatic_SVs/plots directory
    # Check if the directory exists first
    if [[ -d ~{pname + "_severus/somatic_SVs/plots"} ]]
      then tar -czvf ~{pname + "_severus/somatic_SVs/" + pname + ".plots.tar.gz"} ~{pname + "_severus/somatic_SVs/plots"}
    fi

    # Rename output vcf
    if [[ -f ~{pname + "_severus/somatic_SVs/severus_somatic.vcf"} ]]; then mv ~{pname + "_severus/somatic_SVs/severus_somatic.vcf"} ~{pname + "_severus/somatic_SVs/" + pname + ".severus_somatic.vcf"}; fi
    if [[ -f ~{pname + "_severus/all_SVs/severus_all.vcf"} ]]; then mv ~{pname + "_severus/all_SVs/severus_all.vcf"} ~{pname + "_severus/all_SVs/" + pname + ".severus_all.vcf"}; fi
    # Rename breakpoint clusters list
    if [[ -f ~{pname + "_severus/somatic_SVs/breakpoint_clusters_list.tsv"} ]]; then mv ~{pname + "_severus/somatic_SVs/breakpoint_clusters_list.tsv"} ~{pname + "_severus/somatic_SVs/" + pname + ".somatic_breakpoint_clusters_list.tsv"}; fi
    if [[ -f ~{pname + "_severus/all_SVs/breakpoint_clusters_list.tsv"} ]]; then mv ~{pname + "_severus/all_SVs/breakpoint_clusters_list.tsv"} ~{pname + "_severus/all_SVs/" + pname + ".all_breakpoint_clusters_list.tsv"}; fi
  >>>

  output {
    File output_vcf = pname + "_severus/somatic_SVs/" + pname + ".severus_somatic.vcf"
    File? output_all_vcf = pname + "_severus/all_SVs/" + pname + ".severus_all.vcf"
    File? output_breakpoint_clusters = pname + "_severus/somatic_SVs/" + pname + ".somatic_breakpoint_clusters_list.tsv"
    File? output_breakpoint_clusters_all = pname + "_severus/all_SVs/" + pname + ".all_breakpoint_clusters_list.tsv"
    File? output_somatic_sv_plots = pname + "_severus/somatic_SVs/" + pname + ".plots.tar.gz"
  }

  runtime {
    docker: "quay.io/biocontainers/severus@sha256:df04bec4a0ae9c55104ff91d6063d8c7c58b05145db04ad481f5d3d9527a9b7d"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

# SAVANA SV calling
task SAVANA_sv {
  input {
    File tumor_bam
    File tumor_bam_index
    File normal_bam
    File normal_bam_index
    File ref_fasta
    File ref_fasta_index
    String pname
    Int threads
    Int? min_supp_reads
    Int svlength = 50
    File? phased_vcf
  }

  command <<<
    set -euxo pipefail

    echo -e "chr1\nchr2\nchr3\nchr4\nchr5\nchr6\nchr7\nchr8\nchr9\nchr10\nchr11\nchr12\nchr13\nchr14\nchr15\nchr16\nchr17\nchr18\nchr19\nchr20\nchr21\nchr22\nchrX\nchrY" > contig.txt
    
    echo "Running SAVANA for ~{pname}"
    savana \
      --tumour ~{tumor_bam} \
      --normal ~{normal_bam} \
      --ref ~{ref_fasta} \
      --threads ~{threads} \
      --outdir ~{pname + "_savana"} \
      --contigs contig.txt \
      --single_bnd \
      --sample ~{pname} \
      --length ~{svlength} \
      --pb \
      ~{"--snp_vcf " + phased_vcf} \
      ~{"--min_support " + min_supp_reads}
  >>>

  output {
    File savana_vcf = pname + "_savana/" + pname + ".sv_breakpoints.vcf"
    File savana_classified_somatic_vcf = pname + "_savana/" + pname + ".classified.somatic.vcf"
    File savana_classified_vcf = pname + "_savana/" + pname + ".classified.vcf"
    File savana_sv_breakpoints_bedpe = pname + "_savana/" + pname + ".sv_breakpoints.bedpe"
    File savana_sv_breakpoints_read_support = pname + "_savana/" + pname + ".sv_breakpoints_read_support.tsv"
    File savana_10kbp_bin_ref_all = pname + "_savana/10kbp_bin_ref_all_" + pname + "_with_SV_breakpoints.bed"
    File savana_allele_counts_hetSNPs = pname + "_savana/" + pname + "_allele_counts_hetSNPs.bed"
    File savana_fitted_purity_ploidy = pname + "_savana/" + pname + "_fitted_purity_ploidy.tsv"
    File savana_inserted_sequences = pname + "_savana/" + pname + ".inserted_sequences.fa"
    File savana_ranked_solutions = pname + "_savana/" + pname + "_ranked_solutions.tsv"
    File savana_raw_read_counts = pname + "_savana/" + pname + "_raw_read_counts.tsv"
    File savana_read_counts_mnorm_log2r_segmented = pname + "_savana/" + pname + "_read_counts_mnorm_log2r_segmented.tsv"
    File savana_segmented_absolute_copy_number = pname + "_savana/" + pname + "_segmented_absolute_copy_number.tsv"
  }

  runtime {
    container: "quay.io/biocontainers/savana@sha256:adaf7806939fae6f656187a137e362b811780e8b420e83e4381657a6f3a40d66"
    cpu: threads
    memory: "~{threads * 4} GB"
  }
}

# Take SV VCF and do circos plot of BND events
task circos_BND {
  input {
    File sv_vcf
    Int max_gene_labels = 100
    String pname

    Int threads = 8
  }

  Float file_size = ceil(size(sv_vcf, "GB") + 10)

  command <<<
    set -euxo pipefail

    echo "Running Circos for ~{pname}"

    python /app/plot_circos.py \
      ~{sv_vcf} \
      ~{pname + ".circos"} \
      /app/MCGENE.TXT.DATA \
      --max-gene-labels ~{max_gene_labels}
  >>>

  output {
    File circos_plot = pname + ".circos.pdf"
    File circos_png = pname + ".circos.png"
    File mitelman_fusions = pname + ".circos.mitelman_fusions.tsv"
    File known_genes_pairs = pname + ".circos.known_gene_pairs.tsv"
  }

  runtime {
      docker: "quay.io/pacbio/somatic_general_tools@sha256:99159e099d9044c52debabdc9491b168487aaa37534c1a748809bc69f169812a"
      cpu: threads
      memory: "~{threads * 4} GB"
      disk: file_size + " GB"
      maxRetries: 2
      preemptible: 1
  }
}

task wakhan {
  input {
    String pname
    File tumor_bam
    File tumor_bam_index
    File ref_fasta
    File ref_fasta_index
    File severus_sv_vcf
    File? normal_germline_vcf
    File tumor_germline_vcf
    Int threads = 16
    String purity_range = "0.2-0.99"
    String ploidy_range = "1.0-6.0"
  }

  Float file_size = ceil(size(tumor_bam, "GB") + size(ref_fasta, "GB") + size(severus_sv_vcf, "GB") + size(normal_germline_vcf, "GB") + 10)

  command <<<
    set -euxo pipefail

    echo "Running Wakhan"

    wakhan \
      --threads ~{threads} \
      --target-bam ~{tumor_bam} \
      --reference ~{ref_fasta} \
      --genome-name ~{pname} \
      --out-dir ~{pname + "_wakhan"} \
      --breakpoints ~{severus_sv_vcf} \
      --loh-enable \
      --use-sv-haplotypes \
      --ploidy-range ~{ploidy_range} \
      --purity-range ~{purity_range} \
      --centromere /opt/wakhan/src/annotations/grch38.cen_coord.curated.bed \
      ~{if (defined(normal_germline_vcf)) then " --normal-phased-vcf " + normal_germline_vcf else " --tumor-phased-vcf " + tumor_germline_vcf + " --hets-ratio 0.25"}

    # Save purity and ploidy
    cd ~{pname + "_wakhan"}
    # Create TSV file with headers
    echo -e "folder_name\tploidy\tpurity\tconfidence" > folder_numbers.tsv

    # Find directories matching the pattern and process them
    find . -type d -regex ".*/[0-9.]+_[0-9.]+_[0-9.]+$" | while read dir; do
        # Get just the folder name without the path
        folder_name=$(basename "$dir")
        
        # Split the folder name into components
        ploidy=$(echo $folder_name | cut -d'_' -f1)
        purity=$(echo $folder_name | cut -d'_' -f2)
        confidence=$(echo $folder_name | cut -d'_' -f3)
        
        # Append to a temporary file
        echo -e "$folder_name\t$ploidy\t$purity\t$confidence"
    done | sort -t$'\t' -k4,4rn > temp.tsv

    # Combine header and sorted content
    (head -n 1 folder_numbers.tsv; cat temp.tsv) > ../purity_ploidy.tsv
    rm temp.tsv folder_numbers.tsv

    cd ..

    # Compress all output
    tar -czvf ~{pname + "_wakhan.tar.gz"} ~{pname + "_wakhan"}

    # Keep the bed output of the best solution
    best_folder=$(head -n 2 purity_ploidy.tsv | tail -n 1 | cut -f1)
    cp -r ~{pname + "_wakhan"}/${best_folder} ~{pname + "_wakhan_best"}
    # Delete variation_plots folder (big!)
    rm -rf ~{pname + "_wakhan_best/variation_plots"}
    # Rename files
    # if there's copy number segments ending with HP_1.bed, rename to HP_1.bed. Also means
    # there must be HP_2.bed as well.
    if [[ -n $(compgen -G "~{pname + "_wakhan_best"}/bed_output/*copynumbers_segments_HP_1.bed") ]]; then
      mv ~{pname + "_wakhan_best"}/bed_output/*copynumbers_segments_HP_1.bed ~{pname + ".copynumbers_segments_HP_1.bed"}
      mv ~{pname + "_wakhan_best"}/bed_output/*copynumbers_segments_HP_2.bed ~{pname + ".copynumbers_segments_HP_2.bed"}
      # Move the subclonal segments
      mv ~{pname + "_wakhan_best"}/bed_output/*subclonal_segments_HP_1.bed ~{pname + ".subclonal_segments_HP_1.bed"}
      mv ~{pname + "_wakhan_best"}/bed_output/*subclonal_segments_HP_2.bed ~{pname + ".subclonal_segments_HP_2.bed"}
    else
      mv ~{pname + "_wakhan_best"}/bed_output/*copynumbers_segments.bed ~{pname + ".copynumbers_segments.bed"}
      mv ~{pname + "_wakhan_best"}/bed_output/*subclonal_segments.bed ~{pname + ".subclonal_segments.bed"}
    fi
    mv ~{pname + "_wakhan_best"}/bed_output/loh_regions.bed ~{pname + ".loh_regions.bed"}
    mv ~{pname + "_wakhan_best"}/bed_output/genes_copynumber_states.bed ~{pname + ".cancer_genes_copynumber.bed"}
    # HTML plot
    mv ~{pname + "_wakhan_best"}/*genome_copynumbers_breakpoints.html ~{pname + ".genome_copynumbers_breakpoints.html"}
    mv ~{pname + "_wakhan_best"}/*genome_copynumbers_details.html ~{pname + ".genome_copynumbers_details.html"}

    rm -rf ~{pname + "_wakhan"}
  >>>

  output {
    File wakhan_tar = pname + "_wakhan.tar.gz"
    File wakhan_purity_ploidy = "purity_ploidy.tsv"
    File? copynumbers_segments = pname + ".copynumbers_segments.bed"
    File? copynumbers_segments_HP_1 = pname + ".copynumbers_segments_HP_1.bed"
    File? copynumbers_segments_HP_2 = pname + ".copynumbers_segments_HP_2.bed"
    File? subclonal_segments_HP_1 = pname + ".subclonal_segments_HP_1.bed"
    File? subclonal_segments_HP_2 = pname + ".subclonal_segments_HP_2.bed"
    File? subclonal_segments = pname + ".subclonal_segments.bed"
    File loh_regions = pname + ".loh_regions.bed"
    File cancer_genes_copynumber = pname + ".cancer_genes_copynumber.bed"
    File genome_copynumbers_breakpoints = pname + ".genome_copynumbers_breakpoints.html"
    File genome_copynumbers_details = pname + ".genome_copynumbers_details.html"
    Array[File] wakhan_output = select_all([wakhan_tar, wakhan_purity_ploidy, copynumbers_segments, copynumbers_segments_HP_1, copynumbers_segments_HP_2, subclonal_segments, subclonal_segments_HP_1, subclonal_segments_HP_2, loh_regions, cancer_genes_copynumber, genome_copynumbers_breakpoints, genome_copynumbers_details])
  }

  runtime {
    docker: "kpinpb/wakhan@sha256:bbd2eb84e203410c7e9efc476ad1940166fc7f0769d93a0d11fbd7338afda847"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}
