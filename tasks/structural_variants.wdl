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
      then tar -czvf ~{pname + "_severus/somatic_SVs/plots.tar.gz"} ~{pname + "_severus/somatic_SVs/plots"}
    fi
  >>>

  output {
    File? output_vcf = pname + "_severus/somatic_SVs/severus_somatic" + ".vcf"
    File? output_all_vcf = pname + "_severus/all_SVs/severus_all.vcf"
    File? output_breakpoint_clusters = pname + "_severus/somatic_SVs/" + "breakpoint_clusters_list.tsv"
    File? output_breakpoint_clusters_all = pname + "_severus/all_SVs/" + "breakpoint_clusters_list.tsv"
    File? output_somatic_sv_plots = pname + "_severus/somatic_SVs/plots.tar.gz"
  }

  runtime {
    docker: "quay.io/biocontainers/severus@sha256:9aaf81a02cc7f8be5dd9bb3bece00b9097f08368cebce297aba9aa7fc932b42e"
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
    File? phased_vcf
    String pname
    Int threads
    Int? min_supp_reads
    Float? min_af
    Int svlength = 50
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
      --pb \
      --single_bnd \
      --no_blacklist \
      --sample ~{pname} \
      --length ~{svlength} \
      ~{"--phased_vcf " + phased_vcf} \
      ~{"--min_support " + min_supp_reads} \
      ~{"--min_af " + min_af}

  >>>

  output {
    File savana_vcf = pname + "_savana/" + pname + ".classified.somatic.vcf"
    File cnv_segments = pname + "_savana/" + pname + "_read_counts_mnorm_log2r_segmented.tsv"
    File cnv_absolute_cn = pname + "_savana/" + pname + "_segmented_absolute_copy_number.tsv"
    File purity_ploidy = pname + "_savana/" + pname + "_fitted_purity_ploidy.tsv"
    File purity_ploidy_solutions = pname + "_savana/" + pname + "_ranked_solutions.tsv"
    Array[File] savana_output = [savana_vcf, cnv_segments, cnv_absolute_cn, purity_ploidy, purity_ploidy_solutions]

  }

  runtime {
    container: "quay.io/biocontainers/savana@sha256:e7560628f7bd3a212520fbd5005589b78e9d95f95255867b42d143dc67a0028b"
    cpu: threads
    memory: "~{threads * 4} GB"
  }
}

# Take SV VCF and do circos plot of BND events
task circos_BND {
  input {
    File sv_vcf
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
      /app/MCGENE.TXT.DATA
  >>>

  output {
    File circos_plot = pname + ".circos.pdf"
    File circos_png = pname + ".circos.png"
  }

  runtime {
      docker: "quay.io/pacbio/somatic_general_tools@sha256:a25a2e62b88c73fa3c18a0297654420a4675224eb0cf39fa4192f8a1e92b30d6"
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
    String purity_range = "0.2-1.0"
    String ploidy_range = "1-6"
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
      --ploidy-range ~{ploidy_range} \
      --purity-range ~{purity_range} \
      --centromere /opt/wakhan/Wakhan/src/annotations/grch38.cen_coord.curated.bed \
      ~{if (defined(normal_germline_vcf)) then " --normal-phased-vcf " + normal_germline_vcf else " --tumor-vcf " + tumor_germline_vcf + " --hets-ratio 0.25"}

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
    mv ~{pname + "_wakhan_best"}/bed_output/*copynumbers_segments.bed ~{pname + ".copynumbers_segments.bed"}
    mv ~{pname + "_wakhan_best"}/bed_output/loh_regions.bed ~{pname + ".loh_regions.bed"}
    mv ~{pname + "_wakhan_best"}/bed_output/cancer_genes_copynumber_states.bed ~{pname + ".cancer_genes_copynumber.bed"}

    rm -rf ~{pname + "_wakhan"}
  >>>

  output {
    File wakhan_tar = pname + "_wakhan.tar.gz"
    File wakhan_purity_ploidy = "purity_ploidy.tsv"
    File copynumbers_segments = pname + ".copynumbers_segments.bed"
    File loh_regions = pname + ".loh_regions.bed"
    File cancer_genes_copynumber = pname + ".cancer_genes_copynumber.bed"
    Array[File] wakhan_output = [wakhan_tar, wakhan_purity_ploidy, copynumbers_segments, loh_regions, cancer_genes_copynumber]
  }

  runtime {
    docker: "mkolmogo/wakhan:dev_c717baa"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}
