version 1.0

# Use Severus to call SV
task Severus_sv {
  input {
    File tumor_bam
    File tumor_bam_index
    File normal_bam
    File normal_bam_index
    File trf_bed
    File? phased_vcf
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
      --control-bam ~{normal_bam} \
      ~{"--phasing-vcf " + phased_vcf} \
      --out-dir ~{pname + "_severus"} \
      -t ~{threads} \
      --vntr-bed ~{trf_bed} \
      --min-support ~{min_supp_reads}
  >>>

  output {
    File output_vcf = pname + "_severus/somatic_SVs/severus_somatic" + ".vcf"
  }

  runtime {
    docker: "quay.io/biocontainers/severus@sha256:5fe46f34a1b39a5d6d9027d7b55723716fedee0c2f9e98ca8ba366f76c67436d"
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