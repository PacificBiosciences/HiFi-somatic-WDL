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
    File output_vcf = pname + "_severus/somatic_SVs/severus_somatic_" + sub(basename(tumor_bam), "\\.bam$", "") + ".vcf"
  }

  runtime {
    docker: "quay.io/biocontainers/severus@sha256:1f7eb598ed310b7d676d9107f1cecdce918acd469c563d4cb47b8ac81cd86048"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}