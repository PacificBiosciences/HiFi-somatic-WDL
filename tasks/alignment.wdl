version 1.1

# Align using pbmm2
task Align {
  input {
    File bam_file
    File ref_fasta
    File? ref_fasta_mmi
    File ref_fasta_index
    String sample_name
    Int threads
  }

  # Use mmi if provided, else use FASTA
  File ref_to_use = select_first([ref_fasta_mmi, ref_fasta])
  String ofile_name = sub(basename(bam_file), "\\.bam$", ".aligned.bam")
  String ofile_name_index = sub(basename(bam_file), "\\.bam$", ".aligned.bam.bai")
  
  command <<<
  set -euxo pipefail

  echo "Aligning ~{bam_file} and ~{ref_fasta} for ~{sample_name} using pbmm2 into ~{ofile_name}"
  pbmm2 align \
    ~{ref_to_use} \
    ~{bam_file} \
    ~{ofile_name} \
    --sample ~{sample_name} \
    --sort -j ~{threads} \
    --unmapped \
    --log-level INFO --log-file pbmm2.log
  >>>

  output {
    File aligned_bam = ofile_name
    File? aligned_bam_index = ofile_name_index
    File? align_log = "pbmm2.log"
  }

  runtime {
    container: "quay.io/biocontainers/pbmm2:1.12.0--h9ee0642_0"
    cpu: threads
    memory: "~{threads * 4} GB"
  }
}

# Merge bams using pbmerge
task MergeBams {
  input {
    String sample_name
    Array[File] bam_files
    Int threads
  }
  
  command <<<
    set -euxo pipefail
    echo "merging ~{sep(' ', bam_files)} into ~{sample_name + ".aligned.bam"}"
    pbmerge -j ~{threads} \
      -o ~{sample_name + ".aligned.bam"} \
      ~{sep(' ', bam_files)}
  >>>

  output {
    File merged_aligned_bam = sample_name + ".aligned.bam"
    File merged_aligned_bam_index = sample_name + ".aligned.bam.pbi"
  }

  runtime {
    container: "quay.io/biocontainers/pbtk:3.1.0--h9ee0642_0"
    cpu: threads
    memory: "~{threads * 4} GB"
  }
}

task IndexBam {
  input {
    File bam
    Int threads
  }
  
  command <<<
    set -euxo pipefail
    echo "indexing ~{bam}"
    samtools index -@~{threads} ~{bam} \
      -o ~{basename(bam) + ".bai"}
  >>>

  output {
    File merged_aligned_bam_index_bai = basename(bam) + ".bai"
  }

  runtime {
    container: "quay.io/biocontainers/samtools:1.17--hd87286a_1"
    cpu: threads
    memory: "~{threads * 4} GB"
  }
}