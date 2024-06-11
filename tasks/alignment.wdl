version 1.0

# Align using pbmm2
task Align {
  input {
    File bam_file
    File ref_fasta
    File ref_fasta_index
    String sample_name
    Int threads
  }

  String ofile_name = sub(basename(bam_file), "\\.bam$", ".aligned.bam")
  String ofile_name_index = sub(basename(bam_file), "\\.bam$", ".aligned.bam.bai")
  Float file_size = ceil(size(bam_file, "GB") * 2 + size(ref_fasta, "GB") + 20)
  
  command <<<
  set -euxo pipefail

  echo "Aligning ~{bam_file} and ~{ref_fasta} for ~{sample_name} using pbmm2 into ~{ofile_name}"
  
  pbmm2 --version

  pbmm2 align \
    ~{ref_fasta} \
    ~{bam_file} \
    ~{ofile_name} \
    --sample ~{sample_name} \
    --sort -j ~{threads} \
    --unmapped \
    --preset HIFI \
    --log-level DEBUG --log-file pbmm2.log
  >>>

  output {
    File aligned_bam = ofile_name
    File? aligned_bam_index = ofile_name_index
    File? align_log = "pbmm2.log"
  }

  runtime {
    docker: "quay.io/biocontainers/pbmm2:1.13.1--h9ee0642_0"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

# Merge bams using pbmerge
task MergeBams {
  input {
    String sample_name
    Array[File] bam_files
    Int threads
  }

  Float file_size = ceil(size(bam_files, "GB") * 2 + 20)
  
  command <<<
    set -euxo pipefail
    echo "merging ~{sep=' ' bam_files} into ~{sample_name + ".aligned.bam"}"
    
    samtools --version

    samtools merge -@~{threads} \
      -o ~{sample_name + ".aligned.bam"} \
      ~{sep=' ' bam_files}

    samtools index -@~{threads} ~{sample_name + ".aligned.bam"}
  >>>

  output {
    File merged_aligned_bam = sample_name + ".aligned.bam"
    File merged_aligned_bam_index = sample_name + ".aligned.bam.bai"
  }

  runtime {
    docker: "quay.io/biocontainers/samtools:1.17--hd87286a_1"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

task IndexBam {
  input {
    File bam
    Int threads
  }

  Float file_size = ceil(size(bam, "GB") + 20)
  
  command <<<
    set -euxo pipefail
    echo "indexing ~{bam}"
    
    samtools --version

    samtools index -@~{threads} ~{bam} \
      -o ~{basename(bam) + ".bai"}
  >>>

  output {
    File merged_aligned_bam_index_bai = basename(bam) + ".bai"
  }

  runtime {
    docker: "quay.io/biocontainers/samtools:1.17--hd87286a_1"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}
