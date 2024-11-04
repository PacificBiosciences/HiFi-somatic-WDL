version 1.0

workflow align_all_bams {

    input {
        String patient
        Array[File] patient_tumor_bam_files
        Array[File] patient_normal_bam_files
        File ref_fasta
        File ref_fasta_index
        Int pbmm2_threads = 64
        Int merge_bam_threads = 8
        Int samtools_threads = 8
        Boolean skip_align = false
        Boolean strip_kinetics = false
        String additional_pbmm2_args = "-A 2"
    }

     if(!skip_align){
      scatter (tumor_bam in patient_tumor_bam_files) {

        call Align as TumorAlign {
          input:
            sample_name = patient + ".tumor",
            bam_file = tumor_bam,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            additional_args = additional_pbmm2_args,
            strip_kinetics = strip_kinetics,
            threads = pbmm2_threads
        }
      }

      scatter (normal_bam in patient_normal_bam_files) {
        call Align as NormalAlign {
          input:
            sample_name = patient + ".normal",
            bam_file = normal_bam,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            additional_args = additional_pbmm2_args,
            strip_kinetics = strip_kinetics,
            threads = pbmm2_threads
        }
      }
    }

    # If more than one BAM, merge them
      if (!skip_align && length(patient_tumor_bam_files) > 1) {
        call MergeBams as MergeTumorAlignBams {
          input:
            sample_name = patient + ".tumor",
            bam_files = select_first([TumorAlign.aligned_bam]),
            threads = merge_bam_threads
        }
      }

      if (!skip_align && length(patient_normal_bam_files) > 1) {
        call MergeBams as MergeNormalAlignBams {
          input:
            sample_name = patient + ".normal",
            bam_files = select_first([NormalAlign.aligned_bam]),
            threads = merge_bam_threads
        }
      }

      # If skip align, merge if more than one BAM
      if (skip_align && length(patient_tumor_bam_files) > 1) {
        call MergeBams as MergeTumorSkipAlignBams {
          input:
            sample_name = patient + ".tumor",
            bam_files = patient_tumor_bam_files,
            threads = merge_bam_threads
        }
      }

      if (skip_align && length(patient_normal_bam_files) > 1) {
        call MergeBams as MergeNormalSkipAlignBams {
          input:
            sample_name = patient + ".normal",
            bam_files = patient_normal_bam_files,
            threads = merge_bam_threads
        }
      }

      # If skip align and one bam, index the bam
      if (skip_align && length(patient_tumor_bam_files) == 1) {
        call IndexBam as indexTumorBam {
          input:
            bam = patient_tumor_bam_files[0],
            threads = samtools_threads
        }
      }

      if (skip_align && length(patient_normal_bam_files) == 1) {
        call IndexBam as indexNormalBam {
          input:
            bam = patient_normal_bam_files[0],
            threads = samtools_threads
        }
      }

      # Array[File?] all_final_tumor_bams = [indexTumorBam.out_bam, MergeTumorSkipAlignBams.merged_aligned_bam, MergeTumorAlignBams.merged_aligned_bam, TumorAlign.aligned_bam[0]]
      # Array[File?] all_final_tumor_bam_indexes = [indexTumorBam.out_bam_index, MergeTumorSkipAlignBams.merged_aligned_bam_index, MergeTumorAlignBams.merged_aligned_bam_index, TumorAlign.aligned_bam_index[0]]
      # Array[File?] all_final_normal_bams = [indexNormalBam.out_bam, MergeNormalSkipAlignBams.merged_aligned_bam, MergeNormalAlignBams.merged_aligned_bam, NormalAlign.aligned_bam[0]]
      # Array[File?] all_final_normal_bam_indexes = [indexNormalBam.out_bam_index, MergeNormalSkipAlignBams.merged_aligned_bam_index, MergeNormalAlignBams.merged_aligned_bam_index, NormalAlign.aligned_bam_index[0]]

      Array[File?] all_final_tumor_bams = if (skip_align) then [indexTumorBam.out_bam, MergeTumorSkipAlignBams.merged_aligned_bam] else [MergeTumorAlignBams.merged_aligned_bam, TumorAlign.aligned_bam[0]]
      Array[File?] all_final_tumor_bam_indexes = if (skip_align) then [indexTumorBam.out_bam_index, MergeTumorSkipAlignBams.merged_aligned_bam_index] else [MergeTumorAlignBams.merged_aligned_bam_index, TumorAlign.aligned_bam_index[0]]
      Array[File?] all_final_normal_bams = if (skip_align) then [indexNormalBam.out_bam, MergeNormalSkipAlignBams.merged_aligned_bam] else [MergeNormalAlignBams.merged_aligned_bam, NormalAlign.aligned_bam[0]]
      Array[File?] all_final_normal_bam_indexes = if (skip_align) then [indexNormalBam.out_bam_index, MergeNormalSkipAlignBams.merged_aligned_bam_index] else [MergeNormalAlignBams.merged_aligned_bam_index, NormalAlign.aligned_bam_index[0]]

    # Note that the index 0 below is in the case where there's only one BAM after alignment (no skipalign), so the bams aren't merged
    output {
      File tumor_bam_final = select_first(all_final_tumor_bams)
      File tumor_bam_final_index = select_first(all_final_tumor_bam_indexes)
      File normal_bam_final = select_first(all_final_normal_bams)
      File normal_bam_final_index = select_first(all_final_normal_bam_indexes)
    }
}

# Align using pbmm2. If prefer to use a different version, change the SHA in the runtime section, e.g.:
# Version 1.12 (minimap2 2.17) SHA: d16e5df5bfab75ff2defd2984ec6cb7665473e383045e2a075ea1261ae188861
# Version 1.14.99 (minimap2 2.26) SHA: 19ca8f306b0c1c61aad0bf914c2a291b9ea9b8437e28dbdb5d76f93b81a0dbdf
# Version 1.16.0 (minimap2 2.26) SHA: 4a01f9a3ede68cd018d8d2e79f8773f331f21b35c764ac48b8595709dbd417f7
task Align {
  input {
    File bam_file
    File ref_fasta
    File ref_fasta_index
    String sample_name
    String? additional_args
    Boolean strip_kinetics
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
    --log-level INFO --log-file pbmm2.log \
    ~{if(strip_kinetics) then "--strip" else ""} \
    ~{additional_args}
  >>>

  output {
    File aligned_bam = ofile_name
    File aligned_bam_index = ofile_name_index
    File align_log = "pbmm2.log"
  }

  runtime {
    docker: "quay.io/biocontainers/pbmm2@sha256:19ca8f306b0c1c61aad0bf914c2a291b9ea9b8437e28dbdb5d76f93b81a0dbdf"
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
    File out_bam = bam
    File out_bam_index = basename(bam) + ".bai"
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