version 1.0


# Call VCF and SNF with sniffles
task sniffles {
  input {
    String pname
    File bam
    File bam_index
    File ref_fasta
    File trf_bed
    Int threads
  }

  Float file_size = ceil(size(bam, "GB") * 2 + size(ref_fasta, "GB") + size(trf_bed, "GB") + 10)

  command <<<
    set -euxo pipefail

    echo "Running sniffles2 with noqc flag for ~{bam}"
    sniffles -t ~{threads} \
      --no-qc --non-germline \
      --input ~{bam} \
      --reference ~{ref_fasta} \
      --tandem-repeats ~{trf_bed} \
      --vcf ~{pname}.sniffles.vcf \
      --snf ~{pname}.sniffles.snf \
      --sample-id ~{pname}
  >>>

  output {
    File output_vcf = pname + ".sniffles.vcf"
    File output_snf = pname + ".sniffles.snf"
  }

  runtime {
    docker: "quay.io/biocontainers/sniffles:2.0.7--pyhdfd78af_0"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

# SNFs of tumor normal pair are fed to this step to do joint-calling
task sniffles_call_snf {
  input {
    String pname
    File normal_snf
    File tumor_snf
    File ref_fasta
    File trf_bed
    Int threads
  }

  Float file_size = ceil(size(normal_snf, "GB") + size(tumor_snf, "GB") + size(ref_fasta, "GB") + size(trf_bed, "GB") + 10)

  command <<<
    set -euxo pipefail

    echo "Running sniffles2 joint call for patient ~{pname}"
    sniffles -t ~{threads} \
      --no-qc --non-germline \
      --combine-output-filtered \
      --combine-low-confidence-abs 1 \
      --combine-low-confidence 0 \
      --input ~{normal_snf} ~{tumor_snf} \
      --reference ~{ref_fasta} \
      --tandem-repeats ~{trf_bed} \
      --vcf ~{pname}.sniffles.vcf 
  >>>

  output {
    File output_vcf = pname + ".sniffles.vcf"
  }

  runtime {
    docker: "quay.io/biocontainers/sniffles:2.0.7--pyhdfd78af_0"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

# Use slivar to filter for somatic variants
# Require >1 reads in tumor and 0 reads in normal and
# >0 reads in normal
task slivar_select_somatic {
  input {
    File sniffles_join_vcf
    String normal_name
    String tumor_name
    Int threads
  }

  Float file_size = ceil(size(sniffles_join_vcf, "GB") + 10)

  command <<<
    set -euxo pipefail
    
    echo -e "#tumor\tnormal\n~{tumor_name}\t~{normal_name}" > tumor_normal.alias
    echo "Running slivar select somatic for ~{sniffles_join_vcf}"
    slivar expr --vcf ~{sniffles_join_vcf} \
      --alias tumor_normal.alias \
      --group-expr "tumor_only:tumor.DV>1 && normal.DV==0 && normal.DR + normal.DV > 0" >\
      ~{sub(basename(sniffles_join_vcf), "\\.vcf", ".slivar.vcf")}
  >>>

  output {
    File output_vcf = sub(basename(sniffles_join_vcf), "\\.vcf", ".slivar.vcf")
  }

  runtime {
    docker: "quay.io/biocontainers/slivar:0.3.0--h4e814b3_2"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

# Use bcftools to filter for somatic variants based on slivar's annotation
task bcftools_filter {
  input {
    File slivar_join_vcf
    File ref_bed
    String tumor_name
    Int threads
  }

  Float file_size = ceil(size(slivar_join_vcf, "GB") + size(ref_bed, "GB") + 10)

  command <<<
    set -euxo pipefail
    
    echo "Running bcftools filter for ~{slivar_join_vcf}"

    bgzip ~{slivar_join_vcf}
    tabix ~{slivar_join_vcf + ".gz"}

    # Sed expression below is because Sniffles uses 0-based coord for chrM
    # for some reasons.
    bcftools filter -i 'tumor_only=="~{tumor_name}"' ~{slivar_join_vcf + ".gz"} \
      -R ~{ref_bed} |\
      sed 's/chrM:0/chrM:1/g' |\
      bcftools filter -i 'SVLEN>=50 || SVLEN<=-50 || SVTYPE=="BND"' \
        -Oz \
        -o ~{sub(basename(slivar_join_vcf), "\\.vcf", ".somatic.vcf.gz")}

    tabix ~{sub(basename(slivar_join_vcf), "\\.vcf", ".somatic.vcf.gz")}
  >>>

  output {
    File output_vcf = sub(basename(slivar_join_vcf), "\\.vcf", ".somatic.vcf.gz")
    File output_vcf_index = sub(basename(slivar_join_vcf), "\\.vcf", ".somatic.vcf.gz.tbi")
  }

  runtime {
    docker: "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

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
    docker: "quay.io/biocontainers/severus:0.1.1--pyhdfd78af_0"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}