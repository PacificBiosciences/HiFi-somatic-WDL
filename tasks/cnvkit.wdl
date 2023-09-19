version 1.1

# Call SNV with CNVKit using 10 kbp segmentation window and CBS algorithm
task cnvkit_tumor {
  input {
    File tumor_bam
    File tumor_bam_index
    File normal_bam
    File normal_bam_index
    File ref_fasta
    File ref_fasta_index
    File refFlat
    String pname
    Int threads
  }

  command <<<
    set -euxo pipefail
    
    echo "Running cnvkit tumor for ~{pname}"
    cnvkit.py batch \
        ~{tumor_bam} \
        --normal ~{normal_bam} \
        --annotate ~{refFlat} \
        -f ~{ref_fasta} \
        --target-avg-size 10000 -m wgs \
        -p ~{threads} \
        --diagram --scatter \
        --output-dir ~{pname + "_cnvkit"} \
        --segment-method cbs
  >>>

  output {
    Array[File]+ cnvkit_output = glob(pname + "_cnvkit/*")
  }

  runtime {
    container: "quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0"
    cpu: threads
    memory: "~{threads * 6} GB"
  }
}