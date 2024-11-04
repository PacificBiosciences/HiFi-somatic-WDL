version 1.0

task hiphase {
    input {
        File bam
        File bam_index
        File vcf
        File vcf_index
        String pname
        File ref_fasta
        File ref_fasta_index
        Int threads
    }

    Float file_size = ceil(size(bam, "GB") + size(vcf, "GB") + size(ref_fasta, "GB") + 20)

    command <<<
    set -euxo pipefail

    echo "Running hiphase for ~{pname}"

    hiphase --version
    
    hiphase --bam ~{bam} \
        -t ~{threads} \
        --output-bam ~{sub(basename(bam), "\\.bam$", ".hiphase.bam")} \
        --vcf ~{vcf} \
        --output-vcf ~{sub(basename(vcf), "\\.vcf.gz", ".hiphase.vcf.gz")} \
        -r ~{ref_fasta} \
        --stats-file ~{sub(basename(bam), "\\.bam$", ".hiphase.stats")} \
        --summary-file ~{sub(basename(bam), "\\.bam$", ".hiphase.summary.tsv")} \
        --ignore-read-groups
    >>>

    output {
        File hiphase_bam = sub(basename(bam), "\\.bam$", ".hiphase.bam")
        File hiphase_bam_index = sub(basename(bam), "\\.bam$", ".hiphase.bam.bai")
        File hiphase_vcf = sub(basename(vcf), "\\.vcf.gz$", ".hiphase.vcf.gz")
        File hiphase_stats = sub(basename(bam), "\\.bam$", ".hiphase.stats")
        File hiphase_summary = sub(basename(bam), "\\.bam$", ".hiphase.summary.tsv")
        Array[File] hiphase_stats_summary = [hiphase_stats, hiphase_summary]
    }

    runtime {
        docker: "quay.io/pacbio/hiphase@sha256:c46c8493be8b308c0433441cbafcc1b6ac999dfa6e85001d466ebd551c4a8cf0"
        cpu: threads
        memory: "~{threads * 6} GB"
        disk: file_size + "GB"
        maxRetries: 2
        preemptible: 1
    }
}

task hiphase_with_somatic {
    input {
        File bam
        File bam_index
        File vcf
        File vcf_index
        File somatic_SNP_indel_vcf
        File somatic_SNP_indel_vcf_index
        String pname
        File ref_fasta
        File ref_fasta_index
        Int threads
    }

    Float file_size = ceil(size(bam, "GB") + size(vcf, "GB") + size(ref_fasta, "GB") + 20)

    command <<<
    set -euxo pipefail

    echo "Running hiphase for ~{pname}"

    hiphase --version
    
    hiphase --bam ~{bam} \
        -t ~{threads} \
        --output-bam ~{sub(basename(bam), "\\.bam$", ".hiphase.bam")} \
        --vcf ~{vcf} \
        --output-vcf ~{sub(basename(vcf), "\\.vcf.gz", ".hiphase.vcf.gz")} \
        --vcf ~{somatic_SNP_indel_vcf} \
        --output-vcf ~{sub(basename(somatic_SNP_indel_vcf), "\\.vcf.gz", ".hiphase.vcf.gz")} \
        -r ~{ref_fasta} \
        --stats-file ~{sub(basename(bam), "\\.bam$", ".hiphase.stats")} \
        --summary-file ~{sub(basename(bam), "\\.bam$", ".hiphase.summary.tsv")} \
        --ignore-read-groups
    >>>

    output {
        File hiphase_bam = sub(basename(bam), "\\.bam$", ".hiphase.bam")
        File hiphase_bam_index = sub(basename(bam), "\\.bam$", ".hiphase.bam.bai")
        File hiphase_vcf = sub(basename(vcf), "\\.vcf.gz$", ".hiphase.vcf.gz")
        File hiphase_somatic_small_variants_vcf = sub(basename(somatic_SNP_indel_vcf), "\\.vcf.gz$", ".hiphase.vcf.gz")
        File hiphase_stats = sub(basename(bam), "\\.bam$", ".hiphase.stats")
        File hiphase_summary = sub(basename(bam), "\\.bam$", ".hiphase.summary.tsv")
        Array[File] hiphase_stats_summary = [hiphase_stats, hiphase_summary]
    }

    runtime {
        docker: "quay.io/pacbio/hiphase@sha256:47fe7d42aea6b1b2e6d3c7401bc35a184464c3f647473d0525c00f3c968b40ad"
        cpu: threads
        memory: "~{threads * 6} GB"
        disk: file_size + "GB"
        maxRetries: 2
        preemptible: 1
    }
}

task longphase_with_somatic {
    input {
        File bam
        File bam_index
        File vcf
        File vcf_index
        File somatic_SNP_indel_vcf
        File somatic_SNP_indel_vcf_index
        String pname
        File ref_fasta
        File ref_fasta_index
        Int threads
    }

    Float file_size = ceil(size(bam, "GB") + size(vcf, "GB") + size(ref_fasta, "GB") + 20)

    command <<<
    set -euxo pipefail

    echo "Running Longphase for ~{pname}"

    bcftools concat -a ~{vcf} ~{somatic_SNP_indel_vcf} | bcftools sort -Oz -o merged.vcf.gz

    longphase phase \
        -r ~{ref_fasta} \
        -s merged.vcf.gz \
        -b ~{bam} \
        -t ~{threads} \
        --pb \
        -o tmp.phased

    longphase haplotag \
        -r ~{ref_fasta} \
        -s tmp.phased.vcf \
        -b ~{bam} \
        -t ~{threads} \
        -o ~{sub(basename(bam), "\\.bam$", ".longphase")}
    
    # index bam file
    samtools index -@~{threads} ~{sub(basename(bam), "\\.bam$", ".longphase.bam")}
    
    # Split into somatic and germline VCF
    bcftools view tmp.phased.vcf |\
        grep -v "NAF" > ~{sub(basename(vcf), "\\.vcf.gz$", ".longphase.vcf")}
    bcftools view tmp.phased.vcf |\
        grep -P "^#|NAF" > ~{sub(basename(somatic_SNP_indel_vcf), "\\.vcf.gz$", ".longphase.vcf")}
    
    rm -f tmp.phased.vcf merged.vcf.gz

    bgzip -@4 ~{sub(basename(vcf), "\\.vcf.gz$", ".longphase.vcf")}
    bgzip -@4 ~{sub(basename(somatic_SNP_indel_vcf), "\\.vcf.gz$", ".longphase.vcf")}
    >>>

    output {
        File longphase_bam = sub(basename(bam), "\\.bam$", ".longphase.bam")
        File longphase_bam_index = sub(basename(bam), "\\.bam$", ".longphase.bam.bai")
        File longphase_vcf = sub(basename(vcf), "\\.vcf.gz$", ".longphase.vcf.gz")
        File longphase_somatic_small_variants_vcf = sub(basename(somatic_SNP_indel_vcf), "\\.vcf.gz$", ".longphase.vcf.gz")
    }

    runtime {
        docker: "quay.io/pacbio/longphase@sha256:2d57c35da90fc74f56574c439e1a44f74c7f787f3812de3c402cf9206c9d4be3"
        cpu: threads
        memory: "~{threads * 6} GB"
        disk: file_size + "GB"
        maxRetries: 2
        preemptible: 1
    }
}

