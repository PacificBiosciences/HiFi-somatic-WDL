version 1.1

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

    command <<<
        set -euxo pipefail

        echo "Running hiphase for ~{pname}"

        hiphase --bam ~{bam} \
            -t ~{threads} \
            --output-bam ~{sub(basename(bam), "\\.bam$", ".hiphase.bam")} \
            --vcf ~{vcf} \
            --output-vcf ~{sub(basename(vcf), "\\.vcf.gz", ".hiphase.vcf.gz")} \
            -r ~{ref_fasta} \
            --stats-file ~{sub(basename(bam), "\\.bam$", ".hiphase.stats")} \
            --ignore-read-groups
    >>>

    output {
        File hiphase_bam = sub(basename(bam), "\\.bam$", ".hiphase.bam")
        File hiphase_vcf = sub(basename(vcf), "\\.vcf.gz$", ".hiphase.vcf.gz")
        File hiphase_stats = sub(basename(bam), "\\.bam$", ".hiphase.stats")
    }

    runtime {
        container: "quay.io/pacbio/hiphase:0.10.2"
        cpu: threads
        memory: "~{threads * 6} GB"
    }
}