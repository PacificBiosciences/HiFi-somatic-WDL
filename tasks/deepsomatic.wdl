version 1.0

task DeepSomatic {
    input {
        File tumor_bam
        File normal_bam
        File tumor_bam_index
        File normal_bam_index
        File ref_fasta
        File ref_fasta_index
        File? contig 
        String pname
        Int threads
    }

    Float file_size = ceil(size(tumor_bam, "GB") * 2 + size(normal_bam, "GB") * 2 + size(ref_fasta, "GB") + size(contig, "GB") + 20)

    command <<<
    set -euxo pipefail

    mkdir deepvariant_~{pname}

    /opt/deepvariant/bin/deepsomatic/run_deepsomatic \
        --model_type=PACBIO \
        --ref=~{ref_fasta} \
        --reads_normal=~{normal_bam} \
        --reads_tumor=~{tumor_bam} \
        --output_vcf=deepvariant_~{pname}/~{pname + "_deepsomatic.vcf.gz"} \
        --output_gvcf=deepvariant_~{pname}/~{pname + "_deepsomatic.g.vcf.gz"} \
        --sample_name_tumor=~{pname + ".tumor"} \
        --sample_name_normal=~{pname + ".normal"} \
        --num_shards=~{threads} \
        --logging_dir=deepvariant_~{pname}/logs \
        ~{"--regions=" + contig}
    >>>

    output {
        File deepsomatic_vcf = "deepvariant_" + pname + "/" + pname + "_deepsomatic.vcf.gz"
        File deepsomatic_vcf_tbi = "deepvariant_" + pname + "/" + pname + "_deepsomatic.vcf.gz.tbi"
        File deepsomatic_gvcf = "deepvariant_" + pname + "/" + pname + "_deepsomatic.g.vcf.gz"
        File deepsomatic_gvcf_tbi = "deepvariant_" + pname + "/" + pname + "_deepsomatic.g.vcf.gz.tbi"
    }

    runtime {
        container: "google/deepsomatic:1.6.0"
        cpu: threads
        memory: "~{threads * 4} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}