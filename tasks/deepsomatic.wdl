version 1.0

workflow run_deepsomatic {
    input {
        File tumor_bam
        File normal_bam
        File tumor_bam_index
        File normal_bam_index
        File ref_fasta
        File ref_fasta_index
        File? contig
        Int threads = 64
        String pname
    }

    call call_DeepSomatic {
        input:
            tumor_bam = tumor_bam,
            normal_bam = normal_bam,
            tumor_bam_index = tumor_bam_index,
            normal_bam_index = normal_bam_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            contig = contig,
            pname = pname,
            threads = threads
    }

    call call_clair3 as clair3_normal {
        input:
            aligned_bam = normal_bam,
            aligned_bam_index = normal_bam_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            contig = contig,
            pname = pname + ".normal",
            threads = threads / 2
    }

    call correct_vcf as correct_clair3_normal {
        input:
            vcf = clair3_normal.clair3_vcf,
            vcf_index = clair3_normal.clair3_vcf_tbi,
            reference = ref_fasta,
            reference_index = ref_fasta_index
    }

    call call_clair3 as clair3_tumor {
        input:
            aligned_bam = tumor_bam,
            aligned_bam_index = tumor_bam_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            contig = contig,
            pname = pname + ".tumor",
            threads = threads / 2
    }

    call correct_vcf as correct_clair3_tumor {
        input:
            vcf = clair3_tumor.clair3_vcf,
            vcf_index = clair3_tumor.clair3_vcf_tbi,
            reference = ref_fasta,
            reference_index = ref_fasta_index
    }

    output {
        File deepsomatic_vcf = call_DeepSomatic.deepsomatic_vcf
        File deepsomatic_vcf_tbi = call_DeepSomatic.deepsomatic_vcf_tbi
        File deepsomatic_gvcf = call_DeepSomatic.deepsomatic_gvcf
        File deepsomatic_gvcf_tbi = call_DeepSomatic.deepsomatic_gvcf_tbi
        File clair3_normal_vcf = correct_clair3_normal.output_vcf
        File clair3_normal_vcf_tbi = correct_clair3_normal.output_vcf_index
        File clair3_tumor_vcf = correct_clair3_tumor.output_vcf
        File clair3_tumor_vcf_tbi = correct_clair3_tumor.output_vcf_index
    }
}

task call_DeepSomatic {
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

    /opt/deepvariant/bin/deepsomatic/run_deepsomatic --version

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

    # Filter for PASS
    bcftools --version
    bcftools view \
        -f PASS -Oz \
        -o deepvariant_~{pname}/~{pname + "_deepsomatic.PASS.vcf.gz"} \
        deepvariant_~{pname}/~{pname + "_deepsomatic.vcf.gz"}
    tabix -p vcf deepvariant_~{pname}/~{pname + "_deepsomatic.PASS.vcf.gz"}
    >>>

    output {
        File deepsomatic_vcf = "deepvariant_" + pname + "/" + pname + "_deepsomatic.PASS.vcf.gz"
        File deepsomatic_vcf_tbi = "deepvariant_" + pname + "/" + pname + "_deepsomatic.PASS.vcf.gz.tbi"
        File deepsomatic_gvcf = "deepvariant_" + pname + "/" + pname + "_deepsomatic.g.vcf.gz"
        File deepsomatic_gvcf_tbi = "deepvariant_" + pname + "/" + pname + "_deepsomatic.g.vcf.gz.tbi"
    }

    runtime {
        docker: "google/deepsomatic:1.6.0"
        cpu: threads
        memory: "~{threads * 6} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}

task call_clair3 {
    input {
        File aligned_bam
        File aligned_bam_index
        File ref_fasta
        File ref_fasta_index
        File? contig 
        String pname
        Int threads
        String clair_model = "hifi_revio"
    }

    Float file_size = ceil(size(aligned_bam, "GB") * 2 + size(ref_fasta, "GB") + size(contig, "GB") + 20)

    command <<<
    set -euxo pipefail

    /opt/bin/run_clair3.sh --version
    mkdir clair3_~{pname}

    /opt/bin/run_clair3.sh \
        --bam_fn=~{aligned_bam} \
        --ref_fn=~{ref_fasta} \
        --threads=~{threads} \
        --platform="hifi" \
        --model_path="/opt/models/~{clair_model}" \
        --output=clair3_~{pname} \
        --sample_name=~{pname}
    
    # Filter for PASS only
    gunzip -c clair3_~{pname}/merge_output.vcf.gz |\
        awk '{if($7!="LowQual"){print $0}}' OFS=$'\t' |\
        bgzip -c > ~{pname}.clair3.small_variants.vcf.gz
    tabix -p vcf ~{pname}.clair3.small_variants.vcf.gz
    >>>

    output {
        File clair3_vcf = pname + ".clair3.small_variants.vcf.gz"
        File clair3_vcf_tbi = pname + ".clair3.small_variants.vcf.gz.tbi"
    }

    runtime {
        docker: "hkubal/clair3@sha256:857af16c759b0893fc757511a17c1efdfe253cbb64dffbcc8eecac0d33a60f60"
        cpu: threads
        memory: "~{threads * 2} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}

# Clair3 output replace all non-ACTG characters with N.
# Fill it back with bcftools
task correct_vcf {
    input {
        File vcf
        File vcf_index
        File reference
        File reference_index
    }

    Float file_size = ceil(size(vcf, "GB") + 10)
    
    command <<<
        set -euxo pipefail

        bcftools --version

        bcftools +fill-from-fasta \
            ~{vcf} \
            -Oz \
            -- -c REF \
            --fasta ~{reference} > ~{basename(vcf, ".vcf.gz") + "correct_ref.vcf.gz"}

        tabix -p vcf ~{basename(vcf, ".vcf.gz") + "correct_ref.vcf.gz"}
    >>>

    output {
        File output_vcf = basename(vcf, ".vcf.gz") + "correct_ref.vcf.gz"
        File output_vcf_index = basename(vcf, ".vcf.gz") + "correct_ref.vcf.gz.tbi"
    }

    runtime {
        docker: "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
        cpu: 2
        memory: "8 GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}