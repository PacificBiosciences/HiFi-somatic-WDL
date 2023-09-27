version 1.0

task vep_annotate {
    input {
        File input_vcf
        File? vep_cache
        File ref_fasta
        File ref_fasta_index
        String pname
        Int threads
    }

    Float file_size = ceil(size(input_vcf, "GB") + size(vep_cache, "GB") + size(ref_fasta, "GB") + size(ref_fasta_index, "GB"))

    command <<<
        set -euxo pipefail

        mkdir -p vep_data/
        # If vep_cache is not provided, fail
        if [ ! -f ~{vep_cache} ]; then
            echo "VEP cache file not found. Please provide a valid cache file."
            exit 1
        fi

        tar -xzvf ~{vep_cache} -C vep_data/
        vep \
            --cache \
            --offline \
            --dir vep_data/ \
            --fasta ~{ref_fasta} \
            --format vcf \
            --fork ~{threads} \
            --species homo_sapiens \
            --assembly GRCh38 \
            --symbol \
            --hgvs \
            --refseq \
            --check_existing \
            --vcf \
            --pick \
            --flag_pick_allele_gene \
            --everything \
            --compress_output bgzip \
            -i ~{input_vcf} \
            -o ~{sub(basename(input_vcf), "\\.vcf.gz$", "")}.vep.vcf.gz

        # Delete cache after annotation
        rm -rf vep_data/
    >>>

    output {
        File vep_annotated_vcf = sub(basename(input_vcf), "\\.vcf.gz$", "") + ".vep.vcf.gz"
    }

    runtime {
        docker: "ensemblorg/ensembl-vep:release_110.1"
        cpu: threads
        memory: "~{threads * 4} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}

task annotsv {
    input {
        File sv_vcf
        File sv_vcf_index
        File? annotsv_cache
        String pname
        Int threads
    }

    Float file_size = ceil(size(sv_vcf, "GB") + size(annotsv_cache, "GB"))

    command <<<
        set -euxo pipefail

        mkdir -p annotsv_cache_dir
        # If annotsv_cache is not provided, fail
        if [ ! -f ~{annotsv_cache} ]; then
            echo "AnnotSV cache file not found. Please provide a valid cache file."
            exit 1
        fi

        tar -xzvf ~{annotsv_cache} -C annotsv_cache_dir/

        AnnotSV \
            -annotationsDir annotsv_cache_dir/AnnotSV/ \
            -SVinputFile ~{sv_vcf} \
            -outputDir . \
            -outputFile ~{sub(basename(sv_vcf), "\\.vcf.gz$", "")}.annotsv.tsv \
            -SVinputInfo 1 \
            -genomeBuild GRCh38

        # Delete cache after annotation
        rm -rf annotsv_cache_dir/
    >>>

    output {
        File annotsv_annotated_tsv = sub(basename(sv_vcf), "\\.vcf.gz$", "") + ".annotsv.tsv"
    }

    runtime {
        docker: "quay.io/biocontainers/annotsv:3.3.6--py311hdfd78af_0"
        cpu: threads
        memory: "~{threads * 4} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}