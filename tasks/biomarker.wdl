version 1.0

workflow owl_msi {
    input {
        File bam_file
        File bam_index
        String sample_name
        File? region_bed
        String? additional_args
        Int minimum_depth = 5
    }

    call owl_profile {
        input:
            bam_file = bam_file,
            bam_index = bam_index,
            sample_name = sample_name,
            region_bed = region_bed,
            additional_args = additional_args
    }

    call owl_score {
        input:
            owl_profile_output = owl_profile.owl_profile_output,
            sample_name = sample_name,
            minimum_depth = minimum_depth
    }

    output {
        File owl_profile_output = owl_profile.owl_profile_output
        File owl_score_output = owl_score.owl_score_output
        File owl_motif_counts = owl_score.owl_motif_counts
    }
}

# Profile samples short repeat regions
task owl_profile {
    input {
        File bam_file
        File bam_index
        String sample_name
        File? region_bed
        String? additional_args
    }

    # If region bed not defined, use default regions
    String region_bed_file = select_first([region_bed, "/opt/owl/data/Simple-repeats-50k.filt.bed.gz"])

    command <<<
    set -euxo pipefail

    echo "Running OWL-MSI profile on ~{bam_file}"

    # If region file ends with .gz, uncompress it
    if [[ "~{region_bed_file}" == *.gz ]]; then
        gunzip -c ~{region_bed_file} > region_to_profile.bed
    else
        cp ~{region_bed_file} region_to_profile.bed
    fi

    owl profile --bam ~{bam_file} \
        --regions region_to_profile.bed \
        --sample ~{sample_name} \
        ~{additional_args} \
        > ~{sample_name}.owl.txt
    >>>

    output {
        File owl_profile_output = "~{sample_name}.owl.txt"
    }

    runtime {
        docker: "kpinpb/owl@sha256:ea9d86643a437774725b017b1a2b432783388888a9b9bd3a2ef28d569b1fd2f7"
        cpu: 2
        memory: "4 GB"
        disk: "2 GB"
        maxRetries: 2
        preemptible: 1
    }
}

# Calculate MSI score from owl profile
task owl_score {
    input {
        File owl_profile_output
        Int minimum_depth = 5
        String sample_name
    }

    command <<<
    set -euxo pipefail

    echo "Calculating MSI score from ~{owl_profile_output}"

    owl score --file ~{owl_profile_output} \
        --prefix ~{sample_name} \
        --min-depth ~{minimum_depth}

    >>>

    output {
        File owl_score_output = "~{sample_name}.owl-scores.txt"
        File owl_motif_counts = "~{sample_name}.owl-motif-counts.txt"
    }

    runtime {
        docker: "kpinpb/owl@sha256:ea9d86643a437774725b017b1a2b432783388888a9b9bd3a2ef28d569b1fd2f7"
        cpu: 2
        memory: "4 GB"
        disk: "2 GB"
        maxRetries: 2
        preemptible: 1
    }
}

task tmb_estimate {
    input {
        File vcf
        File coverage
        String pname
        File? region_bed
        Int min_depth = 10
        Int min_coverage = 10
        Float min_vaf = 0.10
        Float gnomad_max_af_cutoff = 0.01
    }
    
    command <<<
    set -euxo pipefail

    echo "Estimating TMB from ~{vcf} and ~{coverage}"

    python /opt/venv/bin/calculate_tmb.py \
        --vcf ~{vcf} \
        --coverage ~{coverage} \
        --min-depth ~{min_depth} \
        --min-coverage ~{min_coverage} \
        --min-vaf ~{min_vaf} \
        --output ~{pname}.tmb_estimate.json \
        --max-af-cutoff ~{gnomad_max_af_cutoff} \
        --debug-tsv ~{pname}.tmb_estimate.tsv \
        ~{"--region_bed " + region_bed}

    # Repeat but with the built-in gencode regions bed
    python /opt/venv/bin/calculate_tmb.py \
        --vcf ~{vcf} \
        --coverage ~{coverage} \
        --min-depth ~{min_depth} \
        --min-coverage ~{min_coverage} \
        --min-vaf ~{min_vaf} \
        --output ~{pname}.tmb_estimate.gencode_coding.json \
        --max-af-cutoff ~{gnomad_max_af_cutoff} \
        --debug-tsv ~{pname}.tmb_estimate.gencode_coding.tsv \
        --region-bed /opt/gencode_46_coding.bed.gz
    >>>
    
    output {
        File tmb_estimate_output = "~{pname}.tmb_estimate.json"
        File tmb_estimate_output_gencode = "~{pname}.tmb_estimate.gencode_coding.json"
    }

    runtime {
        docker: "kpinpb/tmb-calculator@sha256:2bfd33da4865ccce6b7705aebce13cc7126fe2fe65db95950275e05136a4cc07"
        cpu: 2
        memory: "4 GB"
        disk: "2 GB"
        maxRetries: 2
        preemptible: 1
    }
}