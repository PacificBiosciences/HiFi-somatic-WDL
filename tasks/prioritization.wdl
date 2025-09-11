version 1.0

task prioritize_dmr_intogen {
    input {
        Array[File] dmr_files
        Int ncg_to_filter = 50
        Int threads
        String pname
    }

    Float file_size = ceil(size(dmr_files[0], "GB") * length(dmr_files) + 10)

    command <<<
    set -euxo pipefail

    csvtk version

    # Join with intogen CCG genes and filter for at 
    # least 100 CGs
    for dmr_file in ~{sep=" " dmr_files}
    do
    # Annotate only if the file is not *_summary.tsv.gz
        if [[ ! ${dmr_file} == *"_summary.tsv.gz" ]]
        then
            csvtk join -t \
                ${dmr_file} \
                /app/Compendium_Cancer_Genes.tsv \
                -f "annot.symbol;SYMBOL" |\
            csvtk filter -t -f 'nCG>~{ncg_to_filter}' >\
                $(basename ${dmr_file%%.tsv.gz})_intogen.nCG~{ncg_to_filter}.tsv
        fi
    done

    # Summarize and collapse genes coordinate + metadata for each DMR entries
    # As each DMR can overlap multiple promoters, the summary function here
    # collapse the metadata for all promoters for each DMR.
    for intogen_file in *_intogen.nCG~{ncg_to_filter}.tsv
    do
        csvtk mutate2 \
            -t -s -n 'promoter_coord' \
            -e '${annot.seqnames} + ":" + ${annot.start} + "-" + ${annot.end}' ${intogen_file} |\
            csvtk summary -t \
            -g 'seqnames,start,end,width,strand,length,nCG,meanMethy1,meanMethy2,diff.Methy,areaStat,annot.gene_id,annot.symbol,annot.type' \
            -f promoter_coord:collapse,CANCER_TYPE:collapse,COHORT:collapse,TRANSCRIPT:collapse,MUTATIONS:collapse,ROLE:collapse,CGC_GENE:collapse,CGC_CANCER_GENE:collapse,DOMAINS:collapse,2D_CLUSTERS:collapse,3D_CLUSTERS:collapse,annot.tx_id:collapse,annot.id:collapse \
            -s ";" |\
            sed 's/:collapse//g' > ${intogen_file%%.tsv}_summary.tsv
    done

    # Compress 
    for summary_file in *_summary.tsv
    do
        gzip "${summary_file}"
    done

    # Clean up
    rm -f ./*_intogen.nCG~{ncg_to_filter}.tsv
    >>>

    output {
        Array[File]+ DMR_nCGFilter_CCG = glob("*_intogen.nCG*_summary.tsv.gz")
        File promoter_file = pname + "_hg38_genes_promoters_dmrs_intogen.nCG50_summary.tsv.gz"
    }

    runtime {
        docker: "quay.io/pacbio/somatic_general_tools@sha256:99159e099d9044c52debabdc9491b168487aaa37534c1a748809bc69f169812a"
        cpu: threads
        memory: "~{threads * 4} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}

task prioritize_sv_intogen {
    input {
        File annotSV_tsv
        Int threads
    }

    Float file_size = ceil(size(annotSV_tsv, "GB") + 10)

    command <<<
    set -euxo pipefail
    
    csvtk version

    # Remove any quote from the file
    sed 's/"//g' ~{annotSV_tsv} > ~{basename(annotSV_tsv)}_noquote.tsv

    # First do an inner join with the Compendium_Cancer_Genes.tsv
    csvtk join -t \
        ~{basename(annotSV_tsv)}_noquote.tsv \
        /app/Compendium_Cancer_Genes.tsv \
        -f "Gene_name;SYMBOL" |\
            csvtk filter2 -t -f '$Annotation_mode == "split"' |\
            csvtk summary -t -g "$(csvtk headers -t ~{annotSV_tsv} | tr '\n' ',' | sed 's/,$//g')" \
                -f CANCER_TYPE:collapse,COHORT:collapse,TRANSCRIPT:collapse,MUTATIONS:collapse,ROLE:collapse,CGC_GENE:collapse,CGC_CANCER_GENE:collapse,DOMAINS:collapse,2D_CLUSTERS:collapse,3D_CLUSTERS:collapse -s ";" |\
                sed 's/:collapse//g'  > ~{sub(basename(annotSV_tsv), "\\.tsv$", "")}_intogenCCG.tsv
    
    # Then, look for mate pair IDs for all BND events
    { \
        csvtk cut -t -f INFO ~{sub(basename(annotSV_tsv), "\\.tsv$", "")}_intogenCCG.tsv | tail -n+2 | awk -F';' '{for (i=1;i<=NF;i++) if ($i ~ /^MATE_ID=/) {split($i,a,"="); print a[2]}}'; \
    } | csvtk uniq -t | sort > BND_ids.txt

    # Compare the IDs in the intogen TSV with the IDs in the BND_ids.txt
    {
        echo -e "ID"; \
        csvtk cut -t -f "ID" ~{sub(basename(annotSV_tsv), "\\.tsv$", "")}_intogenCCG.tsv | sort | uniq | comm -13 - BND_ids.txt; \
    } > BND_ids_to_recover.txt

    # Filter original rows to the retained IDs so both mates are kept
    csvtk join -t \
        ~{basename(annotSV_tsv)}_noquote.tsv \
        BND_ids_to_recover.txt \
        -f "ID;ID" |\
        csvtk cut -t -f "$(csvtk headers -t ~{annotSV_tsv} | tr '\n' ',' | sed 's/,$//g')" \
        > recovered_BNDs.tsv

    # Concatenate the original TSV with the recovered BNDs
    csvtk concat -t \
        ~{sub(basename(annotSV_tsv), "\\.tsv$", "")}_intogenCCG.tsv \
        recovered_BNDs.tsv \
        > recovered_BNDs_intogenCCG.tsv

    mv recovered_BNDs_intogenCCG.tsv ~{sub(basename(annotSV_tsv), "\\.tsv$", "")}_intogenCCG.tsv

    rm -f ~{basename(annotSV_tsv)}_noquote.tsv ~{basename(annotSV_tsv)}_ccg_id_info.tsv ~{basename(annotSV_tsv)}_retain_ids.tsv ~{basename(annotSV_tsv)}_noquote_retained.tsv BND_ids.txt BND_ids_to_recover.txt recovered_BNDs.tsv
    >>>

    output {
        File annotSV_intogen_tsv = sub(basename(annotSV_tsv), "\\.tsv$", "") + "_intogenCCG.tsv"
    }

    runtime {
        docker: "quay.io/pacbio/somatic_general_tools@sha256:99159e099d9044c52debabdc9491b168487aaa37534c1a748809bc69f169812a"
        cpu: threads
        memory: "~{threads * 4} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}

task prioritize_small_variants {
    input {
        File vep_annotated_vcf
        Int threads
        String pname
    }

    Float file_size = ceil(size(vep_annotated_vcf, "GB") + 10)
    String fname = sub(basename(vep_annotated_vcf), "\\.vcf.gz", "") + ".tsv"
    String fname2 = sub(basename(vep_annotated_vcf), "\\.vcf.gz", "") + "_intogenCCG.tsv"

    command <<<
    set -euxo pipefail

    csvtk version

    echo -e "CHROM\tPOS\tREF\tALT\tFORMAT\t~{pname}\t$(bcftools +split-vep ~{vep_annotated_vcf} -l | cut -f2 | tr '\n' '\t' | sed 's/\t$//g')" > ~{fname}
    bcftools +split-vep ~{vep_annotated_vcf} -A tab -f '%CHROM\t%POS\t%REF\t%ALT\t%FORMAT\t%CSQ\n' >> ~{fname}

    csvtk join -t \
        ~{fname} \
        <(sed 's/DOMAINS/CCG_DOMAINS/g' /app/Compendium_Cancer_Genes.tsv) \
        -f SYMBOL |\
            csvtk summary -t -g "$(csvtk headers -t ~{fname} | tr '\n' ',' | sed 's/,$//g')" \
                -f CANCER_TYPE:collapse,COHORT:collapse,TRANSCRIPT:collapse,MUTATIONS:collapse,ROLE:collapse,CGC_GENE:collapse,CGC_CANCER_GENE:collapse,CCG_DOMAINS:collapse,2D_CLUSTERS:collapse,3D_CLUSTERS:collapse -s ";" |\
                sed 's/:collapse//g' > ~{fname2}
    >>>

    output {
        File vep_annotated_tsv = fname
        File vep_annotated_tsv_intogenCCG = fname2
    }

    runtime {
        docker: "quay.io/pacbio/somatic_general_tools@sha256:99159e099d9044c52debabdc9491b168487aaa37534c1a748809bc69f169812a"
        cpu: threads
        memory: "~{threads * 4} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}

task report_sample {
    input {
        File annotated_small_variant_tsv
        File intogen_small_var_tsv
        File sv_intogen_tsv
        File sv_vcf
        File cnv 
        File purity_ploidy
        File mosdepth_tumor
        File mosdepth_normal
        File mutsig_tsv
        File mut_reconstructed_sigs
        File dmr_intogen_tsv
        File chord_file
        File circos_png
        File tmb_estimate_file
        File tmb_estimate_gencode_file
        String pname
        File vis_file
        String cnv_tool
        File owl_score
        File sv_pairs_file
    }

    Float file_size = ceil(size(annotated_small_variant_tsv, "GB") + size(intogen_small_var_tsv, "GB") + size(sv_intogen_tsv, "GB") + size(sv_vcf, "GB") + size(cnv, "GB") + size(purity_ploidy, "GB") + size(mosdepth_tumor, "GB") + size(mosdepth_normal, "GB") + size(mutsig_tsv, "GB") + size(mut_reconstructed_sigs, "GB") + size(dmr_intogen_tsv, "GB") + size(tmb_estimate_file, "GB") + size(tmb_estimate_gencode_file, "GB") + 10)

    command <<<
    set -euxo pipefail

    # Copy vis_file to the current directory
    cp ~{vis_file} visualize_hifisomatic.qmd

    cp ~{circos_png} circos.png
    
    Rscript -e \
    'quarto::quarto_render(
        "visualize_hifisomatic.qmd",
        output_format = "dashboard",
        execute_params = list(
            vcf_file = "~{annotated_small_variant_tsv}",
            intogen_smallvar = "~{intogen_small_var_tsv}",
            sv_file = "~{sv_intogen_tsv}",
            sv_vcf_file = "~{sv_vcf}",
            cnv_file = "~{cnv}",
            purity_ploidy_file = "~{purity_ploidy}",
            mosdepth_tumor_file = "~{mosdepth_tumor}",
            mosdepth_normal_file = "~{mosdepth_normal}",
            mut_sig_file = "~{mutsig_tsv}",
            mut_context_file = "~{mut_reconstructed_sigs}",
            dmr_promoter_file = "~{dmr_intogen_tsv}",
            chord_file = "~{chord_file}",
            sample_name = "~{pname}",
            cnv_tool = "~{cnv_tool}",
            owl_score = "~{owl_score}",
            tmb_estimate_file = "~{tmb_estimate_file}",
            tmb_estimate_gencode_file = "~{tmb_estimate_gencode_file}",
            sv_pairs_file = "~{sv_pairs_file}"
        )
    )'

    mv visualize_hifisomatic.html ~{pname}_summary_report.html
    >>>

    output {
        File summary_report = pname + "_summary_report.html"
    }

    runtime {
        docker: "quay.io/pacbio/somatic_r_tools@sha256:68dc04908a37e26b30dc9795fa6cc0e85a238c8695afe805ad164a071193fb48"
        cpu: 8
        memory: "32 GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    } 
}

task report_sample_TO {
    input {
        File annotated_small_variant_tsv
        File intogen_small_var_tsv
        File sv_intogen_tsv
        File sv_vcf
        File mosdepth_tumor
        File purity_ploidy
        File cnv_hap1
        File cnv_hap2
        File mutsig_tsv
        File mut_reconstructed_sigs
        File circos_png
        File tmb_estimate_file
        File tmb_estimate_gencode_file
        String pname
        File vis_file
        File owl_score
        File sv_pairs_file
    }

    Float file_size = ceil(size(annotated_small_variant_tsv, "GB") + size(intogen_small_var_tsv, "GB") + size(sv_intogen_tsv, "GB") + size(sv_vcf, "GB") + size(mosdepth_tumor, "GB") + size(cnv_hap1, "GB") + size(cnv_hap2, "GB") + size(purity_ploidy, "GB") + size(mutsig_tsv, "GB") + size(mut_reconstructed_sigs, "GB") + size(tmb_estimate_file, "GB") + size(tmb_estimate_gencode_file, "GB") + 10)

    command <<<
    set -euxo pipefail

    # Copy vis_file to the current directory
    cp ~{vis_file} visualize_hifisomatic.qmd

    cp ~{circos_png} circos.png
    
    Rscript -e \
    'quarto::quarto_render(
        "visualize_hifisomatic.qmd",
        output_format = "dashboard",
        execute_params = list(
            vcf_file = "~{annotated_small_variant_tsv}",
            intogen_smallvar = "~{intogen_small_var_tsv}",
            sv_file = "~{sv_intogen_tsv}",
            sv_vcf_file = "~{sv_vcf}",
            cnv_hap1_file = "~{cnv_hap1}",
            cnv_hap2_file = "~{cnv_hap2}",
            purity_ploidy_file = "~{purity_ploidy}",
             mosdepth_tumor_file = "~{mosdepth_tumor}",
            mut_sig_file = "~{mutsig_tsv}",
            mut_context_file = "~{mut_reconstructed_sigs}",
            sample_name = "~{pname}",
            owl_score = "~{owl_score}",
            tmb_estimate_file = "~{tmb_estimate_file}",
            tmb_estimate_gencode_file = "~{tmb_estimate_gencode_file}",
            sv_pairs_file = "~{sv_pairs_file}"
        )
    )'

    mv visualize_hifisomatic.html ~{pname}_summary_report.html
    >>>

    output {
        File summary_report = pname + "_summary_report.html"
    }

    runtime {
        docker: "quay.io/pacbio/somatic_r_tools@sha256:68dc04908a37e26b30dc9795fa6cc0e85a238c8695afe805ad164a071193fb48"
        cpu: 8
        memory: "32 GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    } 
}
