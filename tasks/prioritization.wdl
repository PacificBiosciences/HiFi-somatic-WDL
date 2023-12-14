version 1.0

task prioritize_dmr_intogen {
    input {
        Array[File] dmr_files
        Int ncg_to_filter = 50
        Int threads
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
    }

    runtime {
        docker: "quay.io/pacbio/somatic_general_tools@sha256:6d8c96585ef32007d7fd375984cb6c85dade282e018982254521ee0c685a0166"
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

    csvtk join -t \
        ~{annotSV_tsv} \
        /app/Compendium_Cancer_Genes.tsv \
        -f "Gene_name;SYMBOL" |\
            csvtk filter2 -t -f '$Annotation_mode == "split"' |\
            csvtk summary -t -g "$(csvtk headers -t ~{annotSV_tsv} | tr '\n' ',' | sed 's/,$//g')" \
                -f CANCER_TYPE:collapse,COHORT:collapse,TRANSCRIPT:collapse,MUTATIONS:collapse,ROLE:collapse,CGC_GENE:collapse,CGC_CANCER_GENE:collapse,DOMAINS:collapse,2D_CLUSTERS:collapse,3D_CLUSTERS:collapse -s ";" |\
                sed 's/:collapse//g'  > ~{sub(basename(annotSV_tsv), "\\.tsv$", "")}_intogenCCG.tsv
    >>>

    output {
        File annotSV_intogen_tsv = sub(basename(annotSV_tsv), "\\.tsv$", "") + "_intogenCCG.tsv"
    }

    runtime {
        docker: "quay.io/pacbio/somatic_general_tools@sha256:6d8c96585ef32007d7fd375984cb6c85dade282e018982254521ee0c685a0166"
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
        docker: "quay.io/pacbio/somatic_general_tools@sha256:6d8c96585ef32007d7fd375984cb6c85dade282e018982254521ee0c685a0166"
        cpu: threads
        memory: "~{threads * 4} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}
