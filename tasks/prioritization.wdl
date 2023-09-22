version 1.1

task prioritize_dmr_intogen {
    input {
        Array[File] dmr_files
        Int threads
    }

    command <<<
    set -euxo pipefail

    # Join with intogen CCG genes and filter for at 
    # least 100 CGs
    for dmr_file in ~{sep(" ", dmr_files)}
    do
    # Annotate only if the file is not *_summary.tsv.gz
        if [[ ! ${dmr_file} == *"_summary.tsv.gz" ]]
        then
            csvtk join -t \
                ${dmr_file} \
                /app/Compendium_Cancer_Genes.tsv \
                -f "annot.symbol;SYMBOL" |\
            csvtk filter -t -f 'nCG>100' >\
                $(basename ${dmr_file%%.tsv.gz})_intogen.nCG100.tsv
        fi
    done

    # Summarize and collapse genes coordinate + metadata for each DMR entries
    # As each DMR can overlap multiple promoters, the summary function here
    # collapse the metadata for all promoters for each DMR.
    for intogen_file in *_intogen.nCG100.tsv
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
    rm -f ./*_intogen.nCG100.tsv
    >>>

    output {
        Array[File]+ DMR_nCG100_CCG = glob("*_intogen.nCG100_summary.tsv.gz")
    }

    runtime {
        docker: "kpinpb/general_tools:v0.1"
        cpu: threads
        memory: "~{threads * 4} GB"
    }
}