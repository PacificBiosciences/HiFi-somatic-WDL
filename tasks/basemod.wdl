version 1.1

task DSS_DMR {
    input {
        File tumor_bed
        File normal_bed
        String pname
        Int threads
    }

    command <<<
        set -euxo pipefail

        cut -f1,2,6,7 ~{tumor_bed} > ~{basename(tumor_bed) + ".tmp"}
        cut -f1,2,6,7 ~{normal_bed} > ~{basename(normal_bed) + ".tmp"}

        Rscript --vanilla /app/DSS_tumor_normal.R \
            ~{basename(tumor_bed) + ".tmp"} \
            ~{basename(normal_bed) + ".tmp"} \
            ~{pname + ".DMR.tsv"} \
            ~{threads}

        rm -f ~{basename(tumor_bed) + ".tmp"} ~{basename(normal_bed) + ".tmp"}
    >>>

    output {
        File output_DMR = pname + ".DMR.tsv"
    }

    runtime {
        docker: "kpinpb/dss:v0.1"
        cpu: threads
        memory: "~{threads * 4} GB"
    }
}

task annotate_DMR {
    input {
        File DMR
        String pname
        Int threads
    }

    command <<<
        set -euxo pipefail

        Rscript --vanilla /app/annotatr_dmr.R \
            ~{DMR} \
            ~{pname} \
            ~{threads}
    >>>

    output {
        Array[File]+ output_DMR_annotated = glob(pname + "*.tsv.gz")
    }

    runtime {
        docker: "kpinpb/annotatr:v0.1"
        cpu: threads
        memory: "~{threads * 4} GB"
    }
}