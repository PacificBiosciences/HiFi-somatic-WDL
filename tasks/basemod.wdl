version 1.0

task DSS_DMR {
    input {
        File tumor_bed
        File normal_bed
        String pname
        Int threads
    }

    Float file_size = ceil(size(tumor_bed, "GB") + size(normal_bed, "GB") + 10)

    command <<<
        set -euxo pipefail
        _PROF_START=$(date +%s)

        # Extract methylation info from bed files for DSS
        gunzip -c ~{tumor_bed} | cut -f1,2,6,7 > ~{basename(tumor_bed) + ".tmp"}
        gunzip -c ~{normal_bed} | cut -f1,2,6,7 > ~{basename(normal_bed) + ".tmp"}

        # Run DSS
        Rscript --vanilla /app/DSS_tumor_normal.R \
            ~{basename(tumor_bed) + ".tmp"} \
            ~{basename(normal_bed) + ".tmp"} \
            ~{pname + ".DMR.tsv"} \
            ~{threads}

        rm -f ~{basename(tumor_bed) + ".tmp"} ~{basename(normal_bed) + ".tmp"}

        # --- Resource profiling output ---
        _PROF_END=$(date +%s)
        _PROF_WALL=$(( _PROF_END - _PROF_START ))
        _VM_PEAK=$(grep VmPeak /proc/self/status 2>/dev/null | awk '{print $2}' || echo "unavailable")
        echo -e "metric\tvalue" > resource_usage.tsv
        echo -e "task_name\tDSS_DMR" >> resource_usage.tsv
        echo -e "wall_clock_seconds\t${_PROF_WALL}" >> resource_usage.tsv
        echo -e "vm_peak_kb\t${_VM_PEAK}" >> resource_usage.tsv
        echo -e "profiling_method\tlightweight" >> resource_usage.tsv
    >>>

    output {
        File output_DMR = pname + ".DMR.tsv"
        File resource_usage = "resource_usage.tsv"
    }

    runtime {
        docker: "quay.io/pacbio/somatic_r_tools@sha256:68dc04908a37e26b30dc9795fa6cc0e85a238c8695afe805ad164a071193fb48"
        cpu: threads
        memory: "~{threads * 4} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}

task annotate_DMR {
    input {
        File DMR
        String pname
        Int threads
    }

    Float file_size = ceil(size(DMR, "GB"))

    command <<<
        set -euxo pipefail
        _PROF_START=$(date +%s)

        # Annotate DMRs
        Rscript --vanilla /app/annotatr_dmr.R \
            ~{DMR} \
            ~{pname} \
            ~{threads}

        # --- Resource profiling output ---
        _PROF_END=$(date +%s)
        _PROF_WALL=$(( _PROF_END - _PROF_START ))
        _VM_PEAK=$(grep VmPeak /proc/self/status 2>/dev/null | awk '{print $2}' || echo "unavailable")
        echo -e "metric\tvalue" > resource_usage.tsv
        echo -e "task_name\tannotate_DMR" >> resource_usage.tsv
        echo -e "wall_clock_seconds\t${_PROF_WALL}" >> resource_usage.tsv
        echo -e "vm_peak_kb\t${_VM_PEAK}" >> resource_usage.tsv
        echo -e "profiling_method\tlightweight" >> resource_usage.tsv
    >>>

    output {
        Array[File]+ output_DMR_annotated = glob(pname + "*.tsv.gz")
        File resource_usage = "resource_usage.tsv"
    }

    runtime {
        docker: "quay.io/pacbio/somatic_r_tools@sha256:68dc04908a37e26b30dc9795fa6cc0e85a238c8695afe805ad164a071193fb48"
        cpu: threads
        memory: "~{threads * 4} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}
