version 1.0

task Amber {
    input {
        String referenceName
        File referenceBam
        File referenceBamIndex
        String tumorName
        File tumorBam
        File tumorBamIndex
        String outputDir = "amber"
        File ensembl_data_dir_tarball
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        String ref_genome_version = "V38"

        Int threads = 2
    }

    String memory = threads * 4 + "GB"
    String javaXmx = threads * 4 - 4 + "G"
    Float file_size = ceil(size(referenceBam, "GB") + size(tumorBam, "GB") + size(referenceFasta, "GB") + 20)

    command {
        set -euxo pipefail

        mkdir ensembl_data_dir 

        tar -xzf ~{ensembl_data_dir_tarball} -C ensembl_data_dir
        
        java -Xmx~{javaXmx} -jar /app/amber.gamma1000.jar \
            -reference ~{referenceName} \
            -reference_bam ~{referenceBam} \
            -tumor ~{tumorName} \
            -tumor_bam ~{tumorBam} \
            -output_dir ~{outputDir} \
            -threads ~{threads} \
            -ref_genome ~{referenceFasta} \
            -ref_genome_version ~{ref_genome_version} \
            -loci ensembl_data_dir/hmf*/dna/copy_number/AmberGermlineSites.*.tsv.gz

        rm -rf ensembl_data_dir
    }

    output {
        File version = "~{outputDir}/amber.version"
        File tumorBafPcf = "~{outputDir}/~{tumorName}.amber.baf.pcf"
        File tumorBafTsv = "~{outputDir}/~{tumorName}.amber.baf.tsv.gz"
        File tumorContaminationVcf = "~{outputDir}/~{tumorName}.amber.contamination.vcf.gz"
        File tumorContaminationVcfIndex = "~{outputDir}/~{tumorName}.amber.contamination.vcf.gz.tbi"
        File tumorContaminationTsv = "~{outputDir}/~{tumorName}.amber.contamination.tsv"
        File tumorQc = "~{outputDir}/~{tumorName}.amber.qc"
        Array[File] outputs = [version, tumorBafPcf, tumorBafTsv,
            tumorContaminationVcf, tumorContaminationVcfIndex, tumorContaminationTsv, tumorQc]
    }

    runtime {
        docker: "quay.io/pacbio/purple@sha256:9074d56ad46f3d6804f1ae6e45a2ac8c300effc5efe120020d3dcd23edb3b063"
        cpu: threads
        memory: memory
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}

task Cobalt {
    input {
        String referenceName
        File referenceBam
        File referenceBamIndex
        String tumorName
        File tumorBam
        File tumorBamIndex
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        String outputDir = "cobalt"
        File ensembl_data_dir_tarball
        Int threads
        Int pcf_gamma = 1000
    }

    String memory = threads * 4 + "GB"
    String javaXmx = threads * 4 - 4 + "G"
    Float file_size = ceil(size(referenceBam, "GB") + size(tumorBam, "GB") + 20)

    command {
        set -euxo pipefail

        mkdir ensembl_data_dir 

        tar -xzf ~{ensembl_data_dir_tarball} -C ensembl_data_dir

        java -Xmx~{javaXmx} -jar /app/cobalt.jar \
            -reference ~{referenceName} \
            -reference_bam ~{referenceBam} \
            -tumor ~{tumorName} \
            -tumor_bam ~{tumorBam} \
            -ref_genome ~{referenceFasta} \
            -output_dir ~{outputDir} \
            -threads ~{threads} \
            -pcf_gamma ~{pcf_gamma} \
            -validation_stringency SILENT \
            -gc_profile ensembl_data_dir/hmf*/dna/copy_number/GC_profile.*.cnp
        
        rm -rf ensembl_data_dir
    }

    output {
        File version = "~{outputDir}/cobalt.version"
        File normalGcMedianTsv = "~{outputDir}/~{referenceName}.cobalt.gc.median.tsv"
        File normalRationMedianTsv = "~{outputDir}/~{referenceName}.cobalt.ratio.median.tsv"
        File normalRationPcf = "~{outputDir}/~{referenceName}.cobalt.ratio.pcf"
        File tumorGcMedianTsv = "~{outputDir}/~{tumorName}.cobalt.gc.median.tsv"
        File tumorRatioPcf = "~{outputDir}/~{tumorName}.cobalt.ratio.pcf"
        File tumorRatioTsv = "~{outputDir}/~{tumorName}.cobalt.ratio.tsv.gz"
        Array[File] outputs = [version, normalGcMedianTsv, normalRationMedianTsv,
            normalRationPcf, tumorGcMedianTsv, tumorRatioPcf, tumorRatioTsv]
    }

    runtime {
        docker: "quay.io/pacbio/purple@sha256:9074d56ad46f3d6804f1ae6e45a2ac8c300effc5efe120020d3dcd23edb3b063"
        cpu: threads
        memory: memory
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}

task Purple {
    input {
        String referenceName
        String tumorName
        String outputDir = "purple"
        Array[File]+ amberOutput
        Array[File]+ cobaltOutput
        File? somaticVcf
        File? germlineVcf
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        File? driverGenePanel
        File? somaticHotspots
        File? germlineHotspots
        Float? highlyDiploidPercentage
        Float? somaticMinPuritySpread
        File ensembl_data_dir_tarball
        Float min_purity = 0.08
        Float max_purity = 0.99
        Float min_ploidy = 1
        Float max_ploidy = 8

        Int threads = 8
    }

    String memory = threads * 4 + "GB"
    String javaXmx = threads * 4 - 4 + "G"
    Float file_size = ceil(size(referenceFasta, "GB") + 20)

    command {
        set -euxo pipefail

        mkdir ensembl_data_dir

        tar -xzf ~{ensembl_data_dir_tarball} -C ensembl_data_dir

        java -Xmx~{javaXmx} -jar /app/purple.jar -version

        vcf_var="~{somaticVcf}"

        # Pre-process somatic VCF
        if [[ ! -z $vcf_var ]]; then
            bash /app/process_clairS_VCF.sh ~{somaticVcf} ~{tumorName} ~{referenceName} tmp.vcf.gz

            java -Xmx~{javaXmx} -jar /app/purple.jar \
                -reference ~{referenceName} \
                ~{"-germline_vcf " + germlineVcf} \
                ~{"-germline_hotspots " + germlineHotspots} \
                -tumor ~{tumorName} \
                -output_dir ~{outputDir} \
                -amber ~{sub(amberOutput[0], basename(amberOutput[0]), "")} \
                -cobalt ~{sub(cobaltOutput[0], basename(cobaltOutput[0]), "")} \
                -gc_profile ensembl_data_dir/hmf*/dna/copy_number/GC_profile.*.cnp \
                -somatic_vcf tmp.vcf.gz \
                -ref_genome ~{referenceFasta} \
                -ref_genome_version 38 \
                -ensembl_data_dir ensembl_data_dir/hmf*/common/ensembl_data \
                ~{"-somatic_hotspots " + somaticHotspots} \
                ~{"-run_drivers -driver_gene_panel " + driverGenePanel} \
                ~{"-highly_diploid_percentage " + highlyDiploidPercentage} \
                ~{"-somatic_min_purity_spread " + somaticMinPuritySpread} \
                -threads ~{threads} \
                -circos /usr/bin/circos \
                -max_purity ~{max_purity} \
                -min_purity ~{min_purity} \
                -min_ploidy ~{min_ploidy} \
                -max_ploidy ~{max_ploidy}
        else
            java -Xmx~{javaXmx} -jar /app/purple.jar \
                -reference ~{referenceName} \
                ~{"-germline_vcf " + germlineVcf} \
                ~{"-germline_hotspots " + germlineHotspots} \
                -tumor ~{tumorName} \
                -output_dir ~{outputDir} \
                -amber ~{sub(amberOutput[0], basename(amberOutput[0]), "")} \
                -cobalt ~{sub(cobaltOutput[0], basename(cobaltOutput[0]), "")} \
                -gc_profile ensembl_data_dir/hmf*/dna/copy_number/GC_profile.*.cnp \
                -ref_genome ~{referenceFasta} \
                -ref_genome_version 38 \
                -ensembl_data_dir ensembl_data_dir/hmf*/common/ensembl_data \
                ~{"-somatic_hotspots " + somaticHotspots} \
                ~{"-run_drivers -driver_gene_panel " + driverGenePanel} \
                ~{"-highly_diploid_percentage " + highlyDiploidPercentage} \
                ~{"-somatic_min_purity_spread " + somaticMinPuritySpread} \
                -threads ~{threads} \
                -circos /usr/bin/circos \
                -max_purity ~{max_purity} \
                -min_purity ~{min_purity} \
                -min_ploidy ~{min_ploidy} \
                -max_ploidy ~{max_ploidy}
        fi

        cut -f1,5 ~{outputDir}/~{tumorName}.purple.purity.tsv | tail -n+2 > purity_ploidy.tsv

        rm -rf ensembl_data_dir
    }

    output {
        File? driverCatalogGermlineTsv = "~{outputDir}/~{tumorName}.driver.catalog.germline.tsv"
        File? driverCatalogSomaticTsv = "~{outputDir}/~{tumorName}.driver.catalog.somatic.tsv"
        File purpleCnvGeneTsv = "~{outputDir}/~{tumorName}.purple.cnv.gene.tsv"
        File purpleCnvSomaticTsv = "~{outputDir}/~{tumorName}.purple.cnv.somatic.tsv"
        File purpleGermlineDeletionTsv = "~{outputDir}/~{tumorName}.purple.germline.deletion.tsv"
        File? purpleGermlineVcf = "~{outputDir}/~{tumorName}.purple.germline.vcf.gz"
        File? purpleGermlineVcfIndex = "~{outputDir}/~{tumorName}.purple.germline.vcf.gz.tbi"
        File purplePurityRangeTsv = "~{outputDir}/~{tumorName}.purple.purity.range.tsv"
        File purplePurityTsv = "~{outputDir}/~{tumorName}.purple.purity.tsv"
        File purity_ploidy = "purity_ploidy.tsv"
        File purpleQc = "~{outputDir}/~{tumorName}.purple.qc"
        File purpleSegmentTsv = "~{outputDir}/~{tumorName}.purple.segment.tsv"
        File purpleSomaticClonalityTsv = "~{outputDir}/~{tumorName}.purple.somatic.clonality.tsv"
        File? purpleSomaticHistTsv = "~{outputDir}/~{tumorName}.purple.somatic.hist.tsv"
        File? purpleSomaticVcf = "~{outputDir}/~{tumorName}.purple.somatic.vcf.gz"
        File? purpleSomaticVcfIndex = "~{outputDir}/~{tumorName}.purple.somatic.vcf.gz.tbi"
        File? purpleSvVcf = "~{outputDir}/~{tumorName}.purple.sv.vcf.gz"
        File? purpleSvVcfIndex = "~{outputDir}/~{tumorName}.purple.sv.vcf.gz.tbi"
        File purpleVersion = "~{outputDir}/purple.version"
        File? circosPlot = "~{outputDir}/plot/~{tumorName}.circos.png"
        File? copynumberPlot = "~{outputDir}/plot/~{tumorName}.copynumber.png"
        File? inputPlot = "~{outputDir}/plot/~{tumorName}.input.png"
        File? mapPlot = "~{outputDir}/plot/~{tumorName}.map.png"
        File purityRangePlot = "~{outputDir}/plot/~{tumorName}.purity.range.png"
        File segmentPlot = "~{outputDir}/plot/~{tumorName}.segment.png"
        File? somaticClonalityPlot = "~{outputDir}/plot/~{tumorName}.somatic.clonality.png"
        File? somaticPlot = "~{outputDir}/plot/~{tumorName}.somatic.png"
        File? somaticRainfallPlot = "~{outputDir}/plot/~{tumorName}.somatic.rainfall.png"
        File? circosNormalRatio = "~{outputDir}/circos/~{referenceName}.ratio.circos"
        File? circosBaf = "~{outputDir}/circos/~{tumorName}.baf.circos"
        File? circosConf = "~{outputDir}/circos/~{tumorName}.circos.conf"
        File? circosCnv = "~{outputDir}/circos/~{tumorName}.cnv.circos"
        File? circosIndel = "~{outputDir}/circos/~{tumorName}.indel.circos"
        File? circosInputConf = "~{outputDir}/circos/~{tumorName}.input.conf"
        File? circosLink = "~{outputDir}/circos/~{tumorName}.link.circos"
        File? circosMap = "~{outputDir}/circos/~{tumorName}.map.circos"
        File? circosTumorRatio = "~{outputDir}/circos/~{tumorName}.ratio.circos"
        File? circosSnp = "~{outputDir}/circos/~{tumorName}.snp.circos"
        File? circosGaps = "~{outputDir}/circos/gaps.txt"
        Array[File] outputs = select_all([driverCatalogSomaticTsv, purpleCnvGeneTsv,
            purpleCnvSomaticTsv, purplePurityRangeTsv, purplePurityTsv, purpleQc,
            purpleSegmentTsv, purpleSomaticClonalityTsv, purpleSomaticHistTsv,
            purpleSomaticVcf, purpleSomaticVcfIndex, purpleSvVcf, purpleSvVcfIndex,
            purpleVersion, purpleGermlineVcf, purpleGermlineVcfIndex, driverCatalogGermlineTsv])
        Array[File] plots = select_all([circosPlot, copynumberPlot, inputPlot, mapPlot, purityRangePlot,
            segmentPlot, somaticClonalityPlot, somaticPlot, somaticRainfallPlot])
        Array[File] circos = select_all([circosNormalRatio, circosConf, circosIndel, circosLink,
            circosTumorRatio, circosGaps, circosBaf, circosCnv, circosInputConf, circosMap,
            circosSnp])
    }

    runtime {
        docker: "quay.io/pacbio/purple@sha256:9074d56ad46f3d6804f1ae6e45a2ac8c300effc5efe120020d3dcd23edb3b063"
        cpu: threads
        memory: memory
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}
