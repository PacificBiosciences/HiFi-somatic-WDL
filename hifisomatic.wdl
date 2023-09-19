version 1.1

import "structs.wdl"
import "tasks/alignment.wdl" as alignment
import "tasks/clairs.wdl" as clairs
import "tasks/cnvkit.wdl" as cnvkit
import "tasks/common.wdl" as common
import "tasks/structural_variants.wdl" as structural_variants
import "tasks/phasing.wdl" as phasing
import "tasks/basemod.wdl" as basemod
import "tasks/annotation.wdl" as annotation

workflow hifisomatic {
  input {
    Cohort cohort
    # Files are not aligned by default. If aligned, set skip_align to true
    Boolean skip_align = false
    File ref_fasta
    # Can also define FASTA MMI for pbmm2
    File? ref_fasta_mmi
    # FASTA index refers to the standard faidx. For mmi index use ref_fasta
    File ref_fasta_index
    File ref_bed
    File cnvkit_refflat
    Int cnvkit_threads = 8
    Int pbmm2_threads = 24
    Int samtools_threads = 8
    Int merge_bam_threads = 4
    Int tumor_pileup_mincov = 5
    Int normal_pileup_mincov = 5
    Int cpg_pileup_threads = 8
    Int dss_threads = 16
    # ClairS platform, default hifi_revio
    Boolean call_small_variants = true
    String clairs_platform = "hifi_revio"
    Int clairs_threads = 16
    Int clairs_snv_qual = 2
    Int clairs_indel_qual = 11
    # Sniffles
    File sniffles_trf_bed
    Int sniffles_threads = 8
    File hprc_sniffles_control_vcf
    File hprc_sniffles_control_vcf_index
    # AnnotSV cache can be downloaded using install script from https://lbgi.fr/AnnotSV/.
    # After the database is downloaded, zip the folder and provide the path to the zip file.
    # E.g. by default this is $ANNOTSV/share/AnnotSV
    # If this is not specified in input, AnnotSV will not be run
    File? annotsv_cache
    Int annotsv_threads = 8
    # Default number of threads for misc tasks (4GB per thread assigned for almost all steps)
    Int def_threads = 2
    # Scatter chunk size, default of 75 Mbp
    # Set to 0 to disable chunking
    Int chunk_size = 75000000
    # Strip kinetics or not
    Boolean strip_kinetics = false
    # Calculate DMR?
    Boolean calculate_DMR = true
    # Annotate VCF. If vep_cache is not specified in JSON, VEP will not be run
    # VEP cache can be downloaded from https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_refseq_vep_110_GRCh38.tar.gz
    File? vep_cache
    Int vep_threads = 8
  }

  call common.splitContigs {
    input:
      ref_fasta_index = ref_fasta_index,
      chunk_size = chunk_size,
      threads = def_threads
  }

  scatter (individual in cohort.patients) {
    String patient = individual.patient_names
    Array[File] patient_tumor_bam_files = individual.tumor_bams
    Array[File] patient_normal_bam_files = individual.normal_bams

    if(!skip_align){
      scatter (tumor_bam in patient_tumor_bam_files) {
        if(strip_kinetics) {
          call common.fgbio_strip as strip_tumor_kinetics {
            input:
              bam = tumor_bam,
              threads = def_threads
          }
        }

        call alignment.Align as TumorAlign {
          input:
            sample_name = patient + ".tumor",
            bam_file = select_first([strip_tumor_kinetics.stripped_bam, tumor_bam]),
            ref_fasta = ref_fasta,
            ref_fasta_mmi = ref_fasta_mmi,
            ref_fasta_index = ref_fasta_index,
            threads = pbmm2_threads
        }
      }

      scatter (normal_bam in patient_normal_bam_files) {
        if(strip_kinetics) {
          call common.fgbio_strip as strip_normal_kinetics {
            input:
              bam = normal_bam,
              threads = def_threads
          }
        }
        call alignment.Align as NormalAlign {
          input:
            sample_name = patient + ".normal",
            bam_file = select_first([strip_normal_kinetics.stripped_bam, normal_bam]),
            ref_fasta = ref_fasta,
            ref_fasta_mmi = ref_fasta_mmi,
            ref_fasta_index = ref_fasta_index,
            threads = pbmm2_threads
        }
      }
    }

    call alignment.MergeBams as MergeTumorBams {
      input:
        sample_name = patient + ".tumor",
        bam_files = select_first([TumorAlign.aligned_bam, patient_tumor_bam_files]),
        threads = merge_bam_threads
    }

    call alignment.IndexBam as IndexTumorBam {
      input:
        bam = MergeTumorBams.merged_aligned_bam,
        threads = samtools_threads
    }

    call alignment.MergeBams as MergeNormalBams {
      input:
        sample_name = patient + ".normal",
        bam_files = select_first([NormalAlign.aligned_bam, patient_normal_bam_files]),
        threads = merge_bam_threads
    }

    call alignment.IndexBam as IndexNormalBam {
      input:
        bam = MergeNormalBams.merged_aligned_bam,
        threads = samtools_threads
    }

    call common.mosdepth as MosdepthTumor {
      input:
        pname = patient + ".tumor",
        bam = MergeTumorBams.merged_aligned_bam,
        bam_index = IndexTumorBam.merged_aligned_bam_index_bai,
        threads = def_threads
    }

    call common.mosdepth as MosdepthNormal {
      input:
        pname = patient + ".normal",
        bam = MergeNormalBams.merged_aligned_bam,
        bam_index = IndexNormalBam.merged_aligned_bam_index_bai,
        threads = def_threads
    }
    
    if (call_small_variants) {
          scatter (ctg in splitContigs.contigs) {
          call clairs.ClairS {
            input:
              pname = patient,
              tumor_bam = MergeTumorBams.merged_aligned_bam,
              tumor_bam_index = IndexTumorBam.merged_aligned_bam_index_bai,
              normal_bam = MergeNormalBams.merged_aligned_bam,
              normal_bam_index = IndexNormalBam.merged_aligned_bam_index_bai,
              contig = ctg,
              ref_fasta = ref_fasta,
              ref_fasta_index = ref_fasta_index,
              platform = clairs_platform,
              threads = clairs_threads
            }
          }
          call clairs.gatherClairS {
            input:
              snv_vcf = ClairS.output_snv_vcf,
              indel_vcf = ClairS.output_indel_vcf,
              snv_vcf_index = ClairS.output_snv_vcf_index,
              indel_vcf_index = ClairS.output_indel_vcf_index,
              ref_fasta_index = ref_fasta_index,
              pname = patient,
              threads = samtools_threads,
              snv_qual = clairs_snv_qual,
              indel_qual = clairs_indel_qual
          }

          # Annotate somatic VCF
          if(size(gatherClairS.output_vcf) > 0 && defined(vep_cache)){
            call annotation.vep_annotate {
              input:
                input_vcf = gatherClairS.output_vcf,
                vep_cache = select_first([vep_cache]),
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                pname = patient,
                threads = vep_threads
            }
          }

          call clairs.gatherClairS_germline {
            input:
              tumor_vcf = ClairS.output_snv_vcf_germline_tumor,
              normal_vcf = ClairS.output_snv_vcf_germline_normal,
              pname = patient,
              threads = samtools_threads
          }
          # Phase only if size of VCF is not zero (no variants)
          if(size(gatherClairS_germline.output_tumor_germline_vcf) > 0){
            call phasing.hiphase as phaseTumorBam {
              input:
                bam = MergeTumorBams.merged_aligned_bam,
                bam_index = IndexTumorBam.merged_aligned_bam_index_bai,
                vcf = gatherClairS_germline.output_tumor_germline_vcf,
                vcf_index = gatherClairS_germline.output_tumor_germline_vcf_index,
                pname = patient,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                threads = samtools_threads
            }
            call alignment.IndexBam as IndexTumorBamPhased {
              input:
                bam = phaseTumorBam.hiphase_bam,
                threads = samtools_threads
            }
          }
          if(size(gatherClairS_germline.output_normal_germline_vcf) > 0) {
            call phasing.hiphase as phaseNormalBam {
              input:
                bam = MergeNormalBams.merged_aligned_bam,
                bam_index = IndexNormalBam.merged_aligned_bam_index_bai,
                vcf = gatherClairS_germline.output_normal_germline_vcf,
                vcf_index = gatherClairS_germline.output_normal_germline_vcf_index,
                pname = patient,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                threads = samtools_threads
            }
            call alignment.IndexBam as IndexNormalBamPhased {
              input:
                bam = phaseNormalBam.hiphase_bam,
                threads = samtools_threads
            }
          }
    }

    # # Simple read length and histogram with cramino
    # # Problem is that cramino counts supp alignments too so you end up with
    # # more reads than it should be since multiple alignments are counted
    # call common.cramino as craminoTumor {
    #   input:
    #     bam = select_first([phaseTumorBam.hiphase_bam, MergeTumorBams.merged_aligned_bam]),
    #     bam_index = select_first([IndexTumorBamPhased.merged_aligned_bam_index_bai, IndexTumorBam.merged_aligned_bam_index_bai]),
    #     threads = samtools_threads
    # }

    # call common.cramino as craminoNormal {
    #   input:
    #     bam = select_first([phaseNormalBam.hiphase_bam, MergeNormalBams.merged_aligned_bam]),
    #     bam_index = select_first([IndexNormalBamPhased.merged_aligned_bam_index_bai, IndexNormalBam.merged_aligned_bam_index_bai]),
    #     threads = samtools_threads
    # }

    # Calculate alignment statistics for tumor
    call common.seqkit_bamstats as bamstatsTumor {
      input:
        bam = select_first([phaseTumorBam.hiphase_bam, MergeTumorBams.merged_aligned_bam]),
        bam_index = select_first([IndexTumorBamPhased.merged_aligned_bam_index_bai, IndexTumorBam.merged_aligned_bam_index_bai]),
        threads = samtools_threads
    }

    # Calculate alignment statistics for normal
    call common.seqkit_bamstats as bamstatsNormal {
      input:
        bam = select_first([phaseNormalBam.hiphase_bam, MergeNormalBams.merged_aligned_bam]),
        bam_index = select_first([IndexNormalBamPhased.merged_aligned_bam_index_bai, IndexNormalBam.merged_aligned_bam_index_bai]),
        threads = samtools_threads
    }

    call common.summarize_seqkit_alignment as summarize_tumor_RL {
      input:
        seqkit_alignment_stats = bamstatsTumor.seqkit_alignment_stats,
        threads = def_threads
    }

    call common.summarize_seqkit_alignment as summarize_normal_RL {
      input:
        seqkit_alignment_stats = bamstatsNormal.seqkit_alignment_stats,
        threads = def_threads
    }

    # Pileup CPG
    call common.cpg_pileup as pileup_tumor {
      input:
        pname = patient + ".tumor",
        bam = select_first([phaseTumorBam.hiphase_bam, MergeTumorBams.merged_aligned_bam]),
        bam_index = select_first([IndexTumorBamPhased.merged_aligned_bam_index_bai, IndexTumorBam.merged_aligned_bam_index_bai]),
        reference = ref_fasta,
        reference_index = ref_fasta_index,
        mincov = tumor_pileup_mincov,
        threads = cpg_pileup_threads
    }

    call common.cpg_pileup as pileup_normal {
      input:
        pname = patient + ".normal",
        bam = select_first([phaseNormalBam.hiphase_bam, MergeNormalBams.merged_aligned_bam]),
        bam_index = select_first([IndexNormalBamPhased.merged_aligned_bam_index_bai, IndexNormalBam.merged_aligned_bam_index_bai]),
        reference = ref_fasta,
        reference_index = ref_fasta_index,
        mincov = normal_pileup_mincov,
        threads = cpg_pileup_threads
    }

    # Do differential methylation
    if(calculate_DMR){
      call basemod.DSS_DMR {
        input:
          tumor_bed = pileup_tumor.pileup_combined_beds,
          normal_bed = pileup_normal.pileup_combined_beds,
          pname = patient,
          threads = dss_threads
      }

      call basemod.annotate_DMR {
        input:
          DMR = DSS_DMR.output_DMR,
          pname = patient,
          threads = samtools_threads
      }
    }

    call structural_variants.sniffles as SnifflesNormal {
      input:
        pname = patient + ".normal",
        bam = MergeNormalBams.merged_aligned_bam,
        bam_index = IndexNormalBam.merged_aligned_bam_index_bai,
        trf_bed = sniffles_trf_bed,
        ref_fasta = ref_fasta,
        threads = sniffles_threads
    }

    call structural_variants.sniffles as SnifflesTumor {
      input:
        pname = patient + ".tumor",
        bam = MergeTumorBams.merged_aligned_bam,
        bam_index = IndexTumorBam.merged_aligned_bam_index_bai,
        trf_bed = sniffles_trf_bed,
        ref_fasta = ref_fasta,
        threads = sniffles_threads
    }

    call structural_variants.sniffles_call_snf as SnifflesJoint {
      input:
        pname = patient,
        normal_snf = SnifflesNormal.output_snf,
        tumor_snf = SnifflesTumor.output_snf,
        trf_bed = sniffles_trf_bed,
        ref_fasta = ref_fasta,
        threads = sniffles_threads
    }

    call structural_variants.slivar_select_somatic {
      input:
        sniffles_join_vcf = SnifflesJoint.output_vcf,
        normal_name = patient + ".normal",
        tumor_name = patient + ".tumor",
        threads = def_threads
    }

    call structural_variants.bcftools_filter {
      input:
        slivar_join_vcf = slivar_select_somatic.output_vcf,
        ref_bed = ref_bed,
        tumor_name = patient + ".tumor",
        threads = samtools_threads
    }

    call common.truvari_filter as filterHPRC_sniffles {
      input:
        vcf = bcftools_filter.output_vcf,
        vcf_index = bcftools_filter.output_vcf_index,
        control_vcf = hprc_sniffles_control_vcf,
        control_vcf_index = hprc_sniffles_control_vcf_index,
        threads = samtools_threads,
        truvari_arguments = "-p 0 -s 0 -S 0 --sizemax 100000000 --dup-to-ins"
    }

    call structural_variants.Severus_sv {
      input:
        tumor_bam = select_first([phaseTumorBam.hiphase_bam, MergeTumorBams.merged_aligned_bam]),
        tumor_bam_index = select_first([IndexTumorBamPhased.merged_aligned_bam_index_bai, IndexTumorBam.merged_aligned_bam_index_bai]),
        normal_bam = select_first([phaseNormalBam.hiphase_bam, MergeNormalBams.merged_aligned_bam]),
        normal_bam_index = select_first([IndexNormalBamPhased.merged_aligned_bam_index_bai, IndexNormalBam.merged_aligned_bam_index_bai]),
        phased_vcf = phaseTumorBam.hiphase_vcf,
        pname = patient,
        threads = sniffles_threads
    }

    call common.tabix_vcf as tabixSeverus {
      input:
        vcf = Severus_sv.output_vcf,
        contig_bed = ref_bed,
        threads = samtools_threads
    }

    call common.truvari_filter as filterHPRC_Severus {
      input:
        vcf = tabixSeverus.output_vcf,
        vcf_index = tabixSeverus.output_vcf_index,
        control_vcf = hprc_sniffles_control_vcf,
        control_vcf_index = hprc_sniffles_control_vcf_index,
        threads = samtools_threads,
        truvari_arguments = "-p 0 -s 0 -S 0 --sizemax 100000000 --dup-to-ins"
    }

    if (defined(annotsv_cache)) {
        call annotation.annotsv as annotateSniffles {
          input:
            sv_vcf = filterHPRC_sniffles.output_vcf,
            sv_vcf_index = filterHPRC_sniffles.output_vcf_index,
            annotsv_cache = select_first([annotsv_cache]),
            pname = patient,
            threads = annotsv_threads
      }

      call annotation.annotsv as annotateSeverus {
          input:
            sv_vcf = filterHPRC_Severus.output_vcf,
            sv_vcf_index = filterHPRC_Severus.output_vcf_index,
            annotsv_cache = select_first([annotsv_cache]),
            pname = patient,
            threads = annotsv_threads
      }
    }

    call cnvkit.cnvkit_tumor {
      input:
        tumor_bam = MergeTumorBams.merged_aligned_bam,
        tumor_bam_index = IndexTumorBam.merged_aligned_bam_index_bai,
        normal_bam = MergeNormalBams.merged_aligned_bam,
        normal_bam_index = IndexNormalBam.merged_aligned_bam_index_bai,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        refFlat = cnvkit_refflat,
        pname = patient,
        threads = cnvkit_threads
    }
  }

  output {
    Array[File?] small_variant_vcf = gatherClairS.output_vcf
    Array[File?] small_variant_vcf_annotated = vep_annotate.vep_annotated_vcf
    Array[File?] tumor_germline_small_variant_vcf = gatherClairS_germline.output_tumor_germline_vcf
    Array[File?] normal_germline_small_variant_vcf = gatherClairS_germline.output_normal_germline_vcf
    Array[File] tumor_bams = MergeTumorBams.merged_aligned_bam
    Array[File] tumor_bams_bai = IndexTumorBam.merged_aligned_bam_index_bai
    Array[File?] tumor_bams_phased = phaseTumorBam.hiphase_bam
    Array[File?] tumor_bams_phase_stats = phaseTumorBam.hiphase_stats
    Array[File?] tumor_bams_phased_index = IndexTumorBamPhased.merged_aligned_bam_index_bai
    Array[File] normal_bams = MergeNormalBams.merged_aligned_bam
    Array[File] normal_bams_bai = IndexNormalBam.merged_aligned_bam_index_bai
    Array[File?] normal_bams_phased = phaseNormalBam.hiphase_bam
    Array[File?] normal_bams_phase_stats = phaseNormalBam.hiphase_stats
    Array[File?] normal_bams_phased_index = IndexNormalBamPhased.merged_aligned_bam_index_bai
    Array[Array[File]] pileup_tumor_bed = pileup_tumor.pileup_beds
    Array[Array[File]] pileup_tumor_bw = pileup_tumor.pileup_bigwigs
    Array[Array[File]] pileup_normal_bed = pileup_normal.pileup_beds
    Array[Array[File]] pileup_normal_bw = pileup_normal.pileup_bigwigs
    Array[File?] DMR_results = DSS_DMR.output_DMR
    Array[Array[File]+?] DMR_annotated = annotate_DMR.output_DMR_annotated
    Array[File] sniffles_normal_vcf = SnifflesNormal.output_vcf
    Array[File] sniffles_tumor_vcf = SnifflesTumor.output_vcf
    Array[File] sniffles_joint_vcf = SnifflesJoint.output_vcf
    Array[File] sniffles_somatic_vcf = bcftools_filter.output_vcf
    Array[File] sniffles_somatic_vcf_index = bcftools_filter.output_vcf_index
    Array[File] sniffles_somatic_vcf_filterHPRC = filterHPRC_sniffles.output_vcf
    Array[File] sniffles_somatic_vcf_filterHPRC_index = filterHPRC_sniffles.output_vcf_index
    Array[File] Severus_vcf = Severus_sv.output_vcf
    Array[File] Severus_filterHPRC_vcf = filterHPRC_Severus.output_vcf
    Array[File] Severus_filterHPRC_vcf_index = filterHPRC_Severus.output_vcf_index
    Array[Array[File]+] cnvkit_output = cnvkit_tumor.cnvkit_output
    Array[File?] AnnotatedSnifflesSV = annotateSniffles.annotsv_annotated_tsv
    Array[File?] AnnotatedSeverusSV = annotateSeverus.annotsv_annotated_tsv
    Array[File] mosdepth_tumor_bed = MosdepthTumor.output_bed
    Array[File] mosdepth_tumor_bed_index = MosdepthTumor.output_bed_index
    Array[File] mosdepth_tumor_summary = MosdepthTumor.output_summary
    Array[File] mosdepth_normal_bed = MosdepthNormal.output_bed
    Array[File] mosdepth_normal_bed_index = MosdepthNormal.output_bed_index
    Array[File] mosdepth_normal_summary = MosdepthNormal.output_summary
    Array[File] overall_tumor_alignment_stats = bamstatsTumor.seqkit_bam_stats
    Array[File] overall_normal_alignment_stats = bamstatsNormal.seqkit_bam_stats
    Array[File] per_alignment_tumor_stats = bamstatsTumor.seqkit_alignment_stats
    Array[File] per_alignment_normal_stats = bamstatsNormal.seqkit_alignment_stats
    Array[File] aligned_RL_summary_tumor = summarize_tumor_RL.output_summary
    Array[File] aligned_RL_summary_normal = summarize_normal_RL.output_summary
    # Array[File] cramino_tumor_bamstats = craminoTumor.output_cramino_stats
    # Array[File] cramino_normal_bamstats = craminoNormal.output_cramino_stats
  }
}