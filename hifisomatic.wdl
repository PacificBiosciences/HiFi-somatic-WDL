version 1.0

import "structs.wdl"
import "tasks/alignment.wdl" as alignment
import "tasks/clairs.wdl" as clairs
import "tasks/cnvkit.wdl" as cnvkit
import "tasks/common.wdl" as common
import "tasks/structural_variants.wdl" as structural_variants
import "tasks/phasing.wdl" as phasing
import "tasks/basemod.wdl" as basemod
import "tasks/annotation.wdl" as annotation
import "tasks/prioritization.wdl" as prioritization
import "tasks/clonality.wdl" as clonality
import "tasks/deepsomatic.wdl" as deepsomatic
import "tasks/biomarker.wdl" as biomarker 

workflow hifisomatic {
  input {
    Cohort cohort
    # Files are not aligned by default. If aligned, set skip_align to true
    Boolean skip_align = false
    # Strip kinetics or not
    Boolean strip_kinetics = false
    File ref_fasta
    File ref_fasta_dict
    File ref_gff
    # FASTA index refers to the standard faidx. For mmi index use ref_fasta
    File ref_fasta_index
    # Can also define FASTA MMI for pbmm2
    File? ref_fasta_mmi
    File ref_bed
    File cnvkit_refflat
    Int cnvkit_threads = 8
    Int pbmm2_threads = 24
    Int samtools_threads = 8
    Int hiphase_threads = 24
    Int merge_bam_threads = 4
    Int tumor_pileup_mincov = 5
    Int normal_pileup_mincov = 5
    Int cpg_pileup_threads = 8
    Int dss_threads = 16
    # Call small variants with deepsomatic
    Boolean use_deepsomatic = true
    # ClairS platform, default hifi_revio
    Boolean call_small_variants = true
    String clairs_platform = "hifi_revio_ssrs"
    Int clairs_threads = 16
    Int clairs_snv_qual = 2
    Int clairs_indel_qual = 11
    ## Deepsomatic threads per task
    Int deepsomatic_threads = 16
    # Mutational signature max delta for MutationalPattern fitting
    Float mutsig_max_delta = 0.004
    # SV-related 
    File trf_bed
    Int sv_threads = 24
    Int wakhan_threads = 24
    File control_vcf
    File control_vcf_index
    File? severus_pon_tsv
    # Minimum reads support for Severus
    Int severus_min_reads = 3
    # AnnotSV cache can be downloaded using install script from https://lbgi.fr/AnnotSV/.
    # After the database is downloaded, zip the folder and provide the path to the zip file.
    # E.g. by default this is $ANNOTSV/share/AnnotSV
    # If this is not specified in input, AnnotSV will not be run
    File? annotsv_cache
    Int annotsv_threads = 8
    # Default number of threads for misc tasks (4GB per thread assigned for almost all steps)
    Int def_threads = 2
    # Scatter small variants calling into equal chunk per chromosome to make use of multiple nodes. Default of 75 Mbp per chromosome (total of 42 chunks for hg38)
    Int chunk_size = 75000000
    # Calculate DMR?
    Boolean calculate_DMR = true
    # Annotate VCF. If vep_cache is not specified in JSON, VEP will not be run
    # VEP cache can be downloaded from https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_refseq_vep_110_GRCh38.tar.gz
    File? vep_cache
    Int vep_threads = 8
    # Annotate germline variants?
    Boolean annotate_germline = false
    # Run SAVANA sv call too?
    Boolean run_savana = false
    # Minimum number of CG to prioritize
    Int ncg_to_filter = 50
    # Amber, Cobalt and Purple
    # Higher value = smoother CNV, but lower resolution
    Int pcf_gamma_value = 1000
    File? ensembl_data_dir_tarball
    Int hmftools_threads = 8
    # Visualization script
    File report_vis_script
  }

  scatter (individual in cohort.patients) {
    String patient = individual.patient_names
    Array[File] patient_tumor_bam_files = individual.tumor_bams
    Array[File] patient_normal_bam_files = individual.normal_bams
    File ref_to_use = select_first([ref_fasta_mmi, ref_fasta])

    call alignment.align_all_bams as align_tumor {
      input:
        patient = patient,
        bam_files = patient_tumor_bam_files,
        ref_fasta = ref_to_use,
        ref_fasta_index = ref_fasta_index,
        pbmm2_threads = pbmm2_threads,
        merge_bam_threads = merge_bam_threads,
        samtools_threads = samtools_threads,
        skip_align = skip_align,
        strip_kinetics = strip_kinetics,
        suffix = "tumor"
    }

    call alignment.align_all_bams as align_normal {
      input:
        patient = patient,
        bam_files = patient_normal_bam_files,
        ref_fasta = ref_to_use,
        ref_fasta_index = ref_fasta_index,
        pbmm2_threads = pbmm2_threads,
        merge_bam_threads = merge_bam_threads,
        samtools_threads = samtools_threads,
        skip_align = skip_align,
        strip_kinetics = strip_kinetics,
        suffix = "normal"
    }

    call common.mosdepth as MosdepthTumor {
      input:
        pname = patient + ".tumor",
        bam = align_tumor.bam_final,
        bam_index = align_tumor.bam_final_index,
        threads = def_threads
    }

    call common.mosdepth as MosdepthNormal {
      input:
        pname = patient + ".normal",
        bam = align_normal.bam_final,
        bam_index = align_normal.bam_final_index,
        threads = def_threads
    }

    call common.split_contigs {
      input:
        ref_fasta_index = ref_fasta_index,
        chunk_size = chunk_size,
        threads = def_threads
    }
    
    if (call_small_variants) {
          if (use_deepsomatic){
              call deepsomatic.run_deepsomatic {
                input:
                  tumor_bam = align_tumor.bam_final,
                  tumor_bam_index = align_tumor.bam_final_index,
                  normal_bam = align_normal.bam_final,
                  normal_bam_index = align_normal.bam_final_index,
                  ref_fasta = ref_fasta,
                  ref_fasta_index = ref_fasta_index,
                  threads = deepsomatic_threads,
                  pname = patient,
                  contigs = split_contigs.contigs
                }
          }

          if (!use_deepsomatic){
            scatter (ctg in split_contigs.contigs) {
              call clairs.ClairS {
                input:
                  pname = patient,
                  tumor_bam = align_tumor.bam_final,
                  tumor_bam_index = align_tumor.bam_final_index,
                  normal_bam = align_normal.bam_final,
                  normal_bam_index = align_normal.bam_final_index,
                  contig = ctg,
                  ref_fasta = ref_fasta,
                  ref_fasta_index = ref_fasta_index,
                  platform = clairs_platform,
                  threads = clairs_threads
                }
            }

            call clairs.gather_ClairS {
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

            call clairs.gather_ClairS_germline {
              input:
                tumor_vcf = ClairS.output_snv_vcf_germline_tumor,
                normal_vcf = ClairS.output_snv_vcf_germline_normal,
                pname = patient,
                threads = samtools_threads
            }
          }

          # Mutational signature
          call common.mutationalpattern {
            input:
              vcf = select_first([run_deepsomatic.deepsomatic_vcf, gather_ClairS.output_vcf]),
              pname = patient,
              max_delta = mutsig_max_delta,
              threads = samtools_threads
          }

          # Phase only if size of VCF is not zero (no variants)
          if(size(select_first([run_deepsomatic.clair3_tumor_vcf, gather_ClairS_germline.output_tumor_germline_vcf])) > 0){
              call phasing.hiphase_with_somatic as phaseTumorBam {
                input:
                  bam = align_tumor.bam_final,
                  bam_index = align_tumor.bam_final_index,
                  vcf = select_first([run_deepsomatic.clair3_tumor_vcf, gather_ClairS_germline.output_tumor_germline_vcf]),
                  vcf_index = select_first([run_deepsomatic.clair3_tumor_vcf_tbi, gather_ClairS_germline.output_tumor_germline_vcf_index]),
                  somatic_SNP_indel_vcf = select_first([run_deepsomatic.deepsomatic_vcf, gather_ClairS.output_vcf]),
                  somatic_SNP_indel_vcf_index = select_first([run_deepsomatic.deepsomatic_vcf_tbi, gather_ClairS.output_vcf_index]),
                  pname = patient,
                  ref_fasta = ref_fasta,
                  ref_fasta_index = ref_fasta_index,
                  threads = hiphase_threads
              }
            
            # Annotate somatic VCF
            if(size(phaseTumorBam.hiphase_somatic_small_variants_vcf) > 0 && defined(vep_cache)){
              call annotation.bcftools_norm as normSomatic {
                input:
                  input_vcf = phaseTumorBam.hiphase_somatic_small_variants_vcf,
                  ref_fasta = ref_fasta,
                  ref_fasta_index = ref_fasta_index,
                  threads = samtools_threads
              }
              call annotation.vep_annotate as annotateSomatic {
                input:
                  input_vcf = normSomatic.norm_vcf,
                  vep_cache = select_first([vep_cache]),
                  ref_fasta = ref_fasta,
                  ref_fasta_index = ref_fasta_index,
                  pname = patient,
                  threads = vep_threads
              }

              # Once annotated, calculate TMB
              call biomarker.tmb_estimate as tmb_estimate {
                input:
                  vcf = annotateSomatic.vep_annotated_vcf,
                  coverage = MosdepthTumor.output_bed,
                  pname = patient
              }

              call prioritization.prioritize_small_variants as prioritizeSomatic {
                input:
                  vep_annotated_vcf = annotateSomatic.vep_annotated_vcf,
                  threads = samtools_threads,
                  pname = patient
              }
            }

          }
          if(size(select_first([run_deepsomatic.clair3_normal_vcf, gather_ClairS_germline.output_normal_germline_vcf])) > 0) {
            call phasing.hiphase as phaseNormalBam {
              input:
                bam = align_normal.bam_final,
                bam_index = align_normal.bam_final_index,
                vcf = select_first([run_deepsomatic.clair3_normal_vcf, gather_ClairS_germline.output_normal_germline_vcf]),
                vcf_index = select_first([run_deepsomatic.clair3_normal_vcf_tbi, gather_ClairS_germline.output_normal_germline_vcf_index]),
                pname = patient,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                threads = hiphase_threads
            }

          # Annotate germline VCF
           if(size(phaseNormalBam.hiphase_vcf) > 0 && defined(vep_cache) && annotate_germline){
              call annotation.vep_annotate as annotateGermline {
                input:
                  input_vcf = phaseNormalBam.hiphase_vcf,
                  vep_cache = select_first([vep_cache]),
                  ref_fasta = ref_fasta,
                  ref_fasta_index = ref_fasta_index,
                  pname = patient,
                  threads = vep_threads
              }
            }
          }
    }
    # Calculate alignment statistics for tumor
    call common.seqkit_bamstats as bamstatsTumor {
      input:
        bam = select_first([phaseTumorBam.hiphase_bam, align_tumor.bam_final]),
        bam_index = select_first([phaseTumorBam.hiphase_bam_index, align_tumor.bam_final_index]),
        threads = samtools_threads
    }

    # Calculate alignment statistics for normal
    call common.seqkit_bamstats as bamstatsNormal {
      input:
        bam = select_first([phaseNormalBam.hiphase_bam, align_normal.bam_final]),
        bam_index = select_first([phaseNormalBam.hiphase_bam_index, align_normal.bam_final_index]),
        threads = samtools_threads
    }

    call common.summarize_seqkit_alignment as summarize_tumor_RL {
      input:
        seqkit_alignment_stats = bamstatsTumor.seqkit_alignment_stats,
        threads = 8
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
        bam = select_first([phaseTumorBam.hiphase_bam, align_tumor.bam_final]),
        bam_index = select_first([phaseTumorBam.hiphase_bam_index, align_tumor.bam_final_index]),
        reference = ref_fasta,
        reference_index = ref_fasta_index,
        mincov = tumor_pileup_mincov,
        threads = cpg_pileup_threads
    }

    call common.cpg_pileup as pileup_normal {
      input:
        pname = patient + ".normal",
        bam = select_first([phaseNormalBam.hiphase_bam, align_normal.bam_final]),
        bam_index = select_first([phaseNormalBam.hiphase_bam_index, align_normal.bam_final_index]),
        reference = ref_fasta,
        reference_index = ref_fasta_index,
        mincov = normal_pileup_mincov,
        threads = cpg_pileup_threads
    }

    # owl MSI score
    call biomarker.owl_msi as msi_score {
      input:
        bam_file = select_first([phaseTumorBam.hiphase_bam, align_tumor.bam_final]),
        bam_index = select_first([phaseTumorBam.hiphase_bam_index, align_tumor.bam_final_index]),
        sample_name = patient
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

      call prioritization.prioritize_dmr_intogen {
        input:
          dmr_files = annotate_DMR.output_DMR_annotated,
          threads = samtools_threads,
          ncg_to_filter = ncg_to_filter,
          pname = patient
      }
    }

    if(call_small_variants){
      call structural_variants.Severus_sv as phased_severus {
        input:
          tumor_bam = select_first([phaseTumorBam.hiphase_bam, align_tumor.bam_final]),
          tumor_bam_index = select_first([phaseTumorBam.hiphase_bam_index, align_tumor.bam_final_index]),
          normal_bam = select_first([phaseNormalBam.hiphase_bam, align_normal.bam_final]),
          normal_bam_index = select_first([phaseNormalBam.hiphase_bam_index, align_normal.bam_final_index]),
          trf_bed = trf_bed,
          phased_vcf = phaseNormalBam.hiphase_vcf,
          pname = patient,
          threads = sv_threads,
          min_supp_reads = severus_min_reads,
          PON_tsv = severus_pon_tsv
      }

      # CNV with Wakhan
      call structural_variants.wakhan {
        input:
          pname = patient,
          tumor_bam = select_first([phaseTumorBam.hiphase_bam, align_tumor.bam_final]),
          tumor_bam_index = select_first([phaseTumorBam.hiphase_bam_index, align_tumor.bam_final_index]),
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          severus_sv_vcf = phased_severus.output_vcf,
          normal_germline_vcf = phaseNormalBam.hiphase_vcf,
          tumor_germline_vcf = select_first([phaseTumorBam.hiphase_vcf]),
          threads = wakhan_threads
      }

      # Chord HRD prediction
      call annotation.chord_hrd {
        input:
          small_variant_vcf = select_first([run_deepsomatic.deepsomatic_vcf, gather_ClairS.output_vcf]),
          sv_vcf = phased_severus.output_vcf,
          pname = patient
      }
      
      # Call SAVANA?
      if(run_savana){
        call structural_variants.SAVANA_sv {
          input:
            pname = patient,
            tumor_bam = select_first([phaseTumorBam.hiphase_bam, align_tumor.bam_final]),
            tumor_bam_index = select_first([phaseTumorBam.hiphase_bam_index, align_tumor.bam_final_index]),
            normal_bam = select_first([phaseNormalBam.hiphase_bam, align_normal.bam_final]),
            normal_bam_index = select_first([phaseNormalBam.hiphase_bam_index, align_normal.bam_final_index]),
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            phased_vcf = select_first([phaseNormalBam.hiphase_vcf]),
            min_supp_reads = 3,
            threads = pbmm2_threads
        }
      }
    }

    if(!call_small_variants){
      call structural_variants.Severus_sv as unphased_severus {
        input:
          tumor_bam = select_first([phaseTumorBam.hiphase_bam, align_tumor.bam_final]),
          tumor_bam_index = select_first([phaseTumorBam.hiphase_bam_index, align_tumor.bam_final_index]),
          normal_bam = select_first([phaseNormalBam.hiphase_bam, align_normal.bam_final]),
          normal_bam_index = select_first([phaseNormalBam.hiphase_bam_index, align_normal.bam_final_index]),
          trf_bed = trf_bed,
          pname = patient,
          threads = sv_threads,
          min_supp_reads = severus_min_reads,
          PON_tsv = severus_pon_tsv
      }
    }

    call common.tabix_vcf as tabixSeverus {
      input:
        vcf = select_first([phased_severus.output_vcf, unphased_severus.output_vcf]),
        contig_bed = ref_bed,
        threads = samtools_threads
    }

    call common.svpack_filter_annotated as filter_Severus {
      input:
        sv_vcf = tabixSeverus.output_vcf,
        population_vcfs = [control_vcf],
        population_vcf_indices = [control_vcf_index],
        gff = ref_gff
    }

    # Recover missing mate
    call common.recover_mate_bnd {
      input:
        sv_vcf_original = select_first([phased_severus.output_vcf, unphased_severus.output_vcf]),
        sv_svpack_filtered = filter_Severus.output_vcf,
    }

    # Do circos plot
    call structural_variants.circos_BND {
      input:
        sv_vcf = recover_mate_bnd.output_vcf,
        pname = patient
    }

    if (defined(annotsv_cache)) {
      call annotation.annotsv as annotateSeverus {
          input:
            sv_vcf = recover_mate_bnd.output_vcf,
            sv_vcf_index = recover_mate_bnd.output_vcf_index,
            annotsv_cache = select_first([annotsv_cache]),
            pname = patient,
            threads = annotsv_threads
      }

      call prioritization.prioritize_sv_intogen as prioritize_Severus {
        input:
          annotSV_tsv = annotateSeverus.annotsv_annotated_tsv,
          threads = samtools_threads
      }
    }

    call cnvkit.cnvkit_tumor {
      input:
        tumor_bam = align_tumor.bam_final,
        tumor_bam_index = align_tumor.bam_final_index,
        normal_bam = align_normal.bam_final,
        normal_bam_index = align_normal.bam_final_index,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        refFlat = cnvkit_refflat,
        pname = patient,
        threads = cnvkit_threads
    }

    if(defined(ensembl_data_dir_tarball)){
      call clonality.Amber {
        input:
          referenceName = patient + ".normal",
          referenceBam = align_normal.bam_final,
          referenceBamIndex = align_normal.bam_final_index,
          tumorName = patient + ".tumor",
          tumorBam = align_tumor.bam_final,
          tumorBamIndex = align_tumor.bam_final_index,
          ensembl_data_dir_tarball = select_first([ensembl_data_dir_tarball]),
          referenceFasta = ref_fasta,
          referenceFastaFai = ref_fasta_index,
          referenceFastaDict = ref_fasta_dict,
          threads = hmftools_threads
      }

      call clonality.Cobalt {
        input:
          referenceName = patient + ".normal",
          referenceBam = align_normal.bam_final,
          referenceBamIndex = align_normal.bam_final_index,
          tumorName = patient + ".tumor",
          tumorBam = align_tumor.bam_final,
          tumorBamIndex = align_tumor.bam_final_index,
          referenceFasta = ref_fasta,
          referenceFastaFai = ref_fasta_index,
          referenceFastaDict = ref_fasta_dict,
          ensembl_data_dir_tarball = select_first([ensembl_data_dir_tarball]),
          threads = hmftools_threads,
          pcf_gamma = pcf_gamma_value
      }

      # DeepSomatic doesn't annotate germline count, so don't supply
      # somatic VCF if using deepsomatic
      if (use_deepsomatic){
        call clonality.Purple as purple_nosomatic {
          input:
            referenceName = patient + ".normal",
            tumorName = patient + ".tumor",
            outputDir = "purple",
            amberOutput = Amber.outputs,
            cobaltOutput = Cobalt.outputs,
            referenceFasta = ref_fasta,
            referenceFastaFai = ref_fasta_index,
            referenceFastaDict = ref_fasta_dict,
            ensembl_data_dir_tarball = select_first([ensembl_data_dir_tarball]),
            threads = hmftools_threads
        }
      }
      if (!use_deepsomatic){
        call clonality.Purple as purple_withsomatic {
          input:
            referenceName = patient + ".normal",
            tumorName = patient + ".tumor",
            outputDir = "purple",
            amberOutput = Amber.outputs,
            cobaltOutput = Cobalt.outputs,
            somaticVcf = select_first([gather_ClairS.output_vcf]),
            referenceFasta = ref_fasta,
            referenceFastaFai = ref_fasta_index,
            referenceFastaDict = ref_fasta_dict,
            ensembl_data_dir_tarball = select_first([ensembl_data_dir_tarball]),
            threads = hmftools_threads
        }
      }

      # Recall CNVKit major and minor CN
      if (call_small_variants && size(select_first([run_deepsomatic.clair3_normal_vcf, gather_ClairS_germline.output_normal_germline_vcf])) > 0 && defined(ensembl_data_dir_tarball)) {
        call cnvkit.merge_germline as mergeGermline {
          input:
            tumor_germline_vcf = select_first([run_deepsomatic.clair3_tumor_vcf, gather_ClairS_germline.output_tumor_germline_vcf]),
            tumor_germline_vcf_tbi = select_first([run_deepsomatic.clair3_tumor_vcf_tbi, gather_ClairS_germline.output_tumor_germline_vcf_index]),
            normal_germline_vcf = select_first([run_deepsomatic.clair3_normal_vcf, gather_ClairS_germline.output_normal_germline_vcf]),
            normal_germline_vcf_tbi = select_first([run_deepsomatic.clair3_normal_vcf_tbi, gather_ClairS_germline.output_normal_germline_vcf_index]),
            pname = patient,
            threads = samtools_threads
        }
        
        call cnvkit.cnvkit_recall {
          input:
            cnvkit_cns = cnvkit_tumor.cnvkit_output_cns,
            merged_germline_heterozygous_vcf = mergeGermline.merged_germline_heterozygous_vcf,
            merged_germline_heterozygous_vcf_tbi = mergeGermline.merged_germline_heterozygous_vcf_tbi,
            threads = samtools_threads,
            purity_ploidy = select_first([purple_withsomatic.purity_ploidy, purple_nosomatic.purity_ploidy]),
            pname = patient
        }

        if(defined(annotsv_cache)) {
          call prioritization.report_sample {
            input:
              annotated_small_variant_tsv = select_first([prioritizeSomatic.vep_annotated_tsv]),
              intogen_small_var_tsv = select_first([prioritizeSomatic.vep_annotated_tsv_intogenCCG]),
              sv_intogen_tsv = select_first([prioritize_Severus.annotSV_intogen_tsv]),
              sv_vcf = recover_mate_bnd.output_vcf,
              cnv = select_first([purple_withsomatic.purpleCnvSomaticTsv, purple_nosomatic.purpleCnvSomaticTsv]),
              purity_ploidy = select_first([purple_withsomatic.purplePurityTsv, purple_nosomatic.purplePurityTsv]),
              mosdepth_tumor = MosdepthTumor.output_summary,
              mosdepth_normal = MosdepthNormal.output_summary,
              mutsig_tsv = select_first([mutationalpattern.mutsig_tsv]),
              mut_reconstructed_sigs = select_first([mutationalpattern.recon_sig]),
              dmr_intogen_tsv = select_first([prioritize_dmr_intogen.promoter_file]),
              chord_file = select_first([chord_hrd.chord_prediction]),
              circos_png = circos_BND.circos_png,
              pname = patient,
              vis_file = report_vis_script,
              cnv_tool = "Purple",
              owl_score = msi_score.owl_score_output,
              tmb_estimate_file = select_first([tmb_estimate.tmb_estimate_output]),
              tmb_estimate_gencode_file = select_first([tmb_estimate.tmb_estimate_output_gencode]),
              sv_pairs_file = circos_BND.known_genes_pairs
          }
      }
      }
    }
  }

  output {
    Array[File?] small_variant_vcf = phaseTumorBam.hiphase_somatic_small_variants_vcf
    Array[File?] small_variant_vcf_annotated = annotateSomatic.vep_annotated_vcf
    Array[File?] small_variant_tsv_annotated = prioritizeSomatic.vep_annotated_tsv
    Array[File?] small_variant_tsv_CCG = prioritizeSomatic.vep_annotated_tsv_intogenCCG
    Array[File?] small_variant_tmb_estimate = tmb_estimate.tmb_estimate_output
    Array[File?] small_variant_tmb_estimate_gencode = tmb_estimate.tmb_estimate_output_gencode
    Array[Array[File]+?] mutsig_SNV = mutationalpattern.mutsig_output
    Array[File?] mutsig_SNV_profile = mutationalpattern.mut_profile
    Array[File?] chord_hrd_prediction = chord_hrd.chord_prediction
    Array[File?] tumor_germline_small_variant_vcf = select_first([phaseTumorBam.hiphase_vcf, run_deepsomatic.clair3_tumor_vcf, gather_ClairS_germline.output_tumor_germline_vcf])
    Array[File?] normal_germline_small_variant_vcf = select_first([phaseNormalBam.hiphase_vcf, run_deepsomatic.clair3_normal_vcf, gather_ClairS_germline.output_normal_germline_vcf])
    Array[File?] normal_germline_small_variant_vcf_annotated = annotateGermline.vep_annotated_vcf
    Array[File] tumor_bams = align_tumor.bam_final
    Array[File] tumor_bams_bai = align_tumor.bam_final_index
    Array[File?] tumor_bams_hiphase = phaseTumorBam.hiphase_bam
    Array[Array[File]+?] tumor_bams_hiphase_stats = phaseTumorBam.hiphase_stats_summary
    Array[File?] tumor_bams_phased_index = select_first([phaseTumorBam.hiphase_bam_index])
    Array[File] normal_bams = align_normal.bam_final
    Array[File] normal_bams_bai = align_normal.bam_final_index
    Array[File?] normal_bams_phased = phaseNormalBam.hiphase_bam
    Array[File?] normal_bams_phase_stats = phaseNormalBam.hiphase_stats
    Array[File?] normal_bams_phased_index = phaseNormalBam.hiphase_bam_index
    Array[Array[File]] pileup_tumor_bed = pileup_tumor.pileup_beds
    Array[Array[File]] pileup_tumor_bw = pileup_tumor.pileup_bigwigs
    Array[Array[File]] pileup_normal_bed = pileup_normal.pileup_beds
    Array[Array[File]] pileup_normal_bw = pileup_normal.pileup_bigwigs
    Array[File] owl_msi_profile = msi_score.owl_profile_output
    Array[File] owl_msi_score = msi_score.owl_score_output
    Array[File] owl_msi_motif_counts = msi_score.owl_motif_counts
    Array[File?] DMR_results = DSS_DMR.output_DMR
    Array[Array[File]+?] DMR_annotated = annotate_DMR.output_DMR_annotated
    Array[Array[File]+?] DMR_annotated_CCG = prioritize_dmr_intogen.DMR_nCGFilter_CCG
    Array[File] Severus_vcf = tabixSeverus.output_vcf
    Array[File] Severus_filtered_vcf = recover_mate_bnd.output_vcf
    Array[File] Severus_filtered_vcf_index = recover_mate_bnd.output_vcf_index
    Array[File] SV_circos = circos_BND.circos_plot
    Array[File] SV_mitelman_fusions = circos_BND.mitelman_fusions
    Array[File] SV_known_genes_pairs = circos_BND.known_genes_pairs
    Array[File?] Severus_breakpoint_cluster = select_first([phased_severus.output_breakpoint_clusters, unphased_severus.output_breakpoint_clusters])
    Array[File?] Severus_breakpoint_cluster_all = select_first([phased_severus.output_breakpoint_clusters_all, unphased_severus.output_breakpoint_clusters_all])
    Array[File?] Severus_cluster_plots = select_first([phased_severus.output_somatic_sv_plots, unphased_severus.output_somatic_sv_plots])
    Array[Array[File]+] cnvkit_output = cnvkit_tumor.cnvkit_output
    Array[File?] cnvkit_cns_with_major_minor_CN = cnvkit_recall.cnvkit_cns_with_major_minor_CN
    Array[File?] AnnotatedSeverusSV = annotateSeverus.annotsv_annotated_tsv
    Array[File?] AnnotatedSeverusSV_intogen = prioritize_Severus.annotSV_intogen_tsv
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
    Array[Array[File]+?] Amber_outputs = Amber.outputs
    Array[Array[File]+?] Cobalt_outputs = Cobalt.outputs
    Array[Array[File]+?] Purple_outputs = select_first([purple_nosomatic.outputs, purple_withsomatic.outputs])
    Array[Array[File]+?] Purple_plots = select_first([purple_nosomatic.plots, purple_withsomatic.plots])
    Array[Array[File]+?] wakhan_cnv = wakhan.wakhan_output
    Array[File?] report = report_sample.summary_report
    Array[File?] SAVANA_sv_vcf = SAVANA_sv.savana_vcf
  }
}
