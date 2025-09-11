version 1.0

task tabix_vcf {
  input {
    File vcf
    File contig_bed
    Int threads
  }

  Float file_size = ceil(size(vcf, "GB") + size(contig_bed, "GB") + 10)
  
  command <<<
    set -euxo pipefail
    echo "indexing ~{vcf}"
    # Get rid of SVLEN=0 and change <INS> to INS
    # to avoid Truvari error problem
    # sed 's/SVLEN=0;//g' ~{vcf} > tmp.vcf
    # sed -i 's/<INS>/INS/g' tmp.vcf

    bcftools --version

    # sort first in case input is not sorted
    bcftools sort \
      -Oz -o tmp.vcf.gz \
      ~{vcf}

    tabix -p vcf tmp.vcf.gz

    # Filter out contigs that are not in the contig bed file
    bcftools view \
      -R ~{contig_bed} \
      -Oz -o ~{basename(vcf) + ".gz"} \
      tmp.vcf.gz

    tabix -p vcf ~{basename(vcf) + ".gz"}
    rm -f tmp.vcf.gz tmp.vcf.gz.tbi
  >>>

  output {
    File output_vcf = basename(vcf) + ".gz"
    File output_vcf_index = basename(vcf) + ".gz.tbi"
  }

  runtime {
    docker: "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

# This is used to filter SV VCF with control VCF (Germline variants
# for all samples called with Sniffles)
task truvari_filter {
  input {
    File vcf
    File vcf_index
    File control_vcf
    File control_vcf_index
    Int threads
    String? truvari_arguments
  }

  Float file_size = ceil(size(vcf, "GB") + size(control_vcf, "GB") + 10)

  command <<<
    set -euxo pipefail

    truvari version

    truvari bench \
        -b ~{control_vcf} \
        -c ~{vcf} \
        -o truvari_filter \
        ~{truvari_arguments}

    mv truvari_filter/fp.vcf.gz ~{sub(basename(vcf), "\\.vcf.gz", ".filtered.vcf.gz")}
    mv truvari_filter/fp.vcf.gz.tbi ~{sub(basename(vcf), "\\.vcf.gz", ".filtered.vcf.gz.tbi")}
  >>>

  output {
    File output_vcf = sub(basename(vcf), "\\.vcf.gz", ".filtered.vcf.gz")
    File output_vcf_index = sub(basename(vcf), "\\.vcf.gz", ".filtered.vcf.gz.tbi")
  }

  runtime {
    docker: "quay.io/biocontainers/truvari:4.0.0--pyhdfd78af_0"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

# Use svpack filtering instead of truvari, copied from
# human-wgs-wdl
task svpack_filter_annotated {
	input {
		File sv_vcf

		Array[File] population_vcfs
		Array[File] population_vcf_indices

    Float? svlen

		File gff
	}

	String sv_vcf_basename = basename(sv_vcf, ".vcf.gz")
	Int file_size = ceil(size(sv_vcf, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		echo "svpack version:"
		cat /opt/svpack/.git/HEAD

		svpack \
			filter \
			--pass-only \
			~{"--min-svlen " + svlen} \
			~{sv_vcf} \
			~{sep=' ' prefix('| svpack match -v - ', population_vcfs)} \
		| svpack \
			consequence \
			- \
			~{gff} \
		| svpack \
			tagzygosity \
			- \
		> ~{sv_vcf_basename}.svpack.vcf

		bgzip --version

		bgzip ~{sv_vcf_basename}.svpack.vcf

		tabix --version

		tabix -p vcf ~{sv_vcf_basename}.svpack.vcf.gz
	>>>

	output {
		File output_vcf = "~{sv_vcf_basename}.svpack.vcf.gz"
		File output_vcf_index = "~{sv_vcf_basename}.svpack.vcf.gz.tbi"
	}

	runtime {
		docker: "quay.io/pacbio/svpack@sha256:a680421cb517e1fa4a3097838719a13a6bd655a5e6980ace1b03af9dd707dd75"
		cpu: 4
    memory: "16 GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
	}
}

# Use bcftools to add missing BND mate back after svpack filtering
task recover_mate_bnd {
  input {
    File sv_vcf_original
    File sv_svpack_filtered
  }

  Float file_size = ceil(size(sv_vcf_original, "GB") + size(sv_svpack_filtered, "GB") + 10)

  command <<<
  set -euxo pipefail

  bcftools --version

  # Copy the VCF here and index them
  cp ~{sv_svpack_filtered} ~{basename(sv_svpack_filtered)}
  tabix ~{basename(sv_svpack_filtered)}

  # For original VCF, bgzip if not already, then index
  cp ~{sv_vcf_original} ~{basename(sv_vcf_original)}
  if [[ ~{basename(sv_vcf_original)} == *.vcf ]];
  then
    bgzip ~{sv_vcf_original}
    tabix ~{sv_vcf_original}.gz
  else
    tabix ~{basename(sv_vcf_original)}
  fi

  comm -13 \
    <(bcftools query -f '%ID\t%MATE_ID\n' ~{basename(sv_svpack_filtered)} | grep -v '\.' | cut -f1 | sort) \
    <(bcftools query -f '%ID\t%MATE_ID\n' ~{basename(sv_svpack_filtered)} | grep -v '\.' | cut -f2 | sort) \
    > missing_mate.txt

  # Extract the missing BNDs
  bcftools view \
    -i 'ID=@missing_mate.txt' \
    ~{basename(sv_vcf_original)} |\
      bcftools sort -Oz -o missing_mate.vcf.gz

  tabix missing_mate.vcf.gz

  # Add the missing BNDs back to the filtered VCF
  bcftools concat -a \
    ~{basename(sv_svpack_filtered)} \
    missing_mate.vcf.gz |\
      bcftools sort -Oz -o tmp.vcf.gz

  mv tmp.vcf.gz ~{basename(sv_svpack_filtered)}
  rm -f ~{basename(sv_svpack_filtered)}.tbi

  tabix -p vcf ~{basename(sv_svpack_filtered)}
  >>>

  output {
    File output_vcf = basename(sv_svpack_filtered)
    File output_vcf_index = basename(sv_svpack_filtered) + ".tbi"
  }

  runtime {
    docker: "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    cpu: 4
    memory: "16 GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

# Use bedtools to split contigs
task split_contigs {
  input {
    File ref_fasta_index
    Int chunk_size = 75000000
    Int threads
  }

  Float file_size = ceil(size(ref_fasta_index, "GB") + 10)

  command <<<
  set -euxo pipefail

  bedtools --version

  echo "Splitting contigs for ~{ref_fasta_index}"
  bedtools makewindows -g ~{ref_fasta_index} -w ~{chunk_size} > contigs.bed
  grep -v -E "random|chrUn|chrM|chrEBV" contigs.bed > noalt.bed
  # Split the contig bed files into one file for each line
  split -l 1 noalt.bed contigs_split.
  # Add .bed to all the contigs_split file
  for file in $(ls contigs_split.*); do mv $file $file.bed; done
  >>>

  output {
    Array[File] contigs = glob("contigs_split.*.bed")
  }

  runtime {
    docker: "quay.io/biocontainers/bedtools:2.31.0--hf5e1c6e_2"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

task mosdepth {
  input {
    String pname
    File bam
    File bam_index
    Int threads
  }

  Float file_size = ceil(size(bam, "GB") * 2)

  command <<<
  set -euxo pipefail

  echo "Running mosdepth for ~{bam}"

  mosdepth --version

  # Run this in case the index doesn't localize properly.
  # Miniwdl v1.1 seems ok, but Cromwell with dev spec will need this
  # to be uncommented otherwise it fails saying it index is not found.
  # ln -s ~{bam_index} ~{bam}.bai
  mosdepth -t ~{threads} \
    --by 500 \
    --no-per-base \
    --use-median \
     ~{pname} \
     ~{bam}
  >>>

  output {
    File output_bed = pname + ".regions.bed.gz"
    File output_bed_index = pname + ".regions.bed.gz.csi"
    File output_summary = pname + ".mosdepth.summary.txt"
  }

  runtime {
    docker: "quay.io/biocontainers/mosdepth:0.3.4--hd299d5a_0"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

task cpg_pileup {
  input {
    File bam
    File bam_index
    File reference
    File reference_index
    String pname
    Int threads
    Int mincov
  }

  String output_prefix = pname + ".cpg"
  Float file_size = ceil(size(bam, "GB") * 2)

  command <<<
  set -euxo pipefail

  echo "Running cpg pileup for ~{pname}"
  
  aligned_bam_to_cpg_scores --version

  aligned_bam_to_cpg_scores \
    --threads ~{threads} \
    --bam ~{bam} \
    --ref ~{reference} \
    --output-prefix ~{output_prefix} \
    --min-mapq 1 \
    --min-coverage ~{mincov} \
    --model "$PILEUP_MODEL_DIR"/pileup_calling_model.v1.tflite

  # gzip all bed
  gzip ~{output_prefix}.combined.bed
  # If hap1 and hap2 beds are present, gzip them
  if [ -f ~{output_prefix}.hap1.bed ]; then
    gzip ~{output_prefix}.hap1.bed
  fi
  if [ -f ~{output_prefix}.hap2.bed ]; then
    gzip ~{output_prefix}.hap2.bed
  fi
  >>>

  output {
    Array[File] pileup_beds = glob("~{output_prefix}.*.bed.gz")
		Array[File] pileup_bigwigs = glob("~{output_prefix}.*.bw")
    File pileup_combined_beds = output_prefix + ".combined.bed.gz"
    File? pileup_hap1_bed = output_prefix + ".hap1.bed.gz"
    File? pileup_hap2_bed = output_prefix + ".hap2.bed.gz"
  }

  runtime {
    docker: "quay.io/pacbio/pb-cpg-tools@sha256:b95ff1c53bb16e53b8c24f0feaf625a4663973d80862518578437f44385f509b"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

task strip_kinetics {
  input {
    File bam
    Int threads
  }

  Float file_size = ceil(size(bam, "GB") * 2)

  command <<<
    set -euxo pipefail

    samtools --version
    
    samtools view \
      -@ ~{threads} \
      -bh \
      ~{bam} \
      -x fi,fp,ri,rp > ~{sub(basename(bam), "\\.bam", ".nokinetics.bam")}
  >>>

  output {
    File stripped_bam = sub(basename(bam), "\\.bam", ".nokinetics.bam")
  }

  runtime {
    docker: "quay.io/biocontainers/samtools:1.17--hd87286a_1"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

task seqkit_bamstats {
  input {
    File bam
    File bam_index
    Int threads
  }

  Float file_size = ceil(size(bam, "GB") * 2)

  command <<<
  set -euxo pipefail

  seqkit version

  seqkit bam -j ~{threads} -s ~{bam} -Z 2> ~{sub(basename(bam), "\\.bam", ".overall.stats.tsv")}
  echo -e '
  Dump:
    Tsv: "~{sub(basename(bam), "\\.bam", ".alignment.stats.tsv")}"
    Fields: ["Read", "Ref", "Pos", "EndPos", "MapQual", "Acc", "Match", "Mismatch", "Ins", "Del", "AlnLen", "ReadLen", "RefLen", "RefAln", "RefCov", "ReadAln", "ReadCov", "Strand", "MeanQual", "LeftClip", "RightClip", "Flags", "IsSec", "  IsSup", "RightSoftClip", "LeftHardClip", "RightHardClip"]
  Sink: True
  ' > seqkit.yaml
  seqkit bam -j ~{threads} -T '{Yaml: "seqkit.yaml"}' ~{bam}
  gzip ~{sub(basename(bam), "\\.bam", ".alignment.stats.tsv")}
  gzip ~{sub(basename(bam), "\\.bam", ".overall.stats.tsv")}
  >>>

  output {
    File seqkit_bam_stats = sub(basename(bam), "\\.bam", ".overall.stats.tsv.gz")
    File seqkit_alignment_stats = sub(basename(bam), "\\.bam", ".alignment.stats.tsv.gz")
  }

  runtime {
    docker: "quay.io/biocontainers/seqkit:2.5.1--h9ee0642_0"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

task summarize_seqkit_alignment {
  input {
    File seqkit_alignment_stats
    Int threads
  }

  Float file_size = ceil(size(seqkit_alignment_stats, "GB") + 10)

  command <<<
    set -euo pipefail

    csvtk version

    echo -e "mean\tmedian\tn50\tsum" > ~{sub(basename(seqkit_alignment_stats), "\\.tsv.gz", ".alignment.summary.tsv")}
    
    # Extract and sort read lengths
    csvtk cut -j ~{threads} -t -fRead,ReadLen ~{seqkit_alignment_stats} | \
      tail -n+2 |
      sort -u -S 4G | \
      cut -f2 | \
      sort -n -S 4G > sorted_lengths.txt
    
    # Calculate basic stats
    total_reads=$(wc -l < sorted_lengths.txt)
    sum=$(awk '{sum += $1} END {print sum}' sorted_lengths.txt)
    mean=$(awk -v sum="$sum" -v count="$total_reads" 'BEGIN {print sum/count}')
    
    # Calculate median
    if [ $((total_reads % 2)) -eq 0 ]; then
      # Even number of reads
      mid1=$((total_reads / 2))
      mid2=$((mid1 + 1))
      median=$(sed -n "${mid1}p;${mid2}p" sorted_lengths.txt | awk '{sum += $1} END {print sum/2}')
    else
      # Odd number of reads
      mid=$(((total_reads + 1) / 2))
      median=$(sed -n "${mid}p" sorted_lengths.txt)
    fi
    
    # Calculate N50
    half_sum=$((sum / 2))
    cumsum=0
    n50=0
    while IFS= read -r length; do
      cumsum=$((cumsum + length))
      if [ $cumsum -ge $half_sum ] && [ $n50 -eq 0 ]; then
        n50=$length
        break
      fi
    done < sorted_lengths.txt
    
    echo -e "${mean}\t${median}\t${n50}\t${sum}" >> ~{sub(basename(seqkit_alignment_stats), "\\.tsv.gz", ".alignment.summary.tsv")}
    
    rm sorted_lengths.txt
  >>>

  output {
    File output_summary = sub(basename(seqkit_alignment_stats), "\\.tsv.gz", ".alignment.summary.tsv")
  }

  runtime {
    docker: "quay.io/biocontainers/csvtk:0.27.2--h9ee0642_0"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

task mutationalpattern {
  input {
    File vcf
    String pname
    Float max_delta
    Int threads
  }

  Float file_size = ceil(size(vcf, "GB") + 10)

  command <<<
  set -euxo pipefail

  Rscript --vanilla /app/mutational_pattern.R \
    ~{vcf} \
    ~{pname} \
    ~{max_delta}
  >>>

  output {
    File mutsig_tsv = pname + ".mut_sigs.tsv"
    File recon_sig = pname + ".reconstructed_sigs.tsv"
    File occurences_tsv = pname + ".type_occurences.tsv"
    File mutsig_bootstrap = pname + ".mut_sigs_bootstrapped.tsv"
    File mut_profile = pname + ".mutation_profile.pdf"
    Array[File] mutsig_output = glob("*.tsv")
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
