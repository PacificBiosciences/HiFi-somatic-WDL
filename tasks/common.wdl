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
  >>>

  output {
    Array[File] pileup_beds = glob("~{output_prefix}.*.bed")
		Array[File] pileup_bigwigs = glob("~{output_prefix}.*.bw")
    File pileup_combined_beds = output_prefix + ".combined.bed"
    File? pileup_hap1_bed = output_prefix + ".hap1.bed"
    File? pileup_hap2_bed = output_prefix + ".hap2.bed"
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
  csvtk cut -j ~{threads} -t -fRead,ReadLen ~{seqkit_alignment_stats} | csvtk uniq -t | csvtk sort -t -k ReadLen:n | awk '{
      sum += $2;         # Sum the values
      values[NR] = $2;  # Store the values in an array
  } 
  END {
      n = length(values);   # Get the total number of values
      # Note that below assume "values" are sorted (using csvtk before input)

      # Calculate and print mean
      mean = sum / n;

      # Calculate and print median
      if (n % 2 == 0) {
          median = (values[n/2] + values[n/2 + 1]) / 2;
      } else {
          median = values[(n + 1) / 2];
      }

      # Calculate N50
      half_sum = sum / 2;
      cumulative_sum = 0;
      for (i = 1; i <= n; i++) {
          cumulative_sum += values[i];
          if (cumulative_sum >= half_sum) {
              n50=values[i];
              break;
          }
      }

      # Print all the values
      print mean, median, n50, sum
  }' OFS=$'\t' >> ~{sub(basename(seqkit_alignment_stats), "\\.tsv.gz", ".alignment.summary.tsv")}
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
      docker: "quay.io/pacbio/somatic_r_tools@sha256:c5e1820b3329e13c092149015d16764e8250755c85a72838dc6ce201b2fa6876"
      cpu: threads
      memory: "~{threads * 4} GB"
      disk: file_size + " GB"
      maxRetries: 2
      preemptible: 1
  }
}
