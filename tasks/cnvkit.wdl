version 1.0

# Call SNV with CNVKit using 10 kbp segmentation window and CBS algorithm
task cnvkit_tumor {
  input {
    File tumor_bam
    File tumor_bam_index
    File normal_bam
    File normal_bam_index
    File ref_fasta
    File ref_fasta_index
    File refFlat
    String pname
    Int threads
  }

  Float file_size = ceil(size(tumor_bam, "GB") * 2 + size(normal_bam, "GB") * 2 + size(ref_fasta, "GB") + 20)

  command <<<
    set -euxo pipefail
    
    echo "Running cnvkit tumor for ~{pname}"

    cnvkit.py version

    cnvkit.py batch \
        ~{tumor_bam} \
        --normal ~{normal_bam} \
        --annotate ~{refFlat} \
        -f ~{ref_fasta} \
        --target-avg-size 10000 -m wgs \
        -p ~{threads} \
        --diagram --scatter \
        --output-dir ~{pname + "_cnvkit"} \
        --segment-method cbs
  >>>

  output {
    Array[File]+ cnvkit_output = glob(pname + "_cnvkit/*")
    File cnvkit_output_cns = pname + "_cnvkit/" + pname + ".tumor.aligned.cns"
  }

  runtime {
    docker: "quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0"
    cpu: threads
    memory: "~{threads * 6} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

# Merge germline variants in tumor and normal for CNVkit
task merge_germline {
  input {
    File tumor_germline_vcf
    File tumor_germline_vcf_tbi
    File normal_germline_vcf
    File normal_germline_vcf_tbi
    String pname
    Int threads
  }

  Float file_size = ceil(size(tumor_germline_vcf, "GB") + size(normal_germline_vcf, "GB") + 10)
  String tumorname = pname + ".tumor"
  String normalname = pname + ".normal"

  command <<<
  set -euxo pipefail

  echo -e "~{normalname}\n~{tumorname}" > sample_names.txt

  bcftools --version

  bcftools merge --force-samples \
    ~{normal_germline_vcf} \
    ~{tumor_germline_vcf} |\
      bcftools reheader -s sample_names.txt -o ~{pname + ".merged_germline.vcf"}

  bgzip ~{pname + ".merged_germline.vcf"}
  
  tabix -p vcf ~{pname + ".merged_germline.vcf.gz"}
  
  bcftools filter ~{pname + ".merged_germline.vcf.gz"} \
    -i '(GT[0]=="0/1" || GT[0]=="1/0") && GT[1]!="./."' \
    -Oz -o ~{pname + ".merged_germline_heterozygous.vcf.gz"}

  tabix -p vcf ~{pname + ".merged_germline_heterozygous.vcf.gz"}

  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%GQ\t%DP\t%AF\t]\n' ~{pname + ".merged_germline_heterozygous.vcf.gz"} | \
  awk -v OFS='\t' '{ 
      n_ref_1 = $7 * (1 - $8); n_alt_1 = $7 * $8;
      ad_1 = int(n_ref_1+0.5) "," int(n_alt_1+0.5);
      
      print $1,$2,ad_1
  }' | bgzip -c > ad_1.txt.gz

  tabix -s1 -b2 -e2 ad_1.txt.gz

  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%GQ\t%DP\t%AF\t]\n' ~{pname + ".merged_germline_heterozygous.vcf.gz"} | \
  awk -v OFS='\t' '{
      n_ref_2 = $11 * (1 - $12); n_alt_2 = $11 * $12;
      ad_2 = int(n_ref_2+0.5) "," int(n_alt_2+0.5);
      
      print $1,$2,ad_2
  }' | bgzip -c > ad_2.txt.gz

  tabix -s1 -b2 -e2 ad_2.txt.gz

  bcftools annotate -s ~{normalname} \
    -a ad_1.txt.gz \
    -c CHROM,POS,FORMAT/AD \
    ~{pname + ".merged_germline_heterozygous.vcf.gz"} |\
      bcftools annotate -s ~{tumorname} \
        -a ad_2.txt.gz \
        -c CHROM,POS,FORMAT/AD \
        -Oz -o ~{pname + ".merged_germline_heterozygous.withAD.vcf.gz"}
  
  tabix -p vcf ~{pname + ".merged_germline_heterozygous.withAD.vcf.gz"}

  >>>

  output {
    File merged_germline_vcf = pname + ".merged_germline.vcf.gz"
    File merged_germline_vcf_tbi = pname + ".merged_germline.vcf.gz.tbi"
    File merged_germline_heterozygous_vcf = pname + ".merged_germline_heterozygous.withAD.vcf.gz"
    File merged_germline_heterozygous_vcf_tbi = pname + ".merged_germline_heterozygous.withAD.vcf.gz.tbi"
  }

  runtime {
    docker: "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    cpu: threads
    memory: "~{threads * 2} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

# Recall with purity and ploidy + BAF from ClairS to obtain
# major and minor copy number for each segment
# Note that purity_ploidy is a headerless TSV file with two columns, 
# first being purity and second being ploidy. 
task cnvkit_recall {
  input {
    File cnvkit_cns
    File merged_germline_heterozygous_vcf
    File merged_germline_heterozygous_vcf_tbi
    Int threads
    File purity_ploidy
    Int min_var_depth = 20
    String pname
  }

  Float file_size = ceil(size(cnvkit_cns, "GB") + size(merged_germline_heterozygous_vcf, "GB") + 10)

  command <<<
  set -euxo pipefail

  # Set integer ploidy
  ploidy=$(cut -f2 ~{purity_ploidy} | awk '{print int($1 + 0.5)}')
  purity=$(cut -f1 ~{purity_ploidy})

  cnvkit.py version

  cnvkit.py call \
    ~{cnvkit_cns} \
    -m clonal \
    --ploidy ${ploidy} \
    --purity ${purity} \
    -o ~{pname + ".CNVKit.with_major_minor_CN.cns"} \
    -v ~{merged_germline_heterozygous_vcf} \
    -i ~{pname + ".tumor"} \
    -n ~{pname + ".normal"} \
    ~{"--min-variant-depth " + min_var_depth}
  >>>

  output {
    File cnvkit_cns_with_major_minor_CN = pname + ".CNVKit.with_major_minor_CN.cns"
  }

  runtime {
    docker: "quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0"
    cpu: threads
    memory: "~{threads * 2} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  } 
}