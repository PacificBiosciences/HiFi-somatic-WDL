# Step by step tutorial to run the workflow on HPC

- [Step by step tutorial to run the workflow on HPC](#step-by-step-tutorial-to-run-the-workflow-on-hpc)
  - [Prerequisites](#prerequisites)
  - [Create a working directory](#create-a-working-directory)
  - [Install miniwdl](#install-miniwdl)
  - [Download the resources and references](#download-the-resources-and-references)
  - [Download the annotation and hmftools databases](#download-the-annotation-and-hmftools-databases)
  - [Download demo dataset (COLO829)](#download-demo-dataset-colo829)
  - [Modify the input JSON file to point to the downloaded files](#modify-the-input-json-file-to-point-to-the-downloaded-files)
  - [Modify the miniwdl config file](#modify-the-miniwdl-config-file)
  - [Run the workflow](#run-the-workflow)
  - [FAQ](#faq)
    - [Can I run with UGE/SGE instead of Slurm via miniwdl?](#can-i-run-with-ugesge-instead-of-slurm-via-miniwdl)
    - [How do I restart the workflow?](#how-do-i-restart-the-workflow)
    - [Is the workflow compatible with Cromwell?](#is-the-workflow-compatible-with-cromwell)
    - [Can I run the workflow if the compute nodes do not have internet access?](#can-i-run-the-workflow-if-the-compute-nodes-do-not-have-internet-access)
    - [I am seeing a lot of false-positives in SNV/INDEL calls. How can I reduce them?](#i-am-seeing-a-lot-of-false-positives-in-snvindel-calls-how-can-i-reduce-them)
    - [Input JSON parameters](#input-json-parameters)

## Prerequisites

An HPC with Singularity (v3 and above) and Anaconda/Miniconda installed.

## Create a working directory

```bash
mkdir -p /path/to/working/directory
cd /path/to/working/directory
# Clone the repository
git clone https://github.com/PacificBiosciences/HiFi-Somatic-WDL
```

## Install miniwdl

```bash
# Install miniwdl with conda
conda create -n miniwdl -c conda-forge -c bioconda miniwdl
conda activate miniwdl
# Install miniwdl-slurm plugin if you are using Slurm
pip install miniwdl-slurm
```

## Download the resources and references

```bash
# Download resource bundle needed for the workflow
curl -LO 'https://zenodo.org/record/14847828/files/hifisomatic_resources.tar.gz'
# Extract the archive
tar -xzvf hifisomatic_resources.tar.gz
# Check that you have the following files:
├── chr.bed
├── GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta
├── GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.fai
├── GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.dict
├── human_GRCh38_no_alt_analysis_set.trf.bed
├── refFlat.hg38.txt
├── severus.jasmine.AN10.AC4.nosample.vcf.gz
├── severus.jasmine.AN10.AC4.nosample.vcf.gz
└── Homo_sapiens.GRCh38.112.chr.reformatted.gff3

# Download 1000G panel of normal from Severus GitHub (optional, please see Severus GitHub for citation)
# this is used for tumor-only calling
wget 'https://github.com/KolmogorovLab/Severus/raw/refs/heads/main/pon/PoN_1000G_hg38_extended.tsv.gz'
```

## Download the annotation and hmftools databases

Annotation is only carried out if `"hifisomatic.vep_cache"` (for SNV and indels) and/or `"hifisomatic.annotsv_cache"` (for structural variants) are specified. If you do not need to annotate the variants, remove the lines containing `cache` in the JSON file.

The `hmftools` resource is used by Amber, Cobalt and Purple for purity and ploidy estimation. If you do not need to estimate purity and ploidy, remove the lines containing `"hifisomatic.ensembl_data_dir_tarball"` in the JSON file. Note that if you are running the downsampled demo dataset, you do not need to download the hmftools resource as the demo dataset does not contain enough data to estimate purity and ploidy and will fail in those steps.

```bash
# Download install script from AnnotSV
wget 'https://github.com/lgmgeo/AnnotSV/raw/master/bin/INSTALL_annotations.sh'
# Download the cache with AnnotSV's script
bash INSTALL_annotations.sh
# The script will create a folder named "AnnotSV_annotations"
# Rename it to just AnnotSV and create a tarball
mv AnnotSV_annotations AnnotSV
tar -czvf annotsv_cache.tar.gz AnnotSV

# Download VEP bundle
wget 'https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_refseq_vep_112_GRCh38.tar.gz'

# Download hmftools resource
wget 'https://storage.googleapis.com/hmf-public/HMFtools-Resources/dna_pipeline/v5_33/38/hmf_dna_pipeline_resources.38_v5.33.tar.gz'
```

## Download demo dataset (COLO829)

Download the full COLO829 and HCC1395 dataset here:

1. COLO829 (60X tumor, 60X normal): <https://downloads.pacbcloud.com/public/revio/2023Q2/COLO829>
2. HCC1395 (60X tumor, 40X normal): <https://downloads.pacbcloud.com/public/revio/2023Q2/HCC1395/>

Or you can download a small demo dataset here that contains the region with truth SVs from Valle-Inclan et. al. 2022.

```bash
# Download tumor demo
curl -LO 'https://zenodo.org/record/10404249/files/COLO829.30X.SV_region.bam'
# Download matched normal demo
curl -LO 'https://zenodo.org/record/10404249/files/COLO829BL.30X.SV_region.bam'
```

## Modify the input JSON file to point to the downloaded files

The input to the workflow is a JSON file describing input parameters and the location of BAM files. An example of an input JSON file can be found at [input.example.json](example_configs/input.example.json). For tumor-only workflow, use [input.example.TO.json](example_configs/input.example.TO.json). The main difference in the tumor-only `json` is the absence of the normal bams and the addition of `severus_pon_tsv` for PON-based filtering in Severus.
  
```bash
cp example_configs/input.example.json input.json
```

Modify the `input.json` file by using your favourite text editor and changing all
`/path/to` to the path of your working directory. The cohort section should be self-explanatory where you can name the patient and modify the path to the tumor and normal BAM file. If you are only running on a single "patient", remove the entire "patient2" block. Remember to delete the comma after the first patient block!

## Modify the miniwdl config file

```bash
cp example_configs/miniwdl.cfg .
```

The miniwdl config file included in the repository assumes you are using Slurm. Change the line:

`container_backend = slurm_singularity`

to

`container_backend = singularity` if you do not need to use Slurm. The workflow also assumes Singularity is installed at `/usr/bin/singularity`. If it is not installed there (find the path by typing `which singularity` on your command line), modify the line:

`exe = ["/usr/bin/singularity"]`

by changing the path to the path of your Singularity installation. Finally, the last line containing `extra_args` lets you specify any additional Slurm arguments. For example, the example JSON assumes submission to the `compute` partition. If you are using a different partition, you can change `compute` to something else.

## Run the workflow

```bash
# Run the tumor/normal workflow
miniwdl run \
  hifisomatic.wdl \
  --input input.json \
  -d /path/to/working/directory/COLO829_demo \
  --cfg miniwdl.cfg

# Run the tumor-only workflow
miniwdl run \
  hifisomatic_tumor_only.wdl \
  --input input.json \
  -d /path/to/working/directory/COLO829_demo \
  --cfg miniwdl.cfg
```

Important outputs are listed in the [README](../README.md#important-outputs-from-workflow).

## FAQ

### Can I run with UGE/SGE instead of Slurm via miniwdl?

To run in other scheduler such as SGE, use the `miniwdl-slurm` package and modify the source code to use `qrsh`

The source file to change is `src/miniwdl_slurm/__init__.py` in the `miniwdl-slurm` repo:

<details>
<summary>Click to expand</summary>

    git clone https://github.com/miniwdl-ext/miniwdl-slurm
    # Edit the following block:
    srun_args = [
        "srun",
        "--job-name", self.run_id,
    ]

    partition = self.runtime_values.get("slurm_partition", None)
    if partition is not None:
        srun_args.extend(["--partition", partition])

    cpu = self.runtime_values.get("cpu", None)
    if cpu is not None:
        srun_args.extend(["--cpus-per-task", str(cpu)])

    memory = self.runtime_values.get("memory_reservation", None)
    if memory is not None:
        # Round to the nearest megabyte.
        srun_args.extend(["--mem", f"{round(memory / (1024 ** 2))}M"])

    gpu = self.runtime_values.get("gpu", None)
    if gpu:
        gpuCount = self.runtime_values.get("gpuCount", 1)
        srun_args.extend(["--gres", f"gpu:{gpuCount}"])

    time_minutes = self.runtime_values.get("time_minutes", None)
    if time_minutes is not None:
        srun_args.extend(["--time", str(time_minutes)])

    slurm_constraint = self.runtime_values.get("slurm_constraint", None)
    if slurm_constraint is not None:
        srun_args.extend(["--constraint", slurm_constraint])

    if self.cfg.has_section("slurm"):
        extra_args = self.cfg.get("slurm", "extra_args")
        if extra_args is not None:
            srun_args.extend(shlex.split(extra_args))
    return srun_args
</details>

After modifying the source file, run `pip install .` (located under your `miniwdl environment`) in the `miniwdl-slurm` directory to install the modified package.

See [example](../example_configs/miniwdl-slurm_init_modified.py) modification (line 93 onwards) to use this with a UGE scheduler that makes use of `-l s_vmem` for memory and `qrsh` for submission. Note that the argument extension cannot have any spaces, so `-now no` will NOT work and you must specify `-now` and `-no` separately in the list to extend.

### How do I restart the workflow?

miniwdl is smart enough to not rerun completed tasks by using call-caching. The cache directory is defined in the `miniwdl.cfg` file as:

`dir = "$PWD/miniwdl_call_cache"`

where by default it will place the call cache directory in the working directory where you submit the miniwdl command. miniwdl will rerun completed tasks if:

1. There is a change in the inputs.
2. There is a change in the task's WDL.
3. The folders/files for the completed tasks are deleted.
4. There is a change in the container image version.

In addition, all the Singularity container images will by default be downloaded according to this line in the config file:

`image_cache = "$PWD/miniwdl_singularity_cache"`

where again by default it will place the image cache directory in the working directory where you submit the miniwdl command. We recommend changing this to a shared directory so that you can reuse the downloaded images for other workflow runs.

### Is the workflow compatible with Cromwell?

We have tested running the workflow using Cromwell version 86. An example Cromwell config file is provided below for Slurm + Singularity (Using Cromwell built-in HSQLDB file-based database). 
There are instructions on using Cromwell on different backends (e.g. LSF, SGE) at the [Cromwell documentation](https://cromwell.readthedocs.io/en/stable/backends/HPC/).

<details>
<summary>cromwell.slurm.conf</summary>

    include required(classpath("application"))


    backend {
    default = slurm

    filesystems {
        local {
            localization: [
                "hard-link", "soft-link", "cached-copy", "copy"
            ]
        }
    }

    database {
        profile = "slick.jdbc.HsqldbProfile$"
        db {
            driver = "org.hsqldb.jdbcDriver"
            url = """
            jdbc:hsqldb:file:cromwell-executions/cromwell-db/cromwell-db;
            shutdown=false;
            hsqldb.default_table_type=cached;hsqldb.tx=mvcc;
            hsqldb.result_max_memory_rows=10000;
            hsqldb.large_data=true;
            hsqldb.applog=1;
            hsqldb.lob_compressed=true;
            hsqldb.script_format=3
            """
            connectionTimeout = 120000
            numThreads = 1
        }
    }

    providers {
        slurm {
        actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"                                                                                     
        config {
            runtime-attributes = """
            Int cpu = 2
            Int memory_gb = 8
            String? docker
            """

            submit = """
                sbatch \
                --wait \
                -J ${job_name} \
                -D ${cwd} \
                -o ${out} \
                -e ${err} \
                --ntasks=1 \
                ${"--cpus-per-task=" + cpu} \
                --mem=${memory_gb}G \
                --wrap "/bin/bash ${script}"
            """

            submit-docker = """
                # Make sure the SINGULARITY_CACHEDIR variable is set. If not use a default
                # based on the users home.
                if [ -z $SINGULARITY_CACHEDIR ]; 
                    then CACHE_DIR=$HOME/.singularity/cache
                    else CACHE_DIR=$SINGULARITY_CACHEDIR
                fi
                # Make sure cache dir exists so lock file can be created by flock
                mkdir -p $CACHE_DIR  
                LOCK_FILE=$CACHE_DIR/singularity_pull_flock
                # Create an exclusive filelock with flock. --verbose is useful for 
                # for debugging, as is the echo command. These show up in `stdout.submit`.
                flock --exclusive --timeout 900 $LOCK_FILE \
                singularity exec --containall docker://${docker} \
                echo "successfully pulled ${docker}!"

                # Submit the script to SLURM
                sbatch \
                --wait \
                -J ${job_name} \
                -D ${cwd} \
                -o ${cwd}/execution/stdout \
                -e ${cwd}/execution/stderr \
                --ntasks=1 \
                ${"--cpus-per-task=" + cpu} \
            --mem=${memory_gb}G \
                --wrap "singularity exec --containall --bind ${cwd}:${docker_cwd} docker://${docker} ${job_shell} ${docker_script}"
            """

            kill = "scancel ${job_id}"
            check-alive = "squeue -j ${job_id}"
            job-id-regex = "Submitted batch job (\\d+).*"
        }
        }
    }
    }
</details>

### Can I run the workflow if the compute nodes do not have internet access?

The workflow will download the Singularity container images from the internet by default. If the compute nodes do not have internet access, 
you can download the Singularity container images on a machine with internet access and then transfer the images to the compute nodes. Then, 
change the `image_cache` option in the miniwdl config file to point to the directory where the Singularity container images are stored. 
Note that miniwdl stores the image with names that replace special `/` and `:` with `_`. 
E.g.: `quay.io/biocontainers/pbmm2:1.12.0--h9ee0642_` will become `docker___quay.io_biocontainers_pbmm2_1.12.0--h9ee0642_0.sif`.

### I am seeing a lot of false-positives in SNV/INDEL calls. How can I reduce them?

You may be using ClairS for SNV/INDEL calling. ClairS was optimized with higher coverage in mind (60X/30X) and the default QUAL model may deviate from what the pipeline filters at (SNV QUAL of 2 and INDEL QUAL of 11). You may want to adjust the QUAL filter value using `hifisomatic.clairs_snv_qual` and `hifisomatic.clairs_indel_qual` in the JSON file.

### Input JSON parameters

- `"hifisomatic.use_deepsomatic"` switches to DeepSomatic for somatic SNV/INDEL calling. DeepSomatic can be computationally expensive - the default uses 64 CPUs and takes ~14-18 hours for a 60X/30X tumor/normal dataset. By default this is set to `true`. Set it to `false` to use ClairS (faster) instead.

- If you already have aligned BAM files, specify `"hifisomatic.skip_align": true` in the JSON file (Note: The pipeline was tested only on `pbmm2`-aligned BAMs). Otherwise, the workflow will align the BAM files using `pbmm2`. If there are multiple BAM files for a sample, they will be merged (post-alignment).
  
- `"hifisomatic.control_vcf"` is used to filter structural variants calls from the workflow for all SV callers. The default control VCF provided in the resource bundle is generated by running Sniffles 2.0.7 over 118 HPRC healthy control samples, then merged. This helps to remove common germline SVs from the somatic SVs called by the workflow. If you have your own control VCF, you can use that instead. The workflow uses `truvari bench` to compare the SV calls from the workflow with the control VCF and takes the `fp.vcf.gz` (false positive calls) as the "filtered" VCF. The pipeline uses `svpack` to filter the VCF.
  
- `threads` for different tasks are suggested based on experience to ensure they do not run out of memory, but may not be optimized for speed and economy. In particular, DeepSomatic can sometimes overload the CPU with more threads than provided. Modify the `threads` parameter to suit your needs.
  
- If `"hifisomatic.strip_kinetics"` is true, the workflow will remove the `fp, ri, rp, ri` tags from the BAM files using `samtools`. This is useful if you want to save disk space and do not need the raw kinetics information.
  
- `"hifisomatic.call_small_variants"` can be used to switch off small variants (SNV/INDEL) calling. By default this is set to `true`. Set it to `false` to skip small variants calling.

- The model for `ClairS` can be specified using `"hifisomatic.clairs_platform"`. The default model is for Revio (`hifi_revio`). There is also a Sequel II model available that can be specified using `hifi_sequel2`. DeepSomatic was trained using Revio data and does not have a Sequel II model.

- By default, the pipeline outputs a set of DMR (differentially methylated region) in the Compendium of Cancer Genes (IntOGen) and filters for at least 50 CpG sites in the DMR region. Change this by modifying `hifisomatic.ncg_to_filter`.

- Annotating germline variants can be slow, so by default it is switched off. You can switch on germline variants annotation using `"hifisomatic.annotate_germline" = true`. Note that only germline variants in the normal sample will be annotated.

- `hifisomatic.uselongphase` will use longphase to phase SNV/INDELs and haplotag BAM files with longphase instead of the default hiphase. This is often unnecessary.

- (Experimental) `hifisomatic.run_savana` will run the SAVANA SV caller in addition to the other SV callers. This is experimental and may not be as well-optimized as the other SV callers.
