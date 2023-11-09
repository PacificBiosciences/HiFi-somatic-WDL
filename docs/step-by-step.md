# Step by step tutorial to run the workflow on HPC

- [Step by step tutorial to run the workflow on HPC](#step-by-step-tutorial-to-run-the-workflow-on-hpc)
  - [Pre-requisites](#pre-requisites)
  - [Create a working directory](#create-a-working-directory)
  - [Install miniwdl](#install-miniwdl)
  - [Download the resources and references](#download-the-resources-and-references)
  - [Download the annotation and hmftools databases](#download-the-annotation-and-hmftools-databases)
  - [Download demo dataset (COLO829)](#download-demo-dataset-colo829)
  - [Modify the input json file to point to the downloaded files](#modify-the-input-json-file-to-point-to-the-downloaded-files)
  - [Modify the miniwdl config file](#modify-the-miniwdl-config-file)
  - [Run the workflow](#run-the-workflow)
  - [FAQ](#faq)
    - [Can I run with UGE/SGE instead of Slurm?](#can-i-run-with-ugesge-instead-of-slurm)
    - [How do I restart the workflow?](#how-do-i-restart-the-workflow)
    - [Input JSON parameters](#input-json-parameters)

## Pre-requisites

A HPC with Singularity (v3 and above) and Anaconda/Miniconda installed.

## Create a working directory

```bash
mkdir -p /path/to/working/directory
cd /path/to/working/directory
# Clone the repository
git clone https://github.com/PacificBiosciences/wdl-hifisomatic
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
curl -LO 'https://zenodo.org/record/10086866/files/hifisomatic_resources.tar.gz'
# Extract the archive
tar -xzvf hifisomatic_resources.tar.gz
# Check that you have the following files:
├── chr.bed
├── GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta
├── GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.fai
├── GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.dict
├── human_GRCh38_no_alt_analysis_set.trf.bed
├── refFlat.hg38.txt
├── sniffles_all_non_germline.nosamples.vcf.gz
└── sniffles_all_non_germline.nosamples.vcf.gz.tbi
```

## Download the annotation and hmftools databases

Annotation is only carried out if `"hifisomatic.vep_cache"` (For SNV and indels) and/or `"hifisomatic.annotsv_cache"` (For structural variants) are specified. If you do not need to annotate the variants, just remove the lines containing `cache` in the JSON file.

`hmftools` resource is used by Amber, Cobalt and Purple for purity and ploidy estimation. If you do not need to estimate purity and ploidy, just remove the lines containing `"hifisomatic.ensembl_data_dir_tarball"` in the JSON file. Note that if you are running the downsampled demo dataset, you do not need to download the hmftools resource as the demo dataset does not contain enough data to estimate purity and ploidy and will fail in those steps.

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
wget 'https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_refseq_vep_110_GRCh38.tar.gz'

# Download hmftools resource
wget 'https://storage.googleapis.com/hmf-public/HMFtools-Resources/dna_pipeline/v5_33/38/hmf_dna_pipeline_resources.38_v5.33.tar.gz'
```

## Download demo dataset (COLO829)

You may download the full COLO829 and HCC1395 dataset here:

1. COLO829 (60X tumor, 60X normal): <https://downloads.pacbcloud.com/public/revio/2023Q2/COLO829>
2. HCC1395 (60X tumor, 40X normal): <https://downloads.pacbcloud.com/public/revio/2023Q2/HCC1395/>

Or you can download a small demo dataset here that contains the region with truth SVs from Valle-Inclan et. al. 2022. Running this dataset should produce 57/62 (90%) of the truth SVs from the paper. Note that because of the sub-sampling, a duplication in chromosome 8 is called as an inversion. For more details on accuracy please refer to [benchmarking](benchmark.md).

```bash
# Download tumor demo
curl -LO 'https://zenodo.org/record/10086866/files/COLO829.30X.SV_region.bam'
# Download matched normal demo
curl -LO 'https://zenodo.org/record/10086866/files/COLO829BL.30X.SV_region.bam'
```

## Modify the input json file to point to the downloaded files

The input to the workflow is a JSON file describing input parameters and location of BAM files. An example of input json file can be found at [input.example.json](test_data/input.example.json).
  
```bash
cp test_data/input.example.json input.json
```

Modify `input.json` file by using your favourite text editor and changing all
`/path/to` to the path of your working directory. The cohort section should be self-explanatory where you can name the patient and modify the path to the tumor and normal BAM file. If you are only running on a single "patient", just remove the entire "patient2" block (remember to delete the comma after the first patient block!).

## Modify the miniwdl config file

```bash
cp test_data/miniwdl.cfg .
```

The miniwdl config file included in the repository assumes you are using Slurm. Change the line:

`container_backend = slurm_singularity`

to

`container_backend = singularity`  if you do not need to use Slurm. The workflow also assumes Singularity is installed at `/usr/bin/singularity`. If it is not (You can find out the path by typing `which singularity` on your command line), please modify the line:

`exe = ["/usr/bin/singularity"]`

by changing the path to the path of your Singularity installation. Finally, the last line containing `extra_args` allow you to specify any additional Slurm arguments. For example the example json assumes submission to the `compute` partition. If you are using a different partition, you can change `compute` to something else.

## Run the workflow

```bash
# Run the workflow
miniwdl run \
  wdl-hifisomatic/hifisomatic.wdl \
  --input input.json \
  -d /path/to/working/directory/COLO829_demo \
  --cfg miniwdl.cfg
```

Important outputs are listed in the [README](../README.md#important-outputs-from-workflow).

## FAQ

### Can I run with UGE/SGE instead of Slurm?

To run in other scheduler such as SGE, we can make use of the `miniwdl-slurm` package and modify the source code to use `qrsh`

The source to change is `src/miniwdl_slurm/__init__.py` in the `miniwdl-slurm` repo:

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

After modifying the source, simply run `pip install .` (Under your `miniwdl environment`) in the `miniwdl-slurm` directory to install the modified package.

See [example](../utility_scripts/miniwdl-slurm_init_modified.py) modification (line 93 onwards) to use this with a UGE scheduler that makes use of `-l s_vmem` for memory and `qrsh` for submission. Note that the argument extension cannot have any space, so e.g. `-now no` won't work and you must specify `-now` and `-no` separately in the list to extend.

### How do I restart the workflow?

miniwdl is smart enough to not rerun completed tasks by using call-caching. The cache directory is defined in the `miniwdl.cfg` file as:

`dir = "$PWD/miniwdl_call_cache"`

where by default it'll place the call cache directory in the working directory where you submit the miniwdl command. miniwdl will rerun completed tasks if:

1. There's a change in the inputs.
2. There's a change in the tasks' WDL.
3. The folders/files for the completed tasks are deleted.

In addition, all the Singularity container images will by default be downloaded according to this line in the config file:

`image_cache = "$PWD/miniwdl_singularity_cache"`

where again by default it'll place the image cache directory in the working directory where you submit the miniwdl command. We recommend changing this to a shared directory so that you can reuse the downloaded images for other workflow run.

### Input JSON parameters

- If you already have aligned BAM files, specify `"hifisomatic.skip_align": true` in the JSON file. Otherwise, the workflow will align the BAM files using `pbmm2`. If there are multiple BAM files for a sample, they will be merged (post-alignment).
  
- `"hifisomatic.hprc_sniffles_control_vcf"` is used to filter structural variants calls from the workflow for all SV callers. The control VCF is generated by running Sniffles 2.0.7 over 118 HPRC healthy control samples, then merged. This helps to remove common germline SVs from the somatic SVs called by the workflow. If you have your own control VCF, you can use that instead. The workflow uses `truvari bench` to compare the SV calls from the workflow with the control VCF and takes the `fp.vcf.gz` (false positive calls) as the "filtered" VCF. The `truvari` command is as follows:
  
  `truvari bench -p 0 -s 0 -S 0 --sizemax 100000000 --dup-to-ins`

  The command ignores percentage identity (requires sequence-resolved SV), minimum size and maximum size, and consider
  duplications to be the same as insertions.
  
- `threads` for different tasks are suggested based on experience to ensure they do not run out of memory, but may not be optimized for speed and economy.
  
- If `"hifisomatic.strip_kinetics"` is true, the workflow will remove the `fp, ri, rp, ri` tags from the BAM files. This is useful if you want to save disk space and do not need the raw kinetics information.
  
- If `"hifisomatic.call_small_variants"` is true, the workflow will call somatic SNVs and indels with [`ClairS`](https://github.com/HKU-BAL/ClairS). Note that without calling small variants, the workflow will not be able to phase the BAM files.

- The model for `ClairS` can be specified with `"hifisomatic.clairs_platform"`. The default model is for Revio (`hifi_revio`). There is also a Sequel II model available that can be specified with `hifi_sequel2`.

- By default, the pipeline output a set of DMR (differentially methylated region) in the Compendium of Cancer Genes (IntOGen) and filter for at least 50 CpG sites in the DMR region. This can be changed by modifying `hifisomatic.ncg_to_filter`.

- Annotating germline variants can be slow, so by default it is switched off. You can switch on germline variants annotation with `"hifisomatic.annotate_germline" = true`. Note that only germline variants in the normal sample will be annotated.

- `hifisomatic.uselongphase` will use longphase to phase SNV/INDELs and haplotag BAM files with longphase instead of the default hiphase. For samples that have extremely high mutation rate you may want to use longphase as hiphase may run out of memory.
