# Json for different testing purposes on PB server

- `input.cancer.json` has COLO829 and HCC1395 at all purities and downsampled coverages. This is a large run to stress test the workflow including
  all features. The BAMs are already aligned to save time (`skip_align` set to true).

- `input.regtest.COLO829_full.json` and `input.regtest.HCC1395_full.json` contains testing of just COLO829 and HCC1395 60X/30X tumor/normal dataset.
  This realigns the BAM file (`skip_align` set to false).

- `input.regtest.downsampled.json` contains testing of downsampled dataset. COLO829 BAMs are regions of truth SVs. HCC1395 BAMs are from chr20. These
  are used for the step-by-step tutorial and are in Zenodo. Note that `ensembl_data_dir` are blank to disable Purple here (purity and ploidy
  estimation).

- `input.regtest.downsampled.skip_smallvariants.json` is similar to `input.regtest.downsampled.json` but with `call_smallvariants` set to false
  to skip SNV/INDELs call. This is to test the logic of calling just SVs.