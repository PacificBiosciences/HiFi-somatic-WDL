# Benchmarking COLO829 (structural variants) and HCC1395 (SNV/INDEL)

In order to simulate tumor purity, we downsample aligned BAM file using samtools and mix the normal and tumor.
For example, a 50% purity 60X tumor consists of 30X tumor merged with 30X normal.

### Structural variants recall for COLO829 (Severus)

Truth sets were obtained from Valle-Inclan et. al. 2022 [(Zenodo record)](https://zenodo.org/record/6426985). Benchmark was carried out
using `truvari` with the parameter of `-p 0 -s 0 -S 0 --sizemax 100000000 --dup-to-ins`. The following results are obtained from Severus' VCF file. All 57 truth SVs called for 60X (100% purity) had also been manually checked to ensure Truvari was correctly matching the SVs.

| Caller  | Normal coverage | Tumor coverage | Purity | Truth SV Called | Total SV Called | Recall (Truth / 62) | Precision |
| ------- | --------------- | -------------- | ------ | --------------- | --------------- | ------------------- | --------- |
| Severus | 30              | 30             | 0.5    | 52              | 97              | 0.84                | 0.54      |
| Severus | 30              | 30             | 1      | 57              | 136             | 0.92                | 0.42      |
| Severus | 30              | 60             | 0.25   | 52              | 118             | 0.84                | 0.44      |
| Severus | 30              | 60             | 0.5    | 57              | 151             | 0.92                | 0.38      |
| Severus | 30              | 60             | 1      | 57              | 246             | 0.92                | 0.23      |
| Severus | 60              | 60             | 0.25   | 52              | 91              | 0.84                | 0.57      |
| Severus | 60              | 60             | 0.5    | 57              | 126             | 0.92                | 0.45      |
| Severus | 60              | 60             | 1      | 57              | 204             | 0.92                | 0.28      |


![Alt text](../figures/severus_recall.png)

There are 4 structural variants that Severus failed to called. `truthset_13` and `truthset_42` did not show any evidence based on manual inspection via IGV. `truthset_32` has just 2 variant reads out of 102 reads, which is 2% VAF. `truthset_51` does not have any coverage in the normal. Finally, `truthset_52` is part of a complex structural variants chain with truthset 19, 20, 53 and 54.

### Small variants recall for HCC1395 (ClairS)

Truth sets were obtained from SEQC2 [website](https://sites.google.com/view/seqc2/home/data-analysis) (v1.2.1) and we compared
the VCFs using [`som.py`](https://github.com/Illumina/hap.py/blob/master/doc/sompy.md). Note that this is different from the benchmarking
method in ClairS preprint where truth variants that do not have any supporting reads or have very low coverage in the tumor BAM are removed. We did not remove any truth variants in this benchmark.

| Tumor Coverage | Normal Coverage | Purity | type     | total.truth | total.query | tp    | fp   | fn   | recall  | precision | F1     |
| -------------- | --------------- | ------ | -------- | ----------- | ----------- | ----- | ---- | ---- | ------- | --------- | ------ |
| 30             | 30              | 1      | indels   | 1602        | 1547        | 1074  | 473  | 528  | 0.67041 | 0.69425   | 68.21% |
| 30             | 30              | 1      | SNVs     | 39447       | 36984       | 34210 | 2774 | 5237 | 0.86724 | 0.92499   | 89.52% |
| 30             | 30              | 1      | combined | 41049       | 38531       | 35284 | 3247 | 5765 | 0.85956 | 0.91573   | 88.68% |
| 60             | 30              | 1      | indels   | 1602        | 1475        | 1103  | 372  | 499  | 0.68851 | 0.7478    | 71.69% |
| 60             | 30              | 1      | SNVs     | 39447       | 38830       | 36249 | 2581 | 3198 | 0.91893 | 0.93353   | 92.62% |
| 60             | 30              | 1      | records  | 41049       | 40305       | 37352 | 2953 | 3697 | 0.90994 | 0.92673   | 91.83% |

![Alt text](../figures/SNV_INDEL_VAF.png)
