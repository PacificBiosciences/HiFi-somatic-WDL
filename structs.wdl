version 1.0

struct Sample {
    String sample_names
    Array[File] normal_bams
    Array[File] tumor_bams

    String? sex
}

struct TumorOnlySample {
    String sample_names
    Array[File] tumor_bams
}

struct TumorOnlyCohort {
    Array[TumorOnlySample] samples
}

struct Cohort {
    Array[Sample] samples
}

struct IndexData {
	File data
	File data_index
}
