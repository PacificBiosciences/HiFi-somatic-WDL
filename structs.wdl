version 1.0

struct Patient {
    String patient_names
    Array[File] normal_bams
    Array[File] tumor_bams

    String? sex
}

struct Cohort {
    Array[Patient] patients
}

struct IndexData {
	File data
	File data_index
}