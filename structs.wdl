version 1.0

struct Patient {
    String patient_names
    Array[File] normal_bams
    Array[File] tumor_bams

    String? sex
}

struct TumorOnlyPatient {
    String patient_names
    Array[File] tumor_bams
}

struct TumorOnlyCohort {
    Array[TumorOnlyPatient] patients
}

struct Cohort {
    Array[Patient] patients
}

struct IndexData {
	File data
	File data_index
}