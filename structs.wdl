version 1.1

struct Patient {
    String patient_names
    Array[File] normal_bams
    Array[File] tumor_bams

    String? sex
}

struct Cohort {
    Array[Patient] patients
}