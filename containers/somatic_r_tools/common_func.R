# Function to count SV type
count_sv <- function(vcf) {
    vcf <- read_tsv(vcf, comment = "#", col_names = FALSE)
    sv_types <- vcf %>%
        mutate(
            sv_type = str_extract(X3, "(?i)SEVERUS_([A-Z]+)\\d+_*", group = 1),
            sv_id = str_extract(X3, "(?i)SEVERUS_([A-Z]+\\d+)_*", group = 1)
        )
    # Group by sv_id
    sv_count <- sv_types %>%
        group_by(sv_id) %>%
        summarise(sv_type = unique(sv_type))
    # Count sv_type
    sv_count <- sv_count %>%
        count(sv_type)
}

make_unique_col <- function(x, y) {
    # Split x and y by semicolon, then
    # paste together and remove duplicates
    x <- str_split(x, ";")
    y <- str_split(y, ";")
    x <- unlist(x)
    y <- unlist(y)
    z <- paste0(x, ":", y)
    z <- unique(z)
    # Remove anything with "False". This is for
    # linking CGC_CANCER_GENE with cancer_type
    z <- z[!str_detect(z, "False")]
    z <- str_c(z, collapse = ";")
    if (z == ":") {
        return("")
    }
    return(z)
}