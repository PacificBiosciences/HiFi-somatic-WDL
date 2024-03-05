#!/bin/bash

# This script takes a ClairS VCF file and outputs a new VCF file with a normal sample field:
# 1. All genotype fields are set to 0/0
# 2. The AD field is set to NAD
# 3. The AF field is set to NAF
# 4. The DP field is set to NDP
# 5. All other fields are set to "."
# The script assumes that the input VCF file is compressed with bgzip and indexed with tabix.
# The script also assumes that the input VCF file has the following fields:

set -euxo pipefail

file=$1
tumorname=$2
normalname=$3
ofilename=$4

gunzip -c ${file} | awk 'BEGIN {FS=OFS="\t"} {
    if ($1 ~ /^#/) {  # Skip lines starting with #
        print;
        next;
    }
    
    # Split the format and data fields
    split($9, formats, ":");
    split($10, values, ":");
    split($10, newValues, ":");

    # Initialize the new values array
    newValues[1] = "0/0";  # Set GT to 0/0

    for (i = 2; i <= length(formats); i++) {
        if (formats[i] == "AD") newValues[i] = values[i+3];  # Set AD to NAD
        else if (formats[i] == "AF") newValues[i] = values[i+2];  # Set AF to NAF
        else if (formats[i] == "DP") newValues[i] = values[i+4];  # Set DP to NDP
        else newValues[i] = ".";  # Dummy value for other fields
    }

    # Construct the new data string
    newData = newValues[1];
    for (i = 2; i <= length(newValues); i++) {
        newData = newData ":" newValues[i];
    }

    # Print the line with the new data
    print $0, newData;
}' | sed "s/FORMAT\tSAMPLE/FORMAT\t${tumorname}\t${normalname}/g" | bgzip -c > ${ofilename}

tabix -p vcf ${ofilename}