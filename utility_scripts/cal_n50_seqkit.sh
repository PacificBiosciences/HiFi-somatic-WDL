#!/bin/bash

# This script calculate N50 for the output per alignment TSV file from seqkit
# output folder per_alignment_(tumor|normal)_stats folder
# Requires csvtk

set -euo pipefail

input_file=$1

echo -e "mean\tmedian\tn50\tsum" 
csvtk uniq -t -f1,12 "${input_file}" | csvtk cut -t -f1,12 | awk '{
    sum += $2;         # Sum the values
    values[NR] = $2;  # Store the values in an array
} 
END {
    n = length(values);   # Get the total number of values
    asort(values, sorted); # Sort the values in ascending order

    # Calculate and print mean
    mean = sum / n;

    # Calculate and print median
    if (n % 2 == 0) {
        median = (sorted[n/2] + sorted[n/2 + 1]) / 2;
    } else {
        median = sorted[(n + 1) / 2];
    }

    # Calculate N50
    half_sum = sum / 2;
    cumulative_sum = 0;
    for (i = 1; i <= n; i++) {
        cumulative_sum += sorted[i];
        if (cumulative_sum >= half_sum) {
            n50=sorted[i];
            break;
        }
    }

    # Print all the values
    print mean, median, n50, sum
}' OFS=$'\t'