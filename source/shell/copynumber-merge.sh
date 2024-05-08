#!/bin/bash

for term in GAIN LOSS; do
    ## Extract lines that contain either "GAIN" or "LOSS"
    entries=$(grep $term $1)

    ## Check if any entries are present and combine overlapping features into a singular, all-encompassing entry in the bed file
    if [ -n "$entries" ]; then
        echo "$entries" | mergeBed -c 4 -o distinct > $2_consensus_merged_$term.bed
    ## If no entries found, create empty file
    else
        echo "No lines containing $term were found."
	touch $2_consensus_merged_$term.bed
    fi
done
