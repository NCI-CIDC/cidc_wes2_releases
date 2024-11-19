#!/bin/bash
# Initialize the report directory and file
# REQUIRES
# /home/pipeline/cidc_wes2/config/tumor-normal-pair-mapping.tsv
# /home/pipeline/cidc_wes2/config/tumor-only-mapping.tsv
# Usage:
# ./wes-schema-mapping.v3.sh

mkdir -p analysis/report
report_file="analysis/report/report.txt"
: > "$report_file"

# Load mapping files into associative arrays
declare -A TN_MAPPING
declare -A TO_MAPPING
declare -A TN_MODIFICATIONS
declare -A TO_MODIFICATIONS

# Load tumor-normal-pair-mapping.tsv
while IFS=$'\t' read -r source target modification; do
    TN_MAPPING["$source"]="$target"
    TN_MODIFICATIONS["$source"]="$modification"
done < /home/pipeline/cidc_wes2/config/tumor-normal-pair-mapping.tsv

# Load tumor-only-mapping.tsv
while IFS=$'\t' read -r source target modification; do
    TO_MAPPING["$source"]="$target"
    TO_MODIFICATIONS["$source"]="$modification"
done < /home/pipeline/cidc_wes2/config/tumor-only-mapping.tsv

# Function to process modifications
process_modification() {
    local src_file="$1"
    local target_file="$2"
    local modification="$3"

    mkdir -p "$(dirname "$target_file")" # Ensure the target directory exists

    if [[ "$modification" == "TOUCH" ]]; then
        touch "$target_file"
        echo "Created dummy file $target_file" | tee -a "$report_file"
    elif [[ "$modification" =~ ^echo ]]; then
        eval "$modification > \"$target_file\""
        if [[ -f "$target_file" ]]; then
            echo "Executed echo command to create $target_file" | tee -a "$report_file"
        else
            echo "Error: Failed to execute echo command for $target_file" | tee -a "$report_file"
        fi
    elif [[ "$modification" == "gunzip vcf.gz file to create .vcf file" ]]; then
        if [[ -f "$src_file" ]]; then
            gunzip -c "$src_file" > "${target_file%.gz}"
            echo "Gunzipped $src_file to ${target_file%.gz}" | tee -a "$report_file"
        else
            echo "Source file $src_file not found for gunzip." | tee -a "$report_file"
        fi
    else
        echo "Unknown modification type: $modification" | tee -a "$report_file"
    fi
}

# Process pairings.csv file
while IFS=',' read -r type run_id normal_id tumor_id; do
    case "$type" in
        "TN")
            echo "Processing Tumor-Normal Pair: Run=$run_id, Normal=$normal_id, Tumor=$tumor_id" | tee -a "$report_file"
            touch maf.vcf.dummy.$run_id.file1.txt
            touch maf.vcf.dummy.$run_id.file2.txt
            touch maf.vcf.dummy.$run_id.file3.txt
            echo 'Dummy tnscope files created'

            for source_path in "${!TN_MAPPING[@]}"; do
                target_path="${TN_MAPPING[$source_path]}"
                modification="${TN_MODIFICATIONS[$source_path]}"

                # Replace placeholders
                source_path="${source_path//\{run\}/$run_id}"
                source_path="${source_path//\{normal_id\}/$normal_id}"
                source_path="${source_path//\{tumor_id\}/$tumor_id}"
                target_path="${target_path//\{run\}/$run_id}"
                target_path="${target_path//\{normal_id\}/$normal_id}"
                target_path="${target_path//\{tumor_id\}/$tumor_id}"

                # Handle file copying or dummy file creation
                if [[ -f "$source_path" ]]; then
                    mkdir -p "$(dirname "$target_path")"
                    cp "$source_path" "$target_path"
                    echo "Copied $source_path to $target_path" | tee -a "$report_file"
                else
                    echo "Source file $source_path not found for mapping." | tee -a "$report_file"
                fi

                # Process modifications
                process_modification "$source_path" "$target_path" "$modification"
            done
            ;;
        "TO")
            echo "Processing Tumor-Only Data: Run=$run_id, Tumor=$tumor_id" | tee -a "$report_file"
            touch maf.vcf.dummy.$run_id.file1.txt
            touch maf.vcf.dummy.$run_id.file2.txt
            touch maf.vcf.dummy.$run_id.file3.txt
            echo 'Dummy tnscope files created'

            for source_path in "${!TO_MAPPING[@]}"; do
                target_path="${TO_MAPPING[$source_path]}"
                modification="${TO_MODIFICATIONS[$source_path]}"

                # Replace placeholders
                source_path="${source_path//\{run\}/$run_id}"
                source_path="${source_path//\{tumor_id\}/$tumor_id}"
                target_path="${target_path//\{run\}/$run_id}"
                target_path="${target_path//\{tumor_id\}/$tumor_id}"

                # Handle file copying or dummy file creation
                if [[ -f "$source_path" ]]; then
                    mkdir -p "$(dirname "$target_path")"
                    cp "$source_path" "$target_path"
                    echo "Copied $source_path to $target_path" | tee -a "$report_file"
                else
                    echo "Source file $source_path not found for mapping." | tee -a "$report_file"
                fi

                # Process modifications
                process_modification "$source_path" "$target_path" "$modification"
            done
            ;;
        *)
            echo "Unknown type: $type" | tee -a "$report_file"
            ;;
    esac
done < /home/pipeline/cidc_wes2/config/pairings.csv

# Compress the report
tar -czf /media/storage/wes2_output/upload/analysis/report.tar.gz -C analysis report
echo "Report generated at analysis/report.tar.gz and printed to screen:" | tee -a "$report_file"
cat "$report_file"

mkdir -p /media/storage/wes2_output/upload/analysis/report/WES_Meta/
echo "Version 2.0 Enhanced Pipeline - Docker Version a0041d20aa0e5" > /media/storage/wes2_output/upload/analysis/report/WES_Meta/02_WES_Run_Version.tsv
echo 'Run Version File Created'