#!/bin/bash


usage() {
    echo "Usage: $0 -i <input_dir> -o <output_dir>"
    echo "  -i: Input directory containing BAM files"
    echo "  -o: Output directory to save SAM files"
    echo "Example: $0 -i ./bams -o ./sams"
    exit 1
}

if ! command -v samtools &> /dev/null; then
    echo "Error: samtools is not installed or not in PATH."
    exit 1
fi


input_dir=""
output_dir=""


while getopts "i:o:h" opt; do
    case $opt in
        i)
            input_dir="$OPTARG"
            ;;
        o)
            output_dir="$OPTARG"
            ;;
        h)
            usage
            ;;
        *)
            usage
            ;;
    esac
done


if [[ -z "$input_dir" ]] || [[ -z "$output_dir" ]]; then
    echo "Error: Both -i and -o are required."
    usage
fi


if [[ ! -d "$input_dir" ]]; then
    echo "Error: Input directory '$input_dir' does not exist."
    exit 1
fi

mkdir -p "$output_dir"


shopt -s nullglob  
bam_files=("$input_dir"/*.bam)

if [[ ${#bam_files[@]} -eq 0 ]]; then
    echo "No BAM files found in '$input_dir'."
    exit 0
fi

echo "Found ${#bam_files[@]} BAM file(s). Converting to SAM..."

for bam in "${bam_files[@]}"; do
    if [[ ! -f "$bam" ]]; then
        continue
    fi

    basename=$(basename "$bam" .bam)
    sam_file="$output_dir/${basename}.sam"

    echo "Converting: $bam → $sam_file"
    samtools view -h "$bam" > "$sam_file"

    if [[ $? -ne 0 ]]; then
        echo "Error converting $bam"
    else
        echo "Successfully converted $bam"
    fi
done

echo "Batch conversion completed."
