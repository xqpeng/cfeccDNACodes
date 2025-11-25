#!/bin/bash
set -euo pipefail

# 参数检查
if [ $# -lt 3 ]; then
    cat << EOF
Uasage: map.sh <input_dir> <output_dir> <genome_index_dir>

Example:
  map.sh /path/to/raw_data /path/to/output /path/to/hg38_index

Attention:
  - fastq filename format: *.R1.fq.gz 和 *.R2.fq.gz
  - bwa and samtools should be installed
EOF
    exit 1
fi

input_dir="$1"
output_dir="$2"
hg_index="$3"

# check bwa and samtools are available
command -v bwa >/dev/null 2>&1 || { echo "error: can not found bwa"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "error: can not found samtools"; exit 1; }

# check input_dir and index file  are available
[ -d "$input_dir" ] || { echo "error: input_dir is not available: $input_dir"; exit 1; }
[ -f "${hg_index}.bwt" ] || { echo "error: BWA index file is not available: ${hg_index}.bwt"; exit 1; }

# build the output_dir 
mkdir -p "$output_dir"
mkdir -p "temp_d"

# process files
for read1_file in "$input_dir"/*.R1.fq.gz; do
    # 
    [ -f "$read1_file" ] || continue
    
    name=$(basename "$read1_file" .R1.fq.gz)
    read2_file="$input_dir/$name.R2.fq.gz"
    
    if [ ! -f "$read2_file" ]; then
        echo "warning: skip $name, can not find R2"
        continue
    fi
    
    echo "====== process sample: $name ======"
    
    # outputfile
    final_bam="$output_dir/${name}.sorted.bam"
    
    # if exist the outputfile, skip 
    if [ -f "$final_bam" ]; then
        echo "skip $name, it has been existed"
        continue
    fi
    
    # pipeline for alignment
    bwa mem -M -t 16 "$hg_index" "$read1_file" "$read2_file" | \
    samtools view -@ 16 -bS - | \
    samtools sort -@ 16 -o "$final_bam" -
    
    echo "完成: $name -> $final_bam"
done

# clear temp files
rm -rf temp_d

echo "====== 所有样本处理完成 ======"