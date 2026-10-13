#!/usr/bin/env bash

set -euo pipefail

# Assemble paired-end RNA-seq reads, reduce redundant contigs, scaffold, and close gaps.
# Run in a directory containing <sample_id>_1.fastq and <sample_id>_2.fastq.
# Usage: bash 1-code_nuclear_assembly_RNASeq.sh -i SAMPLE_ID -s SPECIES -t THREADS -l READ_LENGTH

sample_id=""
species=""
threads=""
read_length=""

usage() {
    echo "Usage: $0 -i SAMPLE_ID -s SPECIES -t THREADS -l READ_LENGTH" >&2
    exit 1
}

while getopts ":i:s:t:l:" opt; do
    case "$opt" in
        i) sample_id="$OPTARG" ;;
        s) species="$OPTARG" ;;
        t) threads="$OPTARG" ;;
        l) read_length="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -n "$sample_id" && -n "$species" && -n "$threads" && -n "$read_length" ]] || usage

work_dir="$(pwd)"
raw_dir="$work_dir/0-rawdata"
assembly_dir="$work_dir/1-assembly"
assembly_out="$assembly_dir/megahit.out"
scaffold_dir="$work_dir/2-scaffolding"
gap_dir="$work_dir/3-gapclosing"
log_file="$work_dir/log.txt"
stats_file="$work_dir/assembly.statistics"

read1="$work_dir/${sample_id}_1.fastq"
read2="$work_dir/${sample_id}_2.fastq"

[[ -s "$read1" ]] || { echo "Error: missing input file: $read1" >&2; exit 1; }
[[ -s "$read2" ]] || { echo "Error: missing input file: $read2" >&2; exit 1; }

for path in "$raw_dir" "$assembly_dir" "$scaffold_dir" "$gap_dir"; do
    [[ ! -e "$path" ]] || { echo "Error: output path already exists: $path" >&2; exit 1; }
done

: > "$log_file"
: > "$stats_file"

# Read processing
mkdir -p "$raw_dir/fastqc"

echo "Running FastQC on raw reads" | tee -a "$log_file"
fastqc "$read1" "$read2" -o "$raw_dir/fastqc"

mv "$read1" "$read2" "$raw_dir/"
read1="$raw_dir/${sample_id}_1.fastq"
read2="$raw_dir/${sample_id}_2.fastq"

adapt1="$raw_dir/${sample_id}_1_adapttrimmed.fq"
adapt2="$raw_dir/${sample_id}_2_adapttrimmed.fq"

fastp \
    -i "$read1" \
    -I "$read2" \
    -o "$adapt1" \
    -O "$adapt2" \
    -Q -L \
    --detect_adapter_for_pe \
    -j "$raw_dir/fastp_report.json" \
    -h "$raw_dir/fastp_report.html"

fastqc "$adapt1" "$adapt2" -o "$raw_dir/fastqc"

clumped1="$raw_dir/1.clumped.fq.gz"
clumped2="$raw_dir/2.clumped.fq.gz"
trimmed1_gz="$raw_dir/1.quatrim.fq.gz"
trimmed2_gz="$raw_dir/2.quatrim.fq.gz"
trimmed1="$raw_dir/1.quatrim.fq"
trimmed2="$raw_dir/2.quatrim.fq"

clumpify.sh \
    in1="$adapt1" \
    in2="$adapt2" \
    out1="$clumped1" \
    out2="$clumped2" \
    -Xmx20g pigz dedupe

bbduk.sh \
    in1="$clumped1" \
    in2="$clumped2" \
    out1="$trimmed1_gz" \
    out2="$trimmed2_gz" \
    pigz ordered \
    qtrim=rl trimq=20 minlen=15 ecco=t maxns=5 \
    trimpolya=10 trimpolyg=10 trimpolyc=10

gzip -d -k "$trimmed1_gz" "$trimmed2_gz"

# Genome assembly
mkdir "$assembly_dir"

echo "Running MEGAHIT" | tee -a "$log_file"
megahit \
    -1 "$trimmed1_gz" \
    -2 "$trimmed2_gz" \
    -o "$assembly_out" \
    --k-list 21,41,61,81,101,121 \
    -t "$threads"

contigs="$assembly_out/final.contigs.fa"
statswrapper.sh in="$contigs" format=6 >> "$stats_file"

# Redundant-contig reduction
redundans.py \
    -v \
    -f "$contigs" \
    -o "$assembly_out/reduced" \
    -t "$threads" \
    --log "$assembly_out/log_redundans.txt" \
    --noscaffolding \
    --nogapclosing \
    --identity 0.7

[[ -s "$assembly_out/log_redundans.txt" ]] && cat "$assembly_out/log_redundans.txt" >> "$log_file"

reduced_contigs="$assembly_out/final.contigs.reduced.fa"
mv "$assembly_out/reduced/scaffolds.reduced.fa" "$reduced_contigs"
rm -rf "$assembly_out/reduced"
rm -f "$assembly_out"/*.fa.fai "$assembly_out/log_redundans.txt" "$assembly_out/reads.list"

# Scaffolding
mkdir "$scaffold_dir"
map_bam="$scaffold_dir/map.bam"

minimap2 \
    -ax sr \
    -t "$threads" \
    "$reduced_contigs" \
    "$trimmed1_gz" \
    "$trimmed2_gz" |
    samtools view -bS -@ "$threads" - |
    samtools sort -@ "$threads" -o "$map_bam" -

samtools index -@ "$threads" "$map_bam"

besst_scaffolds="$scaffold_dir/BESST_output/pass1/Scaffolds_pass1.fa"

run_besst() {
    runBESST \
        -c "$reduced_contigs" \
        -f "$map_bam" \
        -o ./ \
        -orientation fr \
        --iter 10000 \
        "$@"
}

cd "$scaffold_dir"
run_besst || true

if [[ ! -s "$besst_scaffolds" ]]; then
    read -r -p "Enter mean insert size and standard deviation, for example 300 50: " insert_mean insert_sd
    printf "mean_insert_size=%s\ninsert_size_sd=%s\n" "$insert_mean" "$insert_sd" >> parameters.cfg
    run_besst -m "$insert_mean" -s "$insert_sd" || true
fi

if [[ ! -s "$besst_scaffolds" ]]; then
    read -r -p "Enter revised mean insert size and standard deviation: " insert_mean insert_sd
    printf "mean_insert_size=%s\ninsert_size_sd=%s\n" "$insert_mean" "$insert_sd" >> parameters.cfg
    run_besst -m "$insert_mean" -s "$insert_sd" || true
fi

[[ -s "$besst_scaffolds" ]] || { echo "Error: BESST scaffolding failed" >&2; exit 1; }
cd "$work_dir"

# Gap closing
mkdir "$gap_dir"
config_file="$gap_dir/gapcloser.config"
gapclosed="$gap_dir/scaffolds.gapcloser.fa"

cat > "$config_file" <<CONFIG
[LIB]
q1=$trimmed1
q2=$trimmed2
CONFIG

GapCloser \
    -a "$besst_scaffolds" \
    -b "$config_file" \
    -o "$gapclosed" \
    -l "$read_length" \
    -t "$threads"

[[ -s "$gapclosed" ]] || { echo "Error: GapCloser failed" >&2; exit 1; }

final_assembly="$work_dir/scaffolds_${species}.fa"
reformat.sh in="$gapclosed" out="$final_assembly" minlength=1

statswrapper.sh \
    in="$reduced_contigs,$besst_scaffolds,$gapclosed,$final_assembly" \
    >> "$stats_file"

rm -f "$map_bam" "$map_bam.bai"

echo "Final assembly: $final_assembly"
