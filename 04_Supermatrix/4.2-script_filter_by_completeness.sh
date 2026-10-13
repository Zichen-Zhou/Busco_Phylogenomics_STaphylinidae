#!/bin/bash

# Filter aligned sequences by per-sequence completeness.
#
# Usage:
#   4.2-script_filter_by_completeness.sh \
#     -i <input_dir> \
#     -o <output_dir> \
#     -l <locus_list> \
#     -c <completeness_percent> \
#     -r <report_dir>

INPUT_DIR=""
OUT_DIR=""
LIST_FILE=""
THRESHOLD=""
REPORT_DIR=""

while getopts ":i:o:l:c:r:" opt; do
    case "$opt" in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        l) LIST_FILE="$OPTARG" ;;
        c) THRESHOLD="$OPTARG" ;;
        r) REPORT_DIR="$OPTARG" ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

if [[ -z "$INPUT_DIR" || -z "$OUT_DIR" || -z "$LIST_FILE" ||
      -z "$THRESHOLD" || -z "$REPORT_DIR" ]]; then
    echo "Usage: $0 -i <input_dir> -o <output_dir> -l <locus_list> -c <completeness_percent> -r <report_dir>" >&2
    exit 1
fi

if [[ ! -f "$LIST_FILE" ]]; then
    echo "Locus list not found: $LIST_FILE" >&2
    exit 1
fi

if ! awk -v t="$THRESHOLD" 'BEGIN {exit !(t >= 0 && t <= 100)}'; then
    echo "Completeness threshold must be between 0 and 100." >&2
    exit 1
fi

THRESHOLD_DECIMAL=$(awk -v t="$THRESHOLD" 'BEGIN {printf "%.6f", t / 100}')

mkdir -p "$OUT_DIR" "$REPORT_DIR"

LOG_FILE="${REPORT_DIR}/completeness_log.csv"
COUNT_FILE="${REPORT_DIR}/seq_count.csv"

echo "Locus,SeqID,AlignmentLength,NonGapCount,Completeness" > "$LOG_FILE"
echo "Locus,SeqCountBefore,SeqCountAfter" > "$COUNT_FILE"

while IFS= read -r LOCUS; do
    LOCUS="${LOCUS%$'\r'}"
    [[ -z "$LOCUS" ]] && continue

    AA_FILE="${INPUT_DIR}/${LOCUS}.aa.fas"
    NUC_FILE="${INPUT_DIR}/${LOCUS}.nuc.fas"

    if [[ -f "$AA_FILE" ]]; then
        FASTA="$AA_FILE"
        EXT="aa.fas"
    elif [[ -f "$NUC_FILE" ]]; then
        FASTA="$NUC_FILE"
        EXT="nuc.fas"
    else
        echo "[Warning] Alignment not found for $LOCUS" >&2
        continue
    fi

    OUTFILE="${OUT_DIR}/${LOCUS}_len${THRESHOLD}.${EXT}"
    : > "$OUTFILE"

    ALN_LEN=$(
        awk '
            BEGIN {RS=">"; ORS=""}
            NR==2 {
                for (i=2; i<=NF; i++) printf "%s", $i
                print ""
            }
        ' "$FASTA" | tr -d '\n' | wc -c
    )

    if [[ "$ALN_LEN" -eq 0 ]]; then
        echo "[Warning] Zero-length alignment for $LOCUS" >&2
        rm -f "$OUTFILE"
        continue
    fi

    COUNT_BEFORE=$(grep -c '^>' "$FASTA")

    awk \
        -v locus="$LOCUS" \
        -v aln_len="$ALN_LEN" \
        -v logfile="$LOG_FILE" \
        -v outfile="$OUTFILE" \
        -v threshold="$THRESHOLD_DECIMAL" '
        BEGIN {
            seq_id = ""
            buffer_index = 0
            nongap = 0
        }

        function process_sequence(    completeness, i) {
            if (seq_id == "")
                return

            completeness = nongap / aln_len

            printf "%s,%s,%d,%d,%.6f\n", \
                locus, seq_id, aln_len, nongap, completeness >> logfile

            if (completeness >= threshold) {
                for (i=0; i<buffer_index; i++)
                    print buffer[i] >> outfile
            }
        }

        /^>/ {
            process_sequence()

            seq_id = substr($0, 2)
            delete buffer
            buffer[0] = $0
            buffer_index = 1
            nongap = 0
            next
        }

        {
            buffer[buffer_index++] = $0
            seq_line = $0
            gsub(/[-?]/, "", seq_line)
            nongap += length(seq_line)
        }

        END {
            process_sequence()
        }
    ' "$FASTA"

    COUNT_AFTER=$(grep -c '^>' "$OUTFILE" 2>/dev/null || true)
    COUNT_AFTER="${COUNT_AFTER:-0}"

    echo "${LOCUS},${COUNT_BEFORE},${COUNT_AFTER}" >> "$COUNT_FILE"

done < "$LIST_FILE"

echo "[Done] Filtered alignments: $OUT_DIR"
echo "[Done] Completeness log: $LOG_FILE"
echo "[Done] Sequence counts: $COUNT_FILE"
