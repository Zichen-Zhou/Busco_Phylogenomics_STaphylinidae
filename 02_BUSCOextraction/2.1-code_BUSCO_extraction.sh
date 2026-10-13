# BUSCO locus extraction and filtering
# Run from the directory containing 0-raw_busco/.
# Required software: TransDecoder

# Generate sample list
ls -1 0-raw_busco > species.list
SPECIES_NAME=$(cat species.list)

# Generate locus list
: > temp.list
for SPECIES in $SPECIES_NAME
do
    cd "0-raw_busco/$SPECIES/busco_result_SF" || continue
    ls *.fna *.faa 2>/dev/null | sed -E 's/\.(fna|faa)$//' >> ../../../temp.list
    cd ../../../ || exit
 done
sort -nu temp.list > loci.list
rm temp.list

# Convert multi-line FASTA sequences to single-line format
for SPECIES in $SPECIES_NAME
do
    cd "0-raw_busco/$SPECIES/busco_result_SF" || continue
    for FILE in *.faa *.fna
    do
        [ -e "$FILE" ] || continue
        awk '/^>/{printf("%s\n",$0);next} {printf("%s",$0)} END{printf("\n")}' "$FILE" > temp.fasta
        mv temp.fasta "$FILE"
    done
    cd ../../../ || exit
 done

# Recover the longest ORF and standardize sequence headers
for SPECIES in $SPECIES_NAME
do
    cd "0-raw_busco/$SPECIES/busco_result_SF" || continue

    sed -i -n '1,2p' *.fna *.faa

    for FNA in *.fna
    do
        [ -e "$FNA" ] || continue
        TransDecoder.LongOrfs -t "$FNA" -m 20
        cp "${FNA}.transdecoder_dir/longest_orfs.cds" "$FNA"
        cp "${FNA}.transdecoder_dir/longest_orfs.pep" "${FNA%.fna}.faa"
        rm -rf "${FNA}.transdecoder_dir" "${FNA}.transdecoder_dir.__checkpoints" *.cmds
    done

    sed -i -n '1,2p' *.fna *.faa
    sed -i "1c >$SPECIES" *.fna *.faa
    sed -i 's/\*//g' *.faa

    cd ../../../ || exit
 done

# Merge sequences from the same BUSCO locus across samples
mkdir -p 1-raw_loci/fna 1-raw_loci/faa
LOCI_NAME=$(cat loci.list)
: > 1-raw_loci/log.txt

for LOCI in $LOCI_NAME
do
    : > "1-raw_loci/fna/$LOCI.fna"
    : > "1-raw_loci/faa/$LOCI.faa"

    for SPECIES in $SPECIES_NAME
    do
        FNA_FILE="0-raw_busco/$SPECIES/busco_result_SF/$LOCI.fna"
        FAA_FILE="0-raw_busco/$SPECIES/busco_result_SF/$LOCI.faa"

        if [[ -f "$FNA_FILE" && -f "$FAA_FILE" ]]; then
            cat "$FNA_FILE" >> "1-raw_loci/fna/$LOCI.fna"
            cat "$FAA_FILE" >> "1-raw_loci/faa/$LOCI.faa"
        else
            echo "$LOCI in $SPECIES does not exist" >> 1-raw_loci/log.txt
        fi
    done
 done

# Retain loci represented by at least three samples
mkdir -p 2-loci_filter/fna 2-loci_filter/faa
: > 2-loci_filter/sequence_number.log
: > 2-loci_filter/loci_name_filter.log

for LOCI in $LOCI_NAME
do
    FNA_SEQ_COUNT=$(grep -c '^>' "1-raw_loci/fna/$LOCI.fna")
    FAA_SEQ_COUNT=$(grep -c '^>' "1-raw_loci/faa/$LOCI.faa")

    printf '%s\t%s\t%s\n' "$LOCI" "$FNA_SEQ_COUNT" "$FAA_SEQ_COUNT" >> 2-loci_filter/sequence_number.log

    if [[ "$FNA_SEQ_COUNT" -ge 3 && "$FAA_SEQ_COUNT" -ge 3 ]]; then
        cp "1-raw_loci/fna/$LOCI.fna" 2-loci_filter/fna/
        cp "1-raw_loci/faa/$LOCI.faa" 2-loci_filter/faa/
        echo "$LOCI" >> 2-loci_filter/loci_name_filter.log
    fi
 done
