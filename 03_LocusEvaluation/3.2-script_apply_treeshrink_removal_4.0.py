#!/usr/bin/env python3

# Apply a TreeShrink/manual removal table to nucleotide and protein locus files.
# Usage:
# python3 3.2-script_apply_treeshrink_removal_4.0.py \
#   -l <loci_list> -r <removal.csv> -i <input_dir> -o <output_dir> \
#   -c <sequence_count_log> -n <loci_name_log>

import argparse
import os
import shutil


def read_fasta(filepath):
    sequences = []
    current_header = None
    current_seq = []

    with open(filepath) as handle:
        for line in handle:
            line = line.strip()

            if line.startswith(">"):
                if current_header is not None:
                    sequences.append((current_header, "".join(current_seq)))

                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_header is not None:
            sequences.append((current_header, "".join(current_seq)))

    return sequences


def write_fasta(filepath, sequences):
    with open(filepath, "w") as handle:
        for header, sequence in sequences:
            handle.write(f">{header}\n")
            for i in range(0, len(sequence), 60):
                handle.write(sequence[i:i + 60] + "\n")


def resolve_output_path(output_dir, path):
    if os.path.isabs(path) or os.path.dirname(path):
        return path
    return os.path.join(output_dir, path)


def main():
    parser = argparse.ArgumentParser(
        description="Apply TreeShrink/manual sequence removal to locus FASTA files."
    )
    parser.add_argument("-l", "--loci_list", required=True,
                        help="Locus list file.")
    parser.add_argument("-r", "--removal_csv", required=True,
                        help="CSV file containing Locus and SeqID columns.")
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Input directory containing faa/ and fna/.")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="Output directory.")
    parser.add_argument("-c", "--count_log", required=True,
                        help="Sequence-count log filename or path.")
    parser.add_argument("-n", "--name_log", required=True,
                        help="Output locus-name log filename or path.")
    args = parser.parse_args()

    src_dir = args.input_dir
    dst_dir = args.output_dir
    count_log = resolve_output_path(dst_dir, args.count_log)
    name_log = resolve_output_path(dst_dir, args.name_log)

    before_log = os.path.join(dst_dir, "sequence_number_before.log")
    removed_log = os.path.join(dst_dir, "sequence_removed.log")

    os.makedirs(os.path.join(dst_dir, "faa"), exist_ok=True)
    os.makedirs(os.path.join(dst_dir, "fna"), exist_ok=True)

    with open(args.loci_list) as handle:
        loci_list = [line.strip() for line in handle if line.strip()]

    removal_dict = {}
    with open(args.removal_csv) as handle:
        next(handle, None)
        for line in handle:
            fields = [field for field in line.strip().split(",") if field]
            if len(fields) < 2:
                continue

            locus, seqid = fields[0], fields[1]
            removal_dict.setdefault(locus, set()).add(seqid)

    extension_candidates = {
        "faa": ["faa", "aa.fas"],
        "fna": ["fna", "nuc.fas"],
    }

    with (
        open(count_log, "w") as count_handle,
        open(before_log, "w") as before_handle,
        open(removed_log, "w") as removed_handle,
    ):
        count_handle.write("Locus\tfaa\tfna\n")
        before_handle.write("Locus\tfaa\tfna\n")
        removed_handle.write("Locus\tfaa_removed\tfna_removed\n")

        for locus in loci_list:
            counts_before = {}
            counts_after = {}
            removed_counts = {}

            for sequence_type in ("faa", "fna"):
                src_file = None
                dst_file = None

                for extension in extension_candidates[sequence_type]:
                    candidate = os.path.join(
                        src_dir, sequence_type, f"{locus}.{extension}"
                    )
                    if os.path.exists(candidate):
                        src_file = candidate
                        dst_file = os.path.join(
                            dst_dir, sequence_type, f"{locus}.{extension}"
                        )
                        break

                if src_file is None:
                    print(
                        f"[Warning] Missing source file for {locus} "
                        f"in {sequence_type}/"
                    )
                    counts_before[sequence_type] = 0
                    counts_after[sequence_type] = 0
                    removed_counts[sequence_type] = 0
                    continue

                shutil.copyfile(src_file, dst_file)
                sequences = read_fasta(dst_file)
                counts_before[sequence_type] = len(sequences)

                to_remove = removal_dict.get(locus, set())
                filtered_sequences = [
                    sequence for sequence in sequences
                    if sequence[0] not in to_remove
                ]

                write_fasta(dst_file, filtered_sequences)

                counts_after[sequence_type] = len(filtered_sequences)
                removed_counts[sequence_type] = (
                    counts_before[sequence_type] - counts_after[sequence_type]
                )

                if to_remove:
                    print(
                        f"[INFO] {locus}.{sequence_type}: removed "
                        f"{removed_counts[sequence_type]} sequences "
                        f"({counts_before[sequence_type]} -> "
                        f"{counts_after[sequence_type]})"
                    )

            before_handle.write(
                f"{locus}\t{counts_before['faa']}\t{counts_before['fna']}\n"
            )
            count_handle.write(
                f"{locus}\t{counts_after['faa']}\t{counts_after['fna']}\n"
            )
            removed_handle.write(
                f"{locus}\t{removed_counts['faa']}\t"
                f"{removed_counts['fna']}\n"
            )

    with open(name_log, "w") as handle:
        for locus in loci_list:
            handle.write(f"{locus}\n")

    print(f"Total loci processed: {len(loci_list)}")


if __name__ == "__main__":
    main()
