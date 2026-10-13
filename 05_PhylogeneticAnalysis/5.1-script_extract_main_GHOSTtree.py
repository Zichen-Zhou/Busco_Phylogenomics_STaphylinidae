# Usage:
# python3 5.1-script_extract_main_GHOSTtree.py <input_treefile> <output_treefile>

import sys
import re

def clean_tip_labels(tree_str):
    # Remove anything inside square brackets (e.g., [0.1/0.2/0.3/0.4])
    return re.sub(r'\[.*?\]', '', tree_str)

def extract_first_tree(treefile_path, output_path):
    with open(treefile_path) as f:
        tree_data = f.read().strip()

    # Trees are usually separated by semicolons
    trees = [t for t in tree_data.split(';') if t.strip()]
    if not trees:
        raise ValueError("No valid trees found in the file.")

    # Extract the first tree and clean tip labels
    first_tree = trees[0].strip() + ';'
    cleaned_tree = clean_tip_labels(first_tree)

    with open(output_path, 'w') as out:
        out.write(cleaned_tree + '\n')
    print(f"[INFO] First tree saved to: {output_path}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 5.1-script_extract_main_GHOSTtree.py <input_treefile> <output_treefile>")
        sys.exit(1)

    input_tree = sys.argv[1]
    output_tree = sys.argv[2]
    extract_first_tree(input_tree, output_tree)
