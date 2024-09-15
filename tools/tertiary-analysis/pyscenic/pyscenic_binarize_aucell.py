import argparse

import pandas as pd
from pyscenic.binarization import binarize

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Binarize AUC matrix")
    parser.add_argument("input_file", help="Input TSV or CSV file")
    parser.add_argument(
        "--threshold-overrides",
        type=str,
        help="Threshold overrides in JSON format",
    )
    parser.add_argument("--seed", type=int, default=None, help="Random seed")
    parser.add_argument(
        "--num-workers", type=int, default=1, help="Number of workers"
    )
    parser.add_argument(
        "--output-prefix", type=str, default="output", help="Output prefix"
    )

    args = parser.parse_args()

    # Read input file
    if args.input_file.endswith(".tsv"):
        auc_mtx = pd.read_csv(args.input_file, sep="\t", index_col=0)
    elif args.input_file.endswith(".csv"):
        auc_mtx = pd.read_csv(args.input_file, index_col=0)
    else:
        raise ValueError("Input file must be a TSV or CSV file")

    auc_mtx.apply(pd.to_numeric)
    # Parse threshold overrides
    threshold_overrides = None
    if args.threshold_overrides:
        import json

        threshold_overrides = json.loads(args.threshold_overrides)

    # Call binarize function
    binarized_mtx, thresholds = binarize(
        auc_mtx, threshold_overrides, args.seed, args.num_workers
    )

    # set column name for thresholds
    thresholds.rename("threshold", inplace=True)

    # Save output files
    binarized_mtx.to_csv(f"{args.output_prefix}/binarized_mtx.tsv", sep="\t")
    thresholds.to_csv(f"{args.output_prefix}/thresholds.tsv", sep="\t")
