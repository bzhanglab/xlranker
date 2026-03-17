"""This script is for generating the PPI DB for generic tables where one column is Symbol and the other is Homolog."""

import argparse
import gzip
import os
import pickle

import polars as pl


def create_ppi_db(input_path, output_path):
    df = pl.read_csv(input_path, separator="\t")
    with gzip.open(output_path, "wb") as f:
        final = dict(
            zip(
                df["Symbol"],
                df["Homolog"],
            )
        )
        # print first 5 entries
        print(list(final.items())[:5])
        pickle.dump(
            final,
            f,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate PPI DB from a generic table with Symbol and Homolog columns."
    )
    parser.add_argument("input", help="Input TSV file path")
    parser.add_argument("output", help="Output PKL.GZ file path")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    create_ppi_db(args.input, args.output)
