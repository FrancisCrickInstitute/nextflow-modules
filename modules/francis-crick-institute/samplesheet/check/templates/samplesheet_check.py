#!/usr/bin/env python3

import sys
import argparse
import shutil


def check_samplesheet(samplesheet_path, output_path):
    """
    Checks input samplesheet
    """

    with open(samplesheet_path, "r", encoding="UTF-8") as fin:
        # Read header
        header = [x.strip('"') for x in fin.readline().strip().split(",")]

        # Cycle through samplesheet line by line
        line_no = 2
        for line in fin:
            data_line = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check if its just a blank line so we dont error
            if line.strip() == "":
                continue

            # Check valid number of columns per row
            if len(data_line) != len(header):
                print(f"Invalid number of columns (found {len(data_line)} should be {len(header)})! - line no. {line_no}")
                sys.exit(1)

            line_no = line_no + 1

    # Copy to output path
    shutil.copy(samplesheet_path, output_path)

def main(args):
    """
    Main input function
    """

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--process_name", default="!{process_name}")
    parser.add_argument("--sample", default="!{sample}")
    parser.add_argument("--output", default="!{output}")
    args = parser.parse_args(args)

    # Run main function
    check_samplesheet(args.sample, args.output)

if __name__ == "__main__":
    main(sys.argv[1:])
