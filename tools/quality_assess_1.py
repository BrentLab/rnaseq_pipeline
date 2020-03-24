#!/usr/bin/env python
import sys
import os
import re
import argparse
from glob import glob
import pandas as pd
from rnaseq_tools.StandardData import StandardData
from rnaseq_tools import utils

# current order is required -- order is hard coded into parseAlignment and parseGeneCount
ALIGN_VARS = ["Read Sequences", "Unique Alignment", "Multi Mapped", "No Mapping Found"]
COUNT_VARS = ["total_mapped_reads", "with_feature", "no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique"]

def main(argv):
    """
    qual_assess_1 main method
    :param argv: cmd line input -ac path to align_counts; -rn run number; -o output directory
    """
    # parse cmd line arguments
    args = parseArgs(argv)
    # create StandardDataFormat object
    sdf = StandardData(align_count_path=args.reports, run_number=utils.getRunNumber(args.reports), output_dir=args.output)

    # create path to new quality_assessment_sheet
    quality_assessment_filename = "run_{}_quality_summary.csv".format(sdf.run_number)
    output_csv = os.path.join(sdf.output_dir, quality_assessment_filename)

    # create dataframes storing the relevant alignment and count metadata from the novoalign and htseq logs
    alignment_files_df = compileData(sdf.align_count_path, "_novoalign.log")
    count_files_df = compileData(sdf.align_count_path, "_read_count.tsv")

    # concat the alignment and count data
    combined_df = pd.concat([alignment_files_df, count_files_df], axis=1, sort=True, join="inner")
    # write to csv
    combined_df.to_csv(output_csv, columns=ALIGN_VARS + COUNT_VARS[1:], index_label="Sample")

    # TODO: make genotype check automatic
    # prompt user to enter gene_list if genotype check fails
    # set current to 0 for ko, <50 percentile for _over
    # make labels both genotype and fastqfilename in R script
    # take browser shot -- build this as class (certain attributes necessary to taking browser shot)

    # genotype check
    if args.genotype_check:
        sdf.query_sheet_path = args.query_sheet_path
        sdf.log2_cpm_path = args.log2cpm
        sdf.experiment_columns = args.exp_columns
        utils.genotypeCheck(sdf)


def parseArgs(argv):
    parser = argparse.ArgumentParser(description="This script summarizes the output from pipeline wrapper.")
    parser.add_argument("-r", "--reports", required=True,
                        help="[REQUIRED] Directory for alignment log files. This currently only works for Novoalign output.")
    parser.add_argument("-o", "--output", required=True,
                        help="[REQUIRED] Suggested Usage: in reports/run_####/{organism}_pipeline_info. Remember that runs with multiple organisms will have different pipeline_info dirs per organism."
                             " File path to the directory you wish to deposit the summary. Note: the summary will be called run_###_summary.csv")
    parser.add_argument("-gc", "--genotype_check", action='store_true',
                        help="path to a text file with a list (see templates/genotype_check.txt for an example) of fastqFileNames to perform a genotype check on")
    parser.add_argument("-qs", "--query_sheet",
                        help="path to query sheet with (only) the samples you wish to genotype check")
    parser.add_argument("-log2cpm", "--log2_cpm",
                        help="log2 cpm of the raw counts (run raw_counts.py and then log2_cpm.R to produce these")
    parser.add_argument("-exp_cols", "--experiment_columns", nargs='+',
                        help="columns to use to create replicate groups. Eg for crypto -exp_cols genotype treatment timepoint"
                             " note: no quotes or commas")
    args = parser.parse_args(argv[1:])
    return args


def compileData(dir_path, suffix):
    df = pd.DataFrame()
    file_paths = glob("{}/*{}".format(dir_path, suffix))
    for file_path in file_paths:
        fastq_filename = re.findall(r'(.+?){0}'.format(suffix), os.path.basename(file_path))[0]
        if "novoalign" in suffix:
            data = parseAlignmentLog(file_path)
        elif "read_count" in suffix:
            data = parseGeneCount(file_path)
        data.update({"Sample": fastq_filename})
        df = df.append(pd.Series(data), ignore_index=True)
    return df.set_index("Sample")


def parseAlignmentLog(alignment_log_file_path):
    with open(alignment_log_file_path, "r") as file:
        alignment_metadata_dict = {}
        for line in file:
            line = line.strip("#").strip()
            for k in ALIGN_VARS:
                if line.startswith(k):
                    v = int(line.split(":")[-1].split("(")[0].strip())
                    alignment_metadata_dict[k] = v if k == ALIGN_VARS[0] else v / float(
                        alignment_metadata_dict[ALIGN_VARS[0]])
    return alignment_metadata_dict


def parseGeneCount(file_path):
    """
    count the gene counts that mapped either to genes or to other features (see COUNT_VARS at top of script for other features)
    :param file_path: a path toa  _read_count.tsv file (htseq output)
    :returns: a dictionary with the features in COUNT_VARS quantified for the inputted file
    """
    with open(file_path, "r") as file:
        # COUNT_VARS[1] is WITH_FEATURE
        out = {COUNT_VARS[1]: 0}
        # walk through the count file
        for line in file:
            # split the two columns and store values as index and value
            index, value = line.strip().split()
            # if the line starts with __, it stores information such as __no_feature
            if line.startswith("__"):
                out[index[2:]] = int(value)
            else:
                # otherwise the line will store individual gene counts. Add these counts to WITH_FEATURE
                out[COUNT_VARS[1]] += int(value)
        # COUNT_VARS[0] = TOTAL_MAPPED_READS and is the sum across the rows of the counts in the various categories
        out[COUNT_VARS[0]] = sum([out[key] for key in COUNT_VARS[1:]])
        # divide each category starting with WITH_FEATURE by the TOTAL_MAPPED_READS
        for key in COUNT_VARS[1:]:
            out[key] /= float(out[COUNT_VARS[0]])
    return out


if __name__ == "__main__":
    main(sys.argv)
