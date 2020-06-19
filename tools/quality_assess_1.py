#!/usr/bin/env python
import sys
import os
import argparse
from rnaseq_tools.QualityAssessmentObject import QualityAssessmentObject
from rnaseq_tools import utils

def main(argv):
    # parse cmd line arguments
    print('...parsing cmd line input')
    args = parseArgs(argv)
    try:
        if not os.path.isdir(args.reports_dir):
            raise NotADirectoryError('OutputDirDoesNotExist')
    except NotADirectoryError:
        print('%s does not lead to a valid directory. Check the path and resubmit with working -r' % args.reports)
    else:
        align_count_path = args.reports_dir
    try:
        if not os.path.isdir(args.output_dir):
            raise NotADirectoryError
    except NotADirectoryError:
        print('%s does not lead to a valid directory. check the path and resubmit with correct -o' % args.output_dir)
    else:
        output_directory = args.output_dir
    try:
        if not os.path.isfile(args.query_sheet_path):
            raise FileNotFoundError('QuerySheetDoesNotExist')
    except FileNotFoundError:
        print('%s does not lead to a valid file. Check and resubmit correct -qs' % args.query_sheet_path)
    except TypeError:
        pass

    # get run number if exists for output naming. if DNE, ask user to provide name to insert after run_<using run_num>_summary.csv
    try:
        run_number = utils.getRunNumber(args.reports_dir)
    except AttributeError:
        run_number = input(
            'No run number detected in input directory name. Enter something to insert in the output directory\n'
            'name: run_<what_you_input>_summary.csv: ')
    # store interactive flag
    try:
        interactive_flag = args.interactive
    except AttributeError:
        interactive_flag = False

    # TODO: ERROR CHECKING FOR GENOTYPE CHECK

    # create filename
    quality_assessment_filename = "run_{}_quality_summary.csv".format(run_number)
    output_path = os.path.join(output_directory, quality_assessment_filename)

    # for ordering columns below. genotype_1_coverage and genotype_2_coverage added if coverage_check is passed
    column_order = ['LIBRARY_SIZE', 'TOTAL_ALIGNMENT', 'UNIQUE_ALIGNMENT', 'MULTI_MAP', 'NO_MAP', 'HOMOPOLY_FILTER',
                    'READ_LENGTH_FILTER', 'NOT_ALIGNED_TOTAL', 'WITH_FEATURE', 'NO_FEATURE', 'FEATURE_ALIGN_NOT_UNIQUE',
                    'AMBIGUOUS_FEATURE', 'TOO_LOW_AQUAL']
    # if coverage_check is passed in cmd line, include query and coverage_check_flag in constructor (automatically sets some values #TODO make this a function with arugmnets to pass so as not to repeat entire constructor)
    if args.coverage_check:
        qa = QualityAssessmentObject(quality_assess_dir_path=align_count_path,
                                     log_suffix_list=["_novoalign.log", "_read_count.tsv"],
                                     run_number=run_number,
                                     output_dir=args.output_dir,
                                     quality_assessment_filename=quality_assessment_filename,
                                     coverage_check_flag=True,
                                     query_path = args.query_sheet_path,
                                     config_file=args.config_file,
                                     interactive=interactive_flag)
        column_order.extend(['GENOTYPE_1_COVERAGE', 'GENOTYPE_2_COVERAGE'])
    else:
        qa = QualityAssessmentObject(quality_assess_dir_path=align_count_path,
                                     log_suffix_list=["_novoalign.log", "_read_count.tsv"],
                                     run_number=run_number,
                                     output_dir=args.output_dir,
                                     quality_assessment_filename=quality_assessment_filename,
                                     config_file=args.config_file,
                                     interactive=interactive_flag)  # note: config_file is error checked in StandardDataObject. not necessary in this script

    print('...compiling alignment information')
    # create dataframes storing the relevant alignment and count metadata from the novoalign and htseq logs
    qual_assess_1_df = qa.compileData()

    # re_order columns
    qual_assess_1_df = qual_assess_1_df[column_order]

    print('writing output to %s' % output_path)
    qual_assess_1_df.to_csv(output_path, index_label="FASTQFILENAME")

def parseArgs(argv):
    parser = argparse.ArgumentParser(description="This script summarizes the output from pipeline wrapper.")
    parser.add_argument("-r", "--reports_dir", required=True,
                        help="[REQUIRED] Directory for novo and htseq log files. This currently only works for Novo and htseq output.")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="[REQUIRED] Strongly suggested Usage: in reports/run_####/{organism}_pipeline_info. Remember that runs with multiple organisms will have different pipeline_info dirs per organism."
                             " File path to the directory you wish to deposit the summary. Note: the summary will be called run_###_summary.csv")
    parser.add_argument("-cc", "--coverage_check", action='store_true',
                        help="[OPTIONAL] For Crypto experiments. Set this flag to add a column to the output dataframe with percent gene coverage")
    parser.add_argument("-qs", "--query_sheet_path",
                        help="[OPTIONAL] But required with -cc is set. Path to query sheet filtered for the files contained in the path passed to -r")
    parser.add_argument('--config_file', default='/see/standard/data/invalid/filepath/set/to/default',
                        help="[OPTIONAL] default is already configured to handle the invalid default path above in StandardDataObject.\n"
                             "Use this flag to replace that config file")
    parser.add_argument('--interactive', action='store_true',
                        help="[OPTIONAL] set this flag (only --interactive, no input necessary) to tell StandardDataObject not\n"
                             "to attempt to look in /lts if on a compute node on the cluster")
    args = parser.parse_args(argv[1:])
    return args


if __name__ == "__main__":
    main(sys.argv)
