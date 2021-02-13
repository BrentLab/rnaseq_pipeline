#!/usr/bin/env python3
"""
   For a given alignment/count job, create a nextflow config file. See templates for example.
   Currently, the expectation is that you run this on the cluster.
   create_nextflow_config.py -qs /path/to/database/sheet/describing/runs/to/process
"""
import sys
import argparse
import os
import numpy as np
import pandas as pd
from rnaseq_tools.DatabaseObject import DatabaseObject
from rnaseq_tools.OrganismDataObject import OrganismData
from rnaseq_tools import utils


def main(argv):
    """ main method
    :param argv: cmd line arguments
    """
    # parse cmd line arguments
    args = parseArgs(argv)
    query_sheet_path = args.query_sheet
    try:
        if not os.path.isfile(query_sheet_path):
            raise FileNotFoundError
    except FileNotFoundError:
        print('Query sheet path not valid. Check and try again.')
    try:
        interactive_flag = args.interactive
    except AttributeError:
        interactive_flag = False

    # instantiate DatabaseObject --> mostly this will be for access to StandardData paths
    db = DatabaseObject(query_sheet_path=query_sheet_path, config_file=args.config_file, interactive=interactive_flag)
    # read in dataframe
    db.query_df = utils.readInDataframe(db.query_sheet_path)
    # add column organism which identifies either KN99 or S288C_R64 depending on whether genotype1 starts with CNAG
    # TODO: this is point of weakness -- need to keep an eye here
    db.query_df['organism'] = np.where(db.query_df['genotype1'].str.startswith('CNAG'), 'KN99', 'S288C_R64')
    # cast libraryDate to datetime format
    db.query_df['libraryDate'] = pd.to_datetime(db.query_df['libraryDate'])
    # create strandedness column based on libraryDate. May change to prep protocol at some point, but for now this is best
    db.query_df['strandedness'] = np.where(db.query_df['libraryDate'] > '2015-10-25', 'reverse', 'no')
    # add leading zero to runNumber, if necessary -- take care of in loop
    db.query_df['runNumber'] = db.query_df['runNumber'].astype(str)
    db.query_df['fastqFileNumber'] = db.query_df['fastqFileNumber'].astype(str)
    # new dictionary to store run_directory in dataframe
    run_directory_list = []
    for index, row in db.query_df.iterrows():
        # some early runs have run numbers that start with zero in /lts. 0s are dropped in df b/c they are read in as ints
        # this step adds the zero and casts the row to str
        run_num_tmp = int(float(row['runNumber']))# TODO: super ugly, needs to be fixed. Not sure why this is now getting read in as 4422.0, eg as of 20200923
        if run_num_tmp in db._run_numbers_with_zeros: # TODO: Probably the best way to is to always read runnumbers as strings -- requires changing _run_num_with_zeros keys to strings, and checking the rest of the codebase that uses this
            run_number = str(db._run_numbers_with_zeros[run_num_tmp])
        else:
            run_number = run_num_tmp
        # create run directory name, eg run_1234_samples
        run_directory = 'run_' + str(run_number) + '_samples'  # SEE TODO above
        # add to list
        run_directory_list.append(run_directory)
        # create fastqfilename path
        try:
            fastq_filename = os.path.basename(row['fastqFileName']).rstrip()
        except TypeError:
            sys.exit("%s <-- not a fastqfilename?" %row['fastqFileName'])
        fastq_scratch_path = os.path.join(db.scratch_sequence, run_directory, fastq_filename)
        # move fastq file to scratch if it is not already tehre
        if not os.path.exists(fastq_scratch_path):
            fastq_lts_path = os.path.join(db.lts_sequence, run_directory, fastq_filename)
            scratch_run_directory_path = os.path.join(db.scratch_sequence, run_directory)
            utils.mkdirp(scratch_run_directory_path)
            print('...moving %s to %s' %(fastq_lts_path, scratch_run_directory_path))
            rsync_cmd = 'rsync -aHv %s %s' %(fastq_lts_path, scratch_run_directory_path)
            utils.executeSubProcess(rsync_cmd)
        # update fastqFileName in query_df
        db.query_df.loc[index, 'fastqFileName'] = fastq_scratch_path
    # add column runDirectory from run_directory_list
    db.query_df['runDirectory'] = run_directory_list

    # use OrganismDataObject to get paths to novoalign_index and annotation files
    kn99_organism_data = OrganismData(organism='KN99')
    kn99_novoalign_index = kn99_organism_data.novoalign_index
    # this is annotations + nc, t, r RNA with nc,t,r RNA annotations overlapping with protein coding ON SAME STRAND removed. rRNA retained
    kn99_annotation_file = kn99_organism_data.annotation_file
    # this is annotations + nc, t, r RNA with nc,t,r RNA annotations overlapping protein coding removed regardless of strand. rRNA retained
    kn99_annotation_file_no_strand = kn99_organism_data.annotation_file_no_strand
    kn99_genome = kn99_organism_data.genome
    s288c_r64_organism_data = OrganismData(organism='S288C_R64')
    s288c_r64_novoalign_index = s288c_r64_organism_data.novoalign_index
    s288c_r64_annotation_file = s288c_r64_organism_data.annotation_file
    s288c_r64_genome = s288c_r64_organism_data.genome

    # filter
    nextflow_fastqfile_df = db.query_df[['runDirectory', 'fastqFileName', 'organism', 'strandedness', 'fastqFileNumber']]
    for index, row in nextflow_fastqfile_df.iterrows():
        try:
            if not os.path.isfile(row['fastqFileName']):
                raise FileNotFoundError('fastqFileNotFoundInScratch')
        except FileNotFoundError:
            print('file %s was not successfully moved from lts to scratch' %row['fastqFileName'])
    print('\nnextflow fastq file .csv head:\n')
    print(nextflow_fastqfile_df.head())
    print('\n')
    # write out
    fastq_file_list_output_path = os.path.join(db.job_scripts,
                                               'nextflow_fastqfile_list' + '_' + args.name + '.csv')
    print('...writing out to %s' % fastq_file_list_output_path)
    nextflow_fastqfile_df.to_csv(fastq_file_list_output_path, index=False)

    # config_header goes at the top of the config -- includes date created and StandardObject instructions
    config_header = "/*\n" \
                    "* -------------------------------------------------\n" \
                    "*  Brentlab nextflow rnaseq_pipeline configuration\n" \
                    "* -------------------------------------------------\n" \
                    "* created with create_nextflow_config.py on %s\n" \
                    "* note: this is for a specific job for a specific user\n" \
                    "* and not intended as a general config file. To re-create\n" \
                    "* this job, you will need to run create_nextflow_config.py\n" \
                    "* with the same query_sheet input\n" \
                    "*/\n\n" % db.year_month_day

    # params section has all relevant path parameters to run the pipeline
    params_section = "// params necessary for the pipeline\n" \
                     "params {\n" \
                     "\tfastq_file_list = \"%s\"\n" \
                     "\tlts_sequence = \"%s\"\n" \
                     "\tscratch_sequence = \"%s\"\n" \
                     "\tlts_align_expr = \"%s\"\n" \
                     "\talign_count_results = \"%s\"\n" \
                     "\tlog_dir = \"%s\"\n" \
                     "\tKN99_novoalign_index = \"%s\"\n" \
                     "\tKN99_annotation_file = \"%s\"\n" \
                     "\tKN99_annotation_file_no_strand = \"%s\"\n" \
                     "\tKN99_genome = \"%s\"\n" \
                     "\tS288C_R64_novoalign_index = \"%s\"\n" \
                     "\tS288C_R64_annotation_file = \"%s\"\n" \
                     "\tS288C_R64_genome = \"%s\"\n" \
                     "}\n\n" % (fastq_file_list_output_path, db.lts_sequence, db.scratch_sequence,
                                db.lts_align_expr, db.align_count_results, db.log_dir, kn99_novoalign_index,
                                kn99_annotation_file, kn99_annotation_file_no_strand, kn99_genome, s288c_r64_novoalign_index, s288c_r64_annotation_file,
                                s288c_r64_genome)

    # write out and submit sbatch script with named/combined output/err

    nextflow_config_path = os.path.join(db.job_scripts, args.name + '_nextflow.config')
    print('...writing nextflow job config file to %s' % nextflow_config_path)
    with open(nextflow_config_path, 'w') as nextflow_config_file:
        nextflow_config_file.write(config_header)
        nextflow_config_file.write(params_section)

    sbatch_script_name = args.name + '_nextflow'
    nextflow_sbatch_path = os.path.join(db.job_scripts, sbatch_script_name + '.sbatch')
    # write sbatch script to submit nextflow job
    print('...writing sbatch script to %s' %nextflow_sbatch_path)
    with open(nextflow_sbatch_path, 'w') as nf_sbatch_file:
        nf_sbatch_file.write('#!/bin/bash\n'
                             '#SBATCH --mem=15G\n'
                             '#SBATCH -o %s/%s.out\n'
                             '#SBATCH -J %s\n\n'
                             'ml rnaseq_pipeline\n\n'
                             'nextflow -C %s run $CODEBASE/tools/align_count_pipeline.nf\n'
                             %(db.sbatch_log, sbatch_script_name, sbatch_script_name, nextflow_config_path))

    sbatch_cmd = 'sbatch %s' %nextflow_sbatch_path
    print('\nsubmitting sbatch script with cmd:\n\t%s' %sbatch_cmd)
    utils.executeSubProcess(sbatch_cmd)

    print('\nCheck progress by entering:\n\ttail %s/%s.out' %(db.sbatch_log, sbatch_script_name))
    print('\nTo run this in an interactive session, do the following:\n\t'
          'interactive\n\tnextflow -C %s run $CODEBASE/tools/align_count_pipeline.nf\n' % nextflow_config_path)
    print('If this job fails or is interrupted, you can resume it from where it failed by adding the flag -r to the nextflow command in the .sbatch file and resubmitting to sbatch')

def parseArgs(argv):
    parser = argparse.ArgumentParser(
        description="For a given alignment/count job, create a nextflow config file")
    parser.add_argument("-qs", "--query_sheet", required=True,
                        help='[REQUIRED] A .csv subset of metadata database describing the set of files you wish to process.\n'
                             'These will be grouped by run number')
    parser.add_argument("-n", "--name", required=True,
                        help="[REQUIRED] the name of the nextflow job -- this will be used to output the job script config file in\n"
                             "$USER/rnaseq_pipeline/jobs_scripts")
    parser.add_argument('--config_file', default='/see/standard/data/invalid/filepath/set/to/default',
                        help="[OPTIONAL] default is already configured to handle the invalid default path above in StandardDataObject.\n"
                             "Use this flag to replace that config file. Note: this is for StandardData, not nextflow")
    parser.add_argument('--interactive', action='store_true',
                        help="[OPTIONAL] set this flag (only --interactive, no input necessary) to tell StandardDataObject not\n"
                             "to attempt to look in /lts if on a compute node on the cluster")

    args = parser.parse_args(argv[1:])
    return args


if __name__ == "__main__":
    main(sys.argv)
