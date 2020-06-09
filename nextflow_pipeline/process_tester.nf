#!/usr/bin/env nextflow

params.fastq_file_list = '/scratch/mblab/chasem/rnaseq_pipeline/job_scripts/nextflow_fastqfile_list_20200608_164325.csv' // this will need to be inputted by user
params.scratch_sequence = '/scratch/mblab/mblab.shared/scratch_sequence'

// split columns/rows of fastq_file_list for processing
Channel
    .fromPath(params.fastq_file_list)
    .splitCsv(header:true)
    .map{ row-> tuple(row.runDirectory, row.fastqFileName, row.organism, row.strandedness) }
    .set { samples_channel }

scratch_sequence = file(params.scratch_sequence)

// this works
process make_scratch_directory {
    cache = 'false'
    executor = 'local'

    input:
        set run_directory, fastq_filepath, organism, strandedness from samples_channel
    output:
        stdout into scratch_run_directory_ch

    script:
        """
        mkdir -p ${scratch_sequence}/${run_directory}
        echo ${realpath ${scratch_sequence}/${run_directory}}
        """
}

process rsync_from_lts {

    input:
        set scratch_run from scratch_run_directory_ch

    script:
        """
            echo ${scratch_run} >> /scratch/mblab/chasem/nextflow_output_tester.txt
        """

}

