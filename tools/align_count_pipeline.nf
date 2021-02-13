#!/usr/bin/env nextflow

// split columns/rows of fastq_file_list for processing
Channel
    .fromPath(params.fastq_file_list)
    .splitCsv(header:true)
    .map{row-> tuple(row.runDirectory, file(row.fastqFileName), row.organism, row.strandedness, row.fastqFileNumber) }
    .set { fastq_filelist }


scratch_sequence = file(params.scratch_sequence)

process novoalign {

    scratch true
    executor "slurm"
    cpus 8
    memory "20G"
    beforeScript "ml novoalign samtools"
    stageInMode "copy"
    stageOutMode "move"
    publishDir "$params.align_count_results/$run_directory/logs", mode:"copy", overwite: true, pattern: "*.log"


    input:
        tuple val(run_directory), file(fastq_file), val(organism), val(strandedness), val(fastq_file_number) from fastq_filelist
    output:
        tuple val(fastq_file_number), val(run_directory), val(fastq_simple_name), val(organism), val(strandedness), file("${fastq_simple_name}_sorted_aligned_reads.bam") into bam_align_ch
        tuple val(run_directory), file("${fastq_simple_name}_novoalign.log"), file("${fastq_simple_name}_novosort.log") into novoalign_log_ch

    script:
        fastq_simple_name = fastq_file.getSimpleName()
        if (organism == 'S288C_R64')
            """
            novoalign -r All \\
                      -c 8 \\
                      -o SAM \\
                      -d ${params.S288C_R64_novoalign_index} \\
                      -f ${fastq_file} 2> ${fastq_simple_name}_novoalign.log | \\
            samtools view -bS | \\
            novosort - \\
                     --threads 8 \\
                     --markDuplicates \\
                     -o ${fastq_simple_name}_sorted_aligned_reads.bam 2> ${fastq_simple_name}_novosort.log

            """
        else if (organism == 'KN99')
            """
            novoalign -r All \\
                      -c 8 \\
                      -o SAM \\
                      -d ${params.KN99_novoalign_index} \\
                      -f ${fastq_file} 2> ${fastq_simple_name}_novoalign.log | \\
            samtools view -bS | \\
            novosort - \\
                     --threads 8 \\
                     --markDuplicates \\
                     --index \\
                     -o ${fastq_simple_name}_sorted_aligned_reads.bam 2> ${fastq_simple_name}_novosort.log
            """
        else if (organism == 'H99')
            """
            novoalign -r All \\
                      -c 8 \\
                      -o SAM \\
                      -d ${params.H99_novoalign_index} \\
                      -f ${fastq_file} \\
                      2> ${fastq_simple_name}_novoalign.log | \\
            samtools view -bS | \\
            novosort - \\
                     --threads 8 \\
                     --markDuplicates \\
                     -o ${fastq_simple_name}_sorted_aligned_reads.bam 2> ${fastq_simple_name}_novosort.log

            """
}

process htseq_count {

    scratch true
    executor "slurm"
    cpus 8
    memory "20G"
    beforeScript "ml samtools htseq"
    stageInMode "copy"
    stageOutMode "move"
    publishDir "$params.align_count_results/$run_directory/logs", mode:"copy", overwite: true, pattern: "*.log"
    publishDir "$params.align_count_results/$run_directory/count", mode:"copy", overwite: true, pattern: "*_read_count.tsv"
    publishDir "$params.align_count_results/$run_directory/align", mode:"copy", overwite: true, pattern: "*_sorted_aligned_reads_with_annote.bam"


    input:
        tuple val(fastq_file_number), val(run_directory), val(fastq_simple_name), val(organism), val(strandedness), file(sorted_bam) from bam_align_ch
    output:
        tuple val(run_directory), val(fastq_simple_name), val(organism), val(strandedness), file("${fastq_simple_name}_sorted_aligned_reads_with_annote.bam") into bam_align_with_htseq_annote_ch
        tuple val(fastq_file_number), val(run_directory), val(fastq_simple_name), file("${fastq_simple_name}_read_count.tsv") into htseq_count_ch
        tuple val(run_directory), val(fastq_simple_name), file("${fastq_simple_name}_htseq.log") into htseq_log_ch
        tuple val(run_directory), val(organism), val(strandedness) into pipeline_info_ch

    script:
        if (organism == 'S288C_R64')
            """
            htseq-count -f bam \\
                        -o ${fastq_simple_name}_htseq_annote.sam \\
                        -s ${strandedness} \\
                        -t gene \\
                        -i ID \\
                        ${sorted_bam} \\
                        ${params.S288C_R64_annotation_file} \\
                        1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log

            sed "s/\t//" ${fastq_simple_name}_htseq_annote.sam > ${fastq_simple_name}_no_tab_sam.sam

            samtools view --threads 8 ${sorted_bam} | \\
            paste - ${fastq_simple_name}_no_tab_sam.sam | \\
            samtools view --threads 8 -bS -T ${params.S288C_R64_genome} > ${fastq_simple_name}_sorted_aligned_reads_with_annote.bam

            """
        else if (organism == 'KN99' && strandedness == 'reverse')
            """
            htseq-count -f bam \\
                        -o ${fastq_simple_name}_htseq_annote.sam \\
                        -s ${strandedness} \\
                        -t exon \\
                        -i gene \\
                        ${sorted_bam} \\
                        ${params.KN99_annotation_file} \\
                        1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log

            sed "s/\t//" ${fastq_simple_name}_htseq_annote.sam > ${fastq_simple_name}_no_tab_sam.sam

            samtools view --threads 8 ${sorted_bam} | \\
            paste - ${fastq_simple_name}_no_tab_sam.sam | \\
            samtools view --threads 8 -bS -T ${params.KN99_genome} > ${fastq_simple_name}_sorted_aligned_reads_with_annote.bam

            """
        else if (organism == 'KN99' && strandedness == 'no')
            """
            htseq-count -f bam \\
                        -o ${fastq_simple_name}_htseq_annote.sam \\
                        -s ${strandedness} \\
                        -t exon \\
                        -i gene \\
                        ${sorted_bam} \\
                        ${params.KN99_annotation_file_no_strand} \\
                        1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log

            sed "s/\t//" ${fastq_simple_name}_htseq_annote.sam > ${fastq_simple_name}_no_tab_sam.sam

            samtools view --threads 8 ${sorted_bam} | \\
            paste - ${fastq_simple_name}_no_tab_sam.sam | \\
            samtools view --threads 8 -bS -T ${params.KN99_genome} > ${fastq_simple_name}_sorted_aligned_reads_with_annote.bam

            """
        else if (organism == 'H99')
            """
            htseq-count -f bam \\
                        -s ${strandedness} \\
                        -t exon \\
                        ${sorted_bam} \\
                        ${params.H99_annotation_file} \\
                        1> ${fastq_simple_name}_read_count.tsv 2> ${fastq_simple_name}_htseq.log

            sed "s/\t//" ${fastq_simple_name}_htseq_annote.sam > ${fastq_simple_name}_no_tab_sam.sam

            samtools view --threads 8 ${sorted_bam} | \\
            paste - ${fastq_simple_name}_no_tab_sam.sam | \\
            samtools view --threads 8 -bS -T ${params.H99_genome} > ${fastq_simple_name}_sorted_aligned_reads_with_annote.bam

            """
}

process postHtseqCountsToDatabase {

    executor "local"
    beforeScript "ml rnaseq_pipeline"

    input:
        tuple val(fastq_file_number), val(run_directory), val(fastq_simple_name), file(read_count_tsv) from htseq_count_ch

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    import requests
    import json
    from urllib.request import HTTPError
    from rnaseq_tools import utils
    from rnaseq_tools.StandardDataObject import StandardData

    #from urllib.request import urlopen, HTTPError
    #try:
    #html = urlopen('http://www.test.com')
    #except HTTPError as e:
    #print(e)
    #else:
    #print('its ok')

    sd = StandardData(interactive=True)
    # TODO: SET PARAMS FOR OVERWRITE = T/F AND ALTERNATE BTWN PUT/POST DEPENDING ON ERROR RESPONSE

    # read count_file in as pandas df
    count_file_df = pd.read_csv("${read_count_tsv}", sep='\t', names=['htseq_col', "${fastq_simple_name}"])

    # remove the htseq stuff at the bottom of the file (it starts with __, eg __ambiguous)
    gene_counts = count_file_df[~count_file_df.htseq_col.str.startswith("__")]

    # maybe push this into its own channel?
    htseq_log_data = count_file_df[count_file_df.htseq_col.str.startswith("__")] # TODO: SEND THIS INTO HTSEQ QC CHANNEL? Handle here?

    # get the count dict in structure {fastqFileName: [counts]}
    count_dict = count_file_df.drop(['htseq_col'], axis=1).to_dict(orient="list")

    # this is the body of the request. fastqFileNumber is the foreign key of Counts
    data = {'fastqFileNumber': "${fastq_file_number}", 'rawCounts': json.dumps(count_dict)}

    # TODO: make this a parameter
    url = 'http://13.59.167.2/api/Counts'

    try:
        utils.postData(url, data, sd.logger)
    except(HTTPError):
        sd.logger.warning('COUNTS http post failed on fastq: %s, fastqfilenumber: %s' %("${fastq_simple_name}", "${fastq_file_number}"))
    """



}

process writePipelineInfo {

    executor "local"
    beforeScript "ml rnaseq_pipeline"

    input:
        tuple val(run_directory), val(organism), val(strandedness) from pipeline_info_ch

    script:
"""
#!/usr/bin/env python

from rnaseq_tools.OrganismDataObject import OrganismData
from rnaseq_tools import utils
import os

# instantiate OrganismDataObject (see brentlab rnaseq_pipeline)
od = OrganismData(organism = "${organism}", interactive=True)

# create pipeline_info subdir of in rnaseq_pipeline/align_count_results/${organism}_pipeline_info
pipeline_info_subdir_path = os.path.join(od.align_count_results, "${run_directory}", "${organism}_pipeline_info")
utils.mkdirp(pipeline_info_subdir_path)

# write version info from the module .lua file (see the .lua whatis statements)
pipeline_info_txt_file_path = os.path.join(pipeline_info_subdir_path, 'pipeline_info.txt')
cmd_pipeline_info = 'module whatis rnaseq_pipeline 2> %s' %pipeline_info_txt_file_path
utils.executeSubProcess(cmd_pipeline_info)

# include the date processed in pipeline_info_subdir_path/pipeline_into.txt
with open(pipeline_info_txt_file_path, "a+") as file:
    file.write('')
    current_datetime = od.year_month_day + '_' + utils.hourMinuteSecond()
    file.write('Date processed: %s' %current_datetime)
    file.write('')

# set annotation_file
if "${organism}" == 'KN99' and "${strandedness}" == 'no':
      annotation_file = od.annotation_file_no_strand
else:
    annotation_file = od.annotation_file
# include the head of the gff/gtf in pipeline_info
cmd_annotation_info = 'head %s >> %s' %(annotation_file, pipeline_info_txt_file_path)
utils.executeSubProcess(cmd_annotation_info)

# TODO: try copying nextflow jobscript to pipeline_info

"""
}
