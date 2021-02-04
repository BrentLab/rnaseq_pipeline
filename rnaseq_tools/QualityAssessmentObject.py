import os
import subprocess
import re
import configparser
import pandas as pd
import sys
from rnaseq_tools import utils
from rnaseq_tools.DatabaseObject import DatabaseObject
from rnaseq_tools.OrganismDataObject import OrganismData
import abc

# turn off SettingWithCopyWarning in pandas
pd.options.mode.chained_assignment = None


class QualityAssessmentObject(OrganismData):

    def __init__(self, expected_attributes=None, **kwargs):
        # add expected attributes to super._attributes
        self._add_expected_attributes = ['bam_file_list', 'count_file_list', 'novoalign_log_list',
                                         'coverage_check_flag', 'query_path', 'standardized_database_df',
                                         'qual_assess_dir_path']
        # This is a method of adding expected attributes to StandardData from StandardData children
        if isinstance(expected_attributes, list):
            self._add_expected_attributes.extend(expected_attributes)
        # initialize Standard data with the extended _attributes
        # recall that this will check for and/or create the directory structure found at
        super(QualityAssessmentObject, self).__init__(self._add_expected_attributes, **kwargs)
        # # set standardDirectory structure
        # self.standardDirectoryStructure() ## should already be done in StandardData constructor
        # overwrite super.self_type with object type of child (this object)
        self.self_type = 'QualityAssessmentObject'
        # create logger
        utils.createStandardObjectChildLogger(self, __name__)
        try:
            self.logger = None
            self.createOrganismDataLogger()
        except NotADirectoryError:
            print('Cannot create OrganismData logger. Check code and config file.')
            exit(1)
        else:
            if not os.path.isfile(self.log_file_path):
                raise FileNotFoundError('LoggerNotSuccessfullyCreated')
        # if query_path passed, read in as dataframe
        try:
            self.query_df = utils.readInDataframe(self.query_path)
            try:
                if not len(self.query_df[self.query_df.fastqFileName.isnull()]) == 0:
                    raise ValueError('EmptyFastqFileNamesIn %s' % self.query_path)
            except ValueError:
                self.logger.critical('Removing rows with null fastqFileNames from %s' % self.query_path)
                self.query_df = self.query_df[~self.query_df.fastqFileName.isnull()]
        except ValueError:
            self.logger.critical('query_path not valid')
        except FileNotFoundError:
            self.logger.critical('%s  --> query_path not valid' % self.query_path)
        except AttributeError:
            pass
        try:
            print('...extracting alignment information from novoalign logs')
            align_log_df = self.parseAlignmentLogs()
        except AttributeError:
            print("no novoalign files found")
        try:
            print('...extracting count information from htseq count files')
            if self.organism == 'KN99':
                count_summary_df = self.parseCountFiles(count_ambiguous_unique=True)
            else:
                count_summary_df = self.parseCountFiles()
        except AttributeError:
            print('no count files found')
        try:
            print('...compiling alignment and count information')
            self.qual_assess_df = self.compileAlignCountMetadata(align_log_df, count_summary_df)
            self.formatLibrarySizeColumns()
        except UnboundLocalError:
            print('no log summary or count summary')

    def formatLibrarySizeColumns(self):
        """
            format columns which are fractions of library size (this will be common to all organisms)
        """
        self.qual_assess_df['MULTI_MAP_PERCENT'] = self.qual_assess_df['MULTI_MAP'] / self.qual_assess_df[
            'LIBRARY_SIZE'].astype(
            'float')
        self.qual_assess_df['NO_MAP_PERCENT'] = self.qual_assess_df['NO_MAP'] / self.qual_assess_df[
            'LIBRARY_SIZE'].astype('float')
        self.qual_assess_df['HOMOPOLY_FILTER_PERCENT'] = self.qual_assess_df['HOMOPOLY_FILTER'] / self.qual_assess_df[
            'LIBRARY_SIZE'].astype('float')
        self.qual_assess_df['READ_LENGTH_FILTER_PERCENT'] = self.qual_assess_df['READ_LENGTH_FILTER'] / \
                                                            self.qual_assess_df[
                                                                'LIBRARY_SIZE'].astype('float')
        # htseq output not_aligned_total_percent is no_map + homopoly_filter + read_length filter. present as fraction of library_size
        self.qual_assess_df['NOT_ALIGNED_TOTAL_PERCENT'] = self.qual_assess_df['NOT_ALIGNED_TOTAL'] / \
                                                           self.qual_assess_df[
                                                               'LIBRARY_SIZE'].astype('float')

    def extractInfoFromQuerySheet(self, sample_name, extract_column):
        """ TODO: REMOVE THIS AND REPLACE WITH FUNCTION IN UTILS IN OTHER CODE
                  PUT INTO UTILS 20201220, BUT HAVE NOT REPLACED THE CODE EVERYWHERE
            //TODO: YES. THAT IS A GOOD IDEA. THIS IS HARD TO FIND WHEN COME ACROSS IN CHILDREN
            //TODO: PASS DF IN THE FUNCTION?
            extract information from query sheet given sample_name from qual_assess_df (which is the basename, no ext, of the fastq.gz)
            :param sample_name: name of sample -- basename of fastq.gz, no containing directory, no extension
            :param extract_column: column from which to extract a value from the query_df based on sample_name
            :returns: value extracted from query_df based on sample name and extract column
        """
        try:
            extract_value = \
            list(self.query_df[self.query_df['fastqFileName'].str.contains(sample_name + '.fastq.gz')][extract_column])[
                0]
        except AttributeError:
            self.logger.error("failure in extractInfoFromQuerySheet. sample_name: %s extract_column: %s" % (
            sample_name, extract_column))
            print('You must pass a query df')

        return str(extract_value)

    def parseAlignmentLogs(self):
        """

        """
        # instantiate dataframe
        align_df = pd.DataFrame()

        # extract metadata from novoalign log files
        for log_file in self.novoalign_log_list:
            # extract fastq filename
            fastq_basename = utils.pathBaseName(log_file).replace('_novoalign', '')
            # set sample name in library_metadata_dict
            library_metadata_dict = {"FASTQFILENAME": fastq_basename}
            print('...extracting information from novoalign log for %s' % fastq_basename)
            library_metadata_dict.update(self.parseAlignmentLog(log_file))
            align_df = align_df.append(pd.Series(library_metadata_dict), ignore_index=True)

        return align_df

    def parseCountFiles(self, count_ambiguous_unique=False):
        """
            get a list of the filenames in the run_#### file that correspond to a given type
            :param count_ambiguous_unique: boolean flag indicating whether to call uniqueAmbiguousProteinCodingCount()
            :returns: a dataframe containing the files according to their suffix
        """
        # instantiate dataframe
        htseq_count_df = pd.DataFrame()

        # extract metadata from count files
        for count_file in self.count_file_list:
            # extract fastq filename
            fastq_basename = utils.pathBaseName(count_file).replace('_read_count', '')
            # set sample name in library_metadata_dict
            library_metadata_dict = {"FASTQFILENAME": fastq_basename}
            print('...extracting count information from count file for %s' % fastq_basename)
            library_metadata_dict.update(self.parseGeneCount(count_file))
            if count_ambiguous_unique:
                library_metadata_dict['AMBIGUOUS_UNIQUE_PROTEIN_CODING_READS'] = self.uniqueAmbiguousProteinCodingCount(
                    fastq_basename)
            htseq_count_df = htseq_count_df.append(pd.Series(library_metadata_dict), ignore_index=True)

        return htseq_count_df

    def compileAlignCountMetadata(self, align_df, htseq_count_df):
        """
            Gather information from the novoalign logs and count files
        """
        # concat df_list dataframes together on the common column
        qual_assess_df = pd.merge(align_df, htseq_count_df, on='FASTQFILENAME')
        return qual_assess_df

    @abc.abstractmethod
    def auditQualAssessDataFrame(self, query_df_path, qual_assess_df, bam_file_list):
        """
            read in appropriate key/values from config file, add status and auto_audit columns to qual_asses_df
            :params qual_assess_df: the quality assessment dataframe with all appropriate columns for values in quality_assess thresholds in config file
            :returns: qual_assess_df with STATUS and AUTO_AUDIT columns added
        """
        raise NotImplementedError('AbstractMethodMustBeOverwrittenByOrganismSpecificQA')

    @abc.abstractmethod
    def quantifyNonCodingRna(self, qual_assess_df):
        """
            quantify amount of non coding RNA in library. Note: this requires annotation of non coding RNA in annotation file
            :param qual_assess_df: qual_assess_df with at least FASTQFILENAME column (or some column used to identify each file uniquely)
            :returns qual_assess_df: noncoding_rna_df with columns 'FASTQFILENAME','TOTAL_rRNA', 'UNIQUE_rRNA', 'UNIQUE_tRNA_ncRNA' added.
            NOTE: FASTQFILENAME may be some other way of consistently identifying a sample uniquely
        """
        raise NotImplementedError('AbstractMethodMustBeOverwrittenByOrganismSpecificQA')

    @abc.abstractmethod
    def formatQualAssessDataFrame(self, qual_assess_df):
        """
            A final step -- format columns in qual_assess_df for final writing
            :params qual_assess_df: a quality_assess_df with all columns
            :returns: qual_assess_df with appropriate formatting for writing out
        """
        raise NotImplementedError('AbstractMethodMustBeOverwrittenByOrganismSpecificQA')

    def parseAlignmentLog(self, alignment_log_file_path):
        """
            parse the information on the alignment out of a novoalign log
            :param alignment_log_file_path: the filepath to a novoalign alignment log
            :returns: a dictionary of the parsed data of the input file
        """
        library_metadata_dict = {}
        alignment_regex_dict = {'LIBRARY_SIZE': r"(?<=Read Sequences:\s)\s*\d*",
                                'UNIQUE_ALIGNMENT': r"(?<=Unique Alignment:\s)\s*\d*",
                                'MULTI_MAP': r"(?<=Multi Mapped:\s)\s*\d*",
                                'NO_MAP': r"(?<=No Mapping Found:\s)\s*\d*",
                                'HOMOPOLY_FILTER': r"(?<=Homopolymer Filter:\s)\s*\d*",
                                'READ_LENGTH_FILTER': r"(?<=Read Length:\s)\s*\d*"}

        # open the log path
        alignment_file = open(alignment_log_file_path, 'r')
        alignment_file_text = alignment_file.read()
        # loop over alignment_regex dict and enter values extracted from alignment_file into alignment_metadata_dict
        for alignment_category, regex_pattern in alignment_regex_dict.items():
            # extract the value corresponding to the alignment_category regex (see alignment_regex_dict)
            try:
                extracted_value = int(re.findall(regex_pattern, alignment_file_text)[0])
            except ValueError:
                msg = 'problem with file %s' % alignment_log_file_path
                self.logger.info(msg)
                print(msg)
            except IndexError:
                print('No %s in %s. Value set to 0' % (alignment_category, alignment_log_file_path))
                extracted_value = 0
            # check that the value is both an int and not 0
            if isinstance(extracted_value, int):
                library_metadata_dict.setdefault(alignment_category, extracted_value)
            else:
                print('cannot find %s in %s' % (alignment_category, alignment_log_file_path))

        # close the alignment_file and return
        alignment_file.close()

        return library_metadata_dict

    @abc.abstractmethod
    def parseGeneCount(self, htseq_counts_path):
        """
            count the gene counts that mapped either to genes (see COUNT_VARS at top of script for other features)
            :param htseq_counts_path: a path to a  _read_count.tsv file (htseq-counts output)
            :returns: a dictionary with the keys FEATURE_ALIGN_NOT_UNIQUE, TOO_LOW_AQUAL, AMBIGUOUS_FEATURE, NO_FEATURE, NOT_ALIGNED_TOTAL (at least, maybe more)
        """
        raise NotImplementedError('AbstractMethodMustBeOverwrittenByOrganismSpecificQA')

    @abc.abstractmethod
    def uniqueAmbiguousProteinCodingCount(self, fastq_simplename):
        """
            intersect bed file with protein coding coords with an alignment file with the htseq annotations added (in this case, grep out __ambiguous)
            :params fastq_simplename: fastq filename minus any path and extention eg /path/to/my_reads_R1.fastq.gz would be my_reads_R1
            :returns: the number of reads (lines) aligning to protein coding coordinates
        """
        raise NotImplementedError('AbstractMethodMustBeOverwrittenByOrganismSpecificQA')

    def totalrRNA(self, bam_path, rRNA_region, strandedness):
        """
            extract number of crypto aligning to rRNA
            :params bam_path: path to a bam file. requires that an index (with samtools index, novoindex, etc) be in same directory (will have .bai extension)
            :params rRNA_region: region of the rRNA in samtools view format -- chromosome:start_coord-stop_coord. eg CP022021.1:204523-235534
            :params strandedness: no or reverse, according to the library prep. This determines whether all
            :returns: total_rRNA is unique_rRNA + num_primary_alignment_rRNA (this is the number of primarily alignments from the multimap), AND unique_rRNA
        """
        # extract number of primary alignments of multimapped reads to rRNA. If reverse, only accept reads mapping in the positive direction
        # NOTE: THIS IS ONLY CORRECT IF THE rRNA IS ON THE FORWARD STRAND
        if strandedness == 'reverse':
            # -F 16 excludes reverse strand reads
            cmd_primary_multi_alignment_rRNA = 'samtools view -F 16 %s %s | grep ZS:Z:R | grep HI:i:1 | grep -v \"HI:i:1[[:digit:]]\" | wc -l' % (
                bam_path, rRNA_region)
        else:
            # grep -v excludes reads with a digit after the 1 (there is a prettier way to do this, im sure)
            cmd_primary_multi_alignment_rRNA = 'samtools view %s %s | grep ZS:Z:R | grep HI:i:1 | grep -v \"HI:i:1[[:digit:]]\" | wc -l' % (
                bam_path, rRNA_region)
        # as long as this is the first function called that needs an index, this will error check that samtools index has been run
        try:
            self.logger.debug(
                'samtools cmd to extract primary multi alignment reads to rRNA: %s' % cmd_primary_multi_alignment_rRNA)
            num_primary_alignment_rRNA = int(subprocess.getoutput(cmd_primary_multi_alignment_rRNA))
        except ValueError:
            sys.exit('You must first index the alignment files with samtools index')

        # extract number of unique alignments to rRNA
        # NOTE: THIS IS ONLY CORRECT IF THE rRNA IS ON THE FORWARD STRAND
        if strandedness == 'reverse':
            cmd_unique_rRNA = 'samtools view -F 16 %s %s | grep -v ZS:Z:R | wc -l' % (bam_path, rRNA_region)
        else:
            cmd_unique_rRNA = 'samtools view %s %s | grep -v ZS:Z:R | wc -l' % (bam_path, rRNA_region)
        self.logger.debug('samtools cmd to extract unique alignment reads to rRNA: %s' % cmd_unique_rRNA)
        unique_rRNA = int(subprocess.getoutput(cmd_unique_rRNA))

        # add for total rRNA
        total_rRNA = num_primary_alignment_rRNA + unique_rRNA

        return total_rRNA, unique_rRNA

    # TODO: CURRENTLY DOES NOT INCLUDE PRIMARY ALIGNMENT OF MULTI MAPS -- SHOULD IT?
    def totaltRNAncRNA(self, bam_path, trna_ncrna_annotation_gff, strandedness):
        """
            count reads that align to provided gff if the annotation overlaps at least 90% of the read uniquely
            :param bam_path: path to bam annotation file
            :param trna_ncrna_annotation_gff: a gff3 format with same chromosome/location notation as genome
            :param strandedness: strandedness of the library ('no' or 'reverse')
            :returns: number of reads fulfilling overlap criteria
        """
        # construct command based on strandedness of library
        if strandedness == 'reverse':
            # note: look up bedtools intersect --help for flags. grep -v returns all lines except those matching the pattern, in this case signifying multimaps
            bedtools_cmd = 'bedtools intersect -s -f .90 -a %s -b %s | samtools view | grep -v ZS:Z:R | wc -l' % (
                bam_path, trna_ncrna_annotation_gff)
        else:
            # no -s means this will count intersects regardless of strand
            bedtools_cmd = 'bedtools intersect -f .90 -a %s -b %s | samtools view | grep -v ZS:Z:R | wc -l' % (
                bam_path, trna_ncrna_annotation_gff)
        self.logger.info('bedtools cmd: %s' % bedtools_cmd)

        # extract unique_alignments to nc and t RNA
        unique_align_tRNA_ncRNA = int(subprocess.getoutput(bedtools_cmd))

        return unique_align_tRNA_ncRNA

    @abc.abstractmethod
    def calculateIntergenicCoverage(self, qual_assess_df):
        """
            calculate coverage of regions in genome between features
            Requires that the number of bases in the intergenic regions of the genome is present in OrganismData_config.ini
            in genome_files/<organism>
            TODO: BETTER NOTES HERE -- SEE CRYPTOQUALITYASSESSMENTOBJECT
        """
        raise NotImplementedError('AbstractMethodMustBeOverwrittenByOrganismSpecificQA')

    @abc.abstractmethod
    def calculateExonicCoverage(self, bam_file, output_directory):
        """
            calculate coverage of exon regions. deposits file in output directory named utils.pathBaseName(bam_file)+'_exonic_coverage.csv'
            :param bam_file: path to bam file
            :param output_directory: path to output directory
        """
        raise NotImplementedError('AbstractMethodMustBeOverwrittenByOrganismSpecificQA')

    def createCoverageBedFile(self, bam_file):
        """
             TODO: coverage is currently calculated separately for every function that needs it. Instead, calculate
             coverage over entire genome, store as bed in rnaseq_pipeline/tmp, and use it as input to any function requiring
             coverage calc
        """
        raise NotImplementedError

    def calculatePercentFeatureCoverage(self, feature, genotype, annotation_path, bam_file, read_strand=None,
                                        num_bases_in_region=None):
        """
            Calculate percent of given feature (regions summed, so all exon in gene, eg) of a gene (exon, CDS, etc) covered by 1 or more reads
            :param feature: annotation feature over which to take percentage, eg all exons in gene, or all CDS
            :param genotype: gene in annotation file
            :param annotation_path: path to annotation file
            :param bam_file: a sorted, indexed alignment file (.bam)
            :param num_bases_in_region: pass number of bases in the region directly, this will skip the step of calculating this number from the annotation file
            :param read_strand: default None, which will calculate coverage of reads on both strands over locus. options are Forward or Reverse. TODO: IMPLEMENT IN FUNCTION
            :returns: the fraction of bases in the (summed over the number of features in the gene) feature region covered by at least one read
        """
        if read_strand not in [None, "forward", "reverse"]:
            raise ValueError('UnrecognizedReadStrand')
        if not num_bases_in_region:
            # extract number of bases in CDS of given gene. Credit: https://www.biostars.org/p/68283/#390427
            num_bases_in_region_cmd = "grep %s %s | grep %s | bedtools merge | awk -F\'\t\' \'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}\'" % (
                genotype, annotation_path, feature)
            self.logger.info(' num bases in region cmd: %s' % num_bases_in_region_cmd)
            num_bases_in_region = int(subprocess.getoutput(num_bases_in_region_cmd))
        # extract number of bases with depth != 0
        num_bases_depth_not_zero_cmd = "grep %s %s | grep %s | gff2bed | samtools depth -aa -Q 10 -b - %s | cut -f3 | grep -v 0 | wc -l" % (
            genotype, annotation_path, feature, bam_file)
        self.logger.info(' num bases depth not zero over region cmd: %s' % num_bases_depth_not_zero_cmd)
        num_bases_in_cds_with_one_or_more_read = int(subprocess.getoutput(num_bases_depth_not_zero_cmd))

        return num_bases_in_cds_with_one_or_more_read / float(num_bases_in_region)

    def extractLog2cpm(self, gene, fastq_simple_name, log2cpm_csv_path):
        """
           extract log2cpm
           :param gene: log2cpm of gene you wish to extract
           :param fastq_simple_name: name of the fastq_file, no ext, no path
           :param log2cpm_csv_path: created by raw_count.py --> log2_cpm.R. index in gene_name, columns are fastq_simple_name_read_count.tsv
           :returns: log2cpm of a gene in a given library
        """
        # read in log2cpm dataframe
        try:
            log2cpm_df = utils.readInDataframe(log2cpm_csv_path)
            log2cpm_df = log2cpm_df.set_index('gene_id')
        except (FileNotFoundError, TypeError):
            error_msg = 'ERROR: log2cpm file not valid in some way check that it both exists and is in right format (eg, that the gene_id column is right): %s' % log2cpm_csv_path
            self.logger.critical(error_msg)
            sys.exit(error_msg)

        # create column name corresponding to log2cpm_df from fastq_simple_name
        column_name = fastq_simple_name + '_read_count.tsv'
        # check that it is actually in the log2cpm_df columns
        try:
            if column_name not in log2cpm_df.columns:
                raise AttributeError('ColumnNameNotInLog2CpmSheet')
        except AttributeError:
            error_msg = '%s not in log2cpm sheet %s' % (column_name, log2cpm_csv_path)
            self.logger.critical(error_msg)
            print(error_msg)
        else:
            # check that the gene is in the gene_id column
            try:
                if gene not in log2cpm_df.index:
                    gene = gene.replace('CNAG', 'CKF44')
                    if gene not in log2cpm_df.index:
                        raise AttributeError('GeneNotInLog2cpmSheet')
            except AttributeError:
                error_msg = '%s not in log2cpm sheet %s' % (gene, log2cpm_csv_path)
                self.logger.critical(error_msg)
                print(error_msg)

            else:
                # extract log2cpm and return
                return log2cpm_df.loc[gene, column_name]

    # TODO: DO NEW LOG2CPM REMAKE BY NEW TREATMENT COLUMNS. INCLUDE STRAINS FOR WILDTYPES
    def foldOverWildtype(self, perturbed_gene, sample_name, log2cpm_path, sample_treatment, sample_timepoint):
        """
            :param perturbed_gene:
            :param sample_name:
            :param log2cpm_path:
            :param sample_treatment: a list containing [medium, temperature, atmosphere] TODO: THIS NEEDS TO BE CHANGED FOR CLARITY SINCE THERE IS ALSO A TREATMENT COLUMN
            :param sample_timepoint:

        """
        # read in median_wt_expression_by_timepoint_treatment_df
        try:
            median_wt_expression_by_timepoint_treatment_df = utils.readInDataframe(
                self.median_wt_expression_by_timepoint_treatment)
            # SET INDEX ON (gene_id, TREATMENT, TIMEPOINT) note: timepoint is read in as an int
            median_wt_expression_by_timepoint_treatment_df = median_wt_expression_by_timepoint_treatment_df.set_index(
                ['gene_id', 'TREATMENT', 'TIMEPOINT'])
        except AttributeError:
            self.logger.critical('genome files config in constructor did not work')
            print('genome files config in constructor did not work')  # set this as attr in crypto organismData
        # extract overexpression log2_cpm # TODO: READ IN LOG2_CPM AS ATTR OF QUAL ASSESS
        overexpression_log2cpm = float(self.extractLog2cpm(perturbed_gene, sample_name, log2cpm_path))
        # get wildtype log2_cpm from median_wt-expression_by_timepoint_treatment
        try:
            wt_log2cpm = float(median_wt_expression_by_timepoint_treatment_df.loc[
                                   (perturbed_gene, sample_treatment, int(sample_timepoint)), 'MEDIAN_LOG2CPM'])
        except KeyError:
            perturbed_gene = perturbed_gene.replace('CNAG', 'CKF44')
            wt_log2cpm = float(median_wt_expression_by_timepoint_treatment_df.loc[
                                   (perturbed_gene, sample_treatment, int(sample_timepoint)), 'MEDIAN_LOG2CPM'])

        overexpression_fow = overexpression_log2cpm - wt_log2cpm

        return overexpression_fow

    def indexBamFileBathScript(self, bam_files_to_index):
        """

        """
        raise NotImplementedError('AbstractMethodMustBeOverwritten')

    def igvShot(self, metadata_df, igv_output_dir=None, manual_audit_flag=False):
        """
            first checks if alignment files are indexed; creates and submit sbatch script to index if not.
            creates a lookup_file.tsv of alignment_file(bam)_path \t bed_path \t igv_genome_path \t igv_output_dir
            sbatch file uses lookup file to input to this command:

            make_IGV_snapshots.py ${bam_file} -bin /opt/apps/igv/2.4.7/igv.jar -nf4 -r ${bed_file} -g ${igv_genome} -fig_format png -o ${igv_outout_dir}

            This is written to be included as part of the rnaseq qc pipeline. alignment files are expected to be in align_count_results/align
            # NOTE: THESE MUST BE INDEXED ALREADY

            :param igv_output_dir: default None, and the output will be deposited in reports/igv_<datetime>. enter this to redirect output elsewhere
            :param metadata_df: a metadata sheet -- note, note a qual assessment. this won't have auto_audit inforamtion
            :param manual_audit_flag: default false. If true, column manualAudit will be used to select samples for igv
        """
        # set output_dir to default if not passed by user
        if not igv_output_dir:
            igv_output_dir = os.path.join(self.reports, 'igv_%s_%s' % (self.year_month_day, utils.hourMinuteSecond()))
            utils.mkdirp(igv_output_dir)

        lookup_file_path = self.createIgvBatchScripts(metadata_df, igv_output_dir)

        # create sbatch script for igv shots
        self.createAndSubmitIgvSbatchScript(lookup_file_path, igv_output_dir)

    def createIgvBatchScripts(self, metadata_df, igv_output_dir, batchscript_lookup_full_path=None, organism="KN99",
                              gene_offset=500, kn99_marker_check=True):
        """
            See templates/igv_batchscript_example.txt for an example batchscript and some notes
            Previously, this repo was used:
                  https://github.com/stevekm/IGV-snapshot-automator
            For documentation on automatic control of the IGV, see:
            # see https://software.broadinstitute.org/software/igv/batch
                  http://software.broadinstitute.org/software/igv/automation
              for coloring the reads, see the issues discussion here: https://github.com/stevekm/IGV-snapshot-automator/issues/14
              which links to the IGV repo here:
                  https://github.com/igvteam/igv/blob/master/src/main/resources/org/broad/igv/prefs/preferences.tab

            The purpose of this script is to create a directory in job_scripts that holds a number of batchscripts (See the SEE note above)
            Each batchscript will tell IGV to take a snapshot of a pane with both the perturbed sample and a WT (perturbed on top, WT on bottom) and separate images of the markers

            These images will be output into a directory with the following structure (this the topmost directory is made in igvShot, the subdirectories for each sample are made here)
            igv_<year_month_day>_<hour_minute_second>
                fastq_simplename1
                    CNAG_12345.batch
                    CNAG_12345.png
                    NAT.png
                    G418.png
                fastq_simplename2
                    ...

            the topmost directory name can be set by the user in igvShots. fastq_simplename1, etc will be taken from the metadata for the perturbed samples.
            Each subdirectory will have the batchscript, and the images.
            The batchscript is worth keeping b/c it can be used to load the locus on your local computer easily.
            It can also obviously be used to create the browser shots again -- all that would need to be updated is the path to the igv_genome and the bam (make sure all the bams are indexed if you do this)
            It would be a good idea to remove the stuff other than the load statement when doing this.

            :params metadata_df: metadata_df for which igv snapshots will be made for perturbed samples
            :params igv_output_dir: where to place the output
            :params batchscript_lookup_full_path: default is None. If None, the filename will be self.job_scripts, igv_batchscript_lookup_<year_month_day>_<hour_minute_second>.txt
            :params organism: right now, only KN99. This is in preparation of generalizing
            :params gene_offset: the amount left/right of the gene to include in the browser shot window
            :params kn99_marker_check: whether or not to take browser shots of the marker, also
            :creates: subdirectories in igv_output_dir, each with a batchscript for that particular perturbed sample (names of subdirectories will be fastq file without .fastq.gz)
            :returns: a lookup to all of the batchscripts in the subdirectories above. This will be in job_scripts
        """
        if batchscript_lookup_full_path is None:
            batchscript_lookup_full_path = os.path.join(self.job_scripts, 'igv_%s_%s.txt' % (
            self.year_month_day, utils.hourMinuteSecond()))
        # get all paths from StandardData, etc
        organism_genome_files = os.path.join(self.genome_files, organism)
        if not os.path.isdir(organism_genome_files):
            self.logger.critical('%s not a in genome_files. Make sure genome_files exists in rnaseq_pipeline, '
                                 'and that the organism is passed as a correctly formatted subdirectory of genome_files' % organism_genome_files)
            raise NotADirectoryError('OrganismDirectoryNotFoundInGenomeFiles')
        organism_config_file = os.path.join(organism_genome_files, 'OrganismData_config.ini')
        if not os.path.isfile(organism_config_file):
            self.logger.critical(
                '%s not found -- check genome_files/organism subdir. Possibly delete genome_files and let the script re-download (make sure it is accessible)' % organism_config_file)
            raise FileNotFoundError('ConfigFileNotFound')
        # read in organism config file
        organism_config = configparser.ConfigParser()
        organism_config.read(organism_config_file)
        organism_config_dict = organism_config['OrganismData']
        # set annotation file
        annotation_file = os.path.join(self.genome_files, organism, organism_config_dict['annotation_file'])
        if not os.path.isfile(annotation_file):
            self.logger.critical('the annotation file %s not valid filepath' % annotation_file)
            raise FileNotFoundError('AnnotationFileNotFound')
        # set igv genome
        try:
            igv_genome = os.path.join(self.genome_files, organism, organism_config_dict['igv_genome'])
        except KeyError:  # this is for KN99
            igv_genome = os.path.join(self.genome_files, organism, organism_config_dict['igv_stranded_genome'])
        if not os.path.isfile(igv_genome):
            self.logger.critical('igv_genome path %s not valid' % igv_genome)
            raise FileNotFoundError('IgvGenomeNotFound')

        marker_list = ["CNAG_NAT", "CNAG_G418"]
        kn99_marker_dict = dict(zip(marker_list, self.createIgvBedLine(marker_list, annotation_file, gene_offset)))

        # the name is a remnant of the usage of the repo specific in the docstring above. To turn this dict into an actual bed entry, just join the list with \t
        # Here, this is used to create the igvBatchScript rather than going through the extra step of writing a bed.
        bed_entry_dict = {}
        for index, row in metadata_df.iterrows():
            genotype_list = utils.extractGenotypeList(row,
                                                      convert_CNAG_to_CKF44=True)  # last argument to convert CNAG to CKF44
            if genotype_list[0] != 'CKF44_00000' and genotype_list[0] is not None:
                run_num = str(self.extractRunNumber(int(float(row[
                                                                  "runNumber"]))))  # TODO: this will hopefully be improved when database is integrated w/sql. no idea why runNumber column now comes out with inconsistent type. should just be string
                self.logger.debug("runnumber extracted by igv func: %s" % run_num)
                fastq_simple_name = utils.pathBaseName(row["fastqFileName"])
                alignment_run_directory_path = os.path.join(self.align_count_results, "run_%s_samples/align" % run_num)
                try:
                    print("...looking for wildtype reference for %s" % fastq_simple_name)
                    wt_reference_list = self.getWildtypeReference(row)
                    wt_reference_bam_path = wt_reference_list[
                        0]  # item 1 is a list of the conditions -- use this to provide information in the subdir?
                except FileNotFoundError:
                    self.logger.critical("No wildtype reference found for %s" % row)
                try:
                    if not os.path.isfile(wt_reference_bam_path):
                        raise FileNotFoundError
                except FileNotFoundError:
                    self.logger.critical("path to wildtype reference DNE: %s" % wt_reference_bam_path)
                try:
                    if not os.path.isdir(alignment_run_directory_path):
                        raise NotADirectoryError("alignment_dir_DNE")
                except NotADirectoryError:
                    self.logger.critical("%s DNE" % alignment_run_directory_path)
                bam_file_path = os.path.join(alignment_run_directory_path,
                                             fastq_simple_name + "_sorted_aligned_reads_with_annote.bam")
                try:
                    if not os.path.isfile(bam_file_path):
                        raise FileNotFoundError("bamDNE")
                except FileNotFoundError:
                    self.logger.critical("%s DNE" % bam_file_path)
                try:
                    if not (os.path.isfile(os.path.join(alignment_run_directory_path,
                                                        fastq_simple_name + "_sorted_aligned_reads_with_annote.bam.bai"))):
                        raise FileNotFoundError("baiDNE")
                except FileNotFoundError:
                    self.logger.critical("%s index DNE" % bam_file_path)

                # create this dictionary to error check and exit if there is a problem. Do not go onto creating batchfiles until this step has passed.
                bed_line_list = self.createIgvBedLine(genotype_list, annotation_file, gene_offset)
                # fill bed_entry_dict
                # TODO: THE FASTQFILES NEED TO BE CHECKED WELL BEFORE THIS -- DO SO AND GET RID OF THE REDUNDANT LOOP HERE
                try:
                    if fastq_simple_name in bed_entry_dict.keys():
                        raise KeyError
                    bed_entry_dict[fastq_simple_name] = {"perturbed_genotype": genotype_list, "igv_genome": igv_genome,
                                                         "perturbed_bam": bam_file_path,
                                                         "wt_bam": wt_reference_bam_path,
                                                         "perturbed_bed_list": bed_line_list}
                except KeyError:
                    sys.exit(
                        "There is a repeated fastqFileName in the metadata. This is a problem -- fix it. Either deduplicate the metadata, or figure out why a fastq is identified in two different records")

        batchfile_list = []
        for sample_name, batch_file_dict in bed_entry_dict.items():
            # make subdir of igv_output_dir
            igv_output_subdir_path = os.path.join(igv_output_dir, sample_name)
            utils.mkdirp(igv_output_subdir_path)
            # make batchfile path
            batchfile_list.extend(
                self.writeIgvBatchfile(igv_output_subdir_path, batch_file_dict, marker_dict=kn99_marker_dict))

        with open(batchscript_lookup_full_path, "w") as batchscript_lookup_file:
            batchscript_lookup_file.write("\n".join(batchfile_list))
            batchscript_lookup_file.write("\n")

        return batchscript_lookup_full_path

    def writeIgvBatchfile(self, output_dir, batch_file_dict, control_image=True, marker_dict=None):
        """
            yet another site to see for port commands for igv: https://software.broadinstitute.org/software/igv/PortCommands
            write a batch file -- see templates/igv_batchfile_example.txt
            :params output_dir: directory in which to write the file
            :params batch_file_dict: {"perturbed_genotype": [genotype1, genotype2,...], "igv_genome": igv_genome, "perturbed_bam": bam_file_path, "wt_bam": wt_reference_bam_path, "perturbed_bed_list": [bed_line_info_genotype1, bed_line_info_genotype2,...]}
                                    note: wt only required if wt_image = True, which is True by default
            :params control_image: True by default. adds a control_image batchscript to take an image of the control at same locus as perturbed
            :params marker_dict: default is None. For KN99, eg, pass {NAT: [bed_line_list], G418: [bed_line_list]} (See createIgvBedLine())
            :returns: path to the file this function writes
        """
        perturbed_filename = ''
        batchscript_cmd_list = ["new", "snapshotDirectory", "load", "goto", "sort", "collapse", "snapshot"]

        # this will be the first line of the batchscript
        # batchfile_text = ["snapshotDirectory %s\n" %output_dir,
        #                   "setSleepInterval 1000\n\n"]
        batchfile_text_dict = {"perturbed": [], "control": []}

        # note -- the value associated with "perturbed_genotype" needs to be passed as a list, even if only one item
        # TODO: ERROR CHECK THAT THERE IS A BED LINE FOR EACH GENOTYPE
        # create perturbed sample batchscript paragraph -- this will have two alignment files loaded. The top is the perturbed
        for genotype in batch_file_dict["perturbed_genotype"]:
            if genotype is not None:
                # get index to extract the correct bed_file_list (this is in the event of multiple genotypes)
                genotype_index = batch_file_dict["perturbed_genotype"].index(genotype)
                bed_line = batch_file_dict["perturbed_bed_list"][genotype_index]
                # note: same format as input to IGV viewer: chromosome:start-stop
                perturbed_locus_bed_line = bed_line[0] + ":" + str(bed_line[1]) + "-" + str(bed_line[2])
                batchfile_text_dict["perturbed"].extend(["new\n",
                                                         "snapshotDirectory %s\n" % output_dir,
                                                         "genome %s\n" % batch_file_dict["igv_genome"],
                                                         "maxPanelHeight 500\n",
                                                         "preference SAM.COLOR_BY READ_STRAND\n",
                                                         "load %s\n" % batch_file_dict["perturbed_bam"],
                                                         "goto %s\n" % perturbed_locus_bed_line,
                                                         "sort position\n",
                                                         "collapse\n",
                                                         "snapshot %s\n" % (genotype + ".png")
                                                         # TODO: ADD CONDITIONS OF SAMPLE
                                                         ])
                # write marker paragraphs for the batch scripts. note: unlike the perturbed sample, only the perturbed alignment file is used
                if marker_dict:
                    for marker, marker_bed_line in marker_dict.items():
                        marker = re.sub("^CNAG_", "",
                                        marker)  # TODO: this is purpose built for KN99 right now -- figure out how to generalize
                        # note: same format as input to IGV viewer: chromosome:start-stop
                        marker_bed_line = marker_bed_line[0] + ":" + str(marker_bed_line[1]) + "-" + str(
                            marker_bed_line[2])
                        batchfile_text_dict["perturbed"].extend(["goto %s\n" % marker_bed_line,
                                                                 "sort position\n",
                                                                 "collapse\n",
                                                                 "snapshot %s\n" % (marker + ".png")
                                                                 ])
                batchfile_text_dict["perturbed"].extend(["exit\n\n"])
                if control_image:
                    batchfile_text_dict["control"].extend(["new\n",
                                                           "snapshotDirectory %s\n" % output_dir,
                                                           "genome %s\n" % batch_file_dict["igv_genome"],
                                                           "maxPanelHeight 500\n",
                                                           "preference SAM.COLOR_BY READ_STRAND\n",
                                                           "load %s\n" % batch_file_dict["wt_bam"],
                                                           "goto %s\n" % perturbed_locus_bed_line,
                                                           "sort position\n",
                                                           "collapse\n",
                                                           "snapshot %s\n" % ("%s_control.png" % genotype)
                                                           # TODO: RENAME WITH CONDITIONS OF WT
                                                           ])
                    batchfile_text_dict["control"].extend(["exit\n\n"])
        # write out
        # TODO: WRITE ANOTHER BATCHFILE TO THE DIRECTORY THAT OMITS THE SNAPSHOT STUFF SO THAT A USER CAN LOAD THE LOCUS DIRECTLY ON THEIR LOCAL WITH THAT NEW BATCHSCRIPT
        try:
            # try to drop any None from list
            sample_genotype = batch_file_dict["perturbed_genotype"].remove(None)
        except ValueError:
            pass
        else:
            sample_genotype = batch_file_dict["perturbed_genotype"]
        # create list to store list of filenames to batchfiles
        igv_batchfilename_list = []
        for key, batchfile_text in batchfile_text_dict.items():  # poorly named -- batchfile_text_dict too close to batch_file_dict
            # as long as there is something to write
            if len(batchfile_text) > 0:
                # create a filename
                igv_batchfilename = os.path.join(output_dir, "".join(sample_genotype) + '_' + key + ".txt")
                # append the filename
                igv_batchfilename_list.append(igv_batchfilename)
                print("...writing batchfile to %s" % igv_batchfilename)
                with open(igv_batchfilename, "w") as batchfile:
                    batchfile.write("".join(batchfile_text))

        return igv_batchfilename_list

    def createIgvBedLine(self, genotype_list, annotation_file, gene_offset):
        """
            TODO: this does just create a bed line in the form of a list -- generalize and put into Utils
                  Rather than making both genotypes, make one genotype at a time. Arguments to "make bed line"
                  should be "locus, annotation_file, gene_offset". To create the output that this currently creates,
                  one would submit it for each gentoype
            The name of this is a remnant of the old pipeline usage of an external repo for the IGV snapshot automation.
            This function returns a list in the following format:
               [chromosome_identifier, start_coord_with_offest, stop_coord_with_offset] eg: [chr1, 1, 1000]
            To turn this into an actual bed line, simply join with tabs, eg "\t".join(bed_line)

            :params genotype_list: genotype list, [genotype1, genotype2]. if single ko, [genotype1, None]
            :returns: a list of lists in the format specified above. If single KO, [line1, None]
        """
        bed_line_list = []
        for gene in genotype_list:
            if gene is not None and gene != "CNAG_00000":
                # get first column corresponding to given gene, take uniq value as chromosome
                extract_chr_cmd = 'grep %s %s | cut -f1 | uniq' % (gene, annotation_file)
                chromosome_identifier = subprocess.getoutput(extract_chr_cmd)

                # get the list of all start coordinates associated with a feature, sort, and take smallest
                extract_start_coord_cmd = 'grep %s %s | cut -f4 | sort -n | head -1' % (gene, annotation_file)
                start_coord = int(subprocess.getoutput(extract_start_coord_cmd))
                start_coord_with_offset = start_coord - gene_offset  # add offset to widen window on igv_viewer
                if start_coord_with_offset < 1:
                    start_coord_with_offset = 1
                # sort stop coordinates of all features of given gene, take largest
                extract_stop_coord_cmd = 'grep %s %s | cut -f5 | sort -n | tail -1' % (gene, annotation_file)
                stop_coord = int(subprocess.getoutput(extract_stop_coord_cmd))
                stop_coord_with_offset = stop_coord + gene_offset  # add offset to widen window on igv_viewer

                # enter {gene: bed_line} to dict -- not the final \t allows to add the bam_file_simplename later
                bed_line = [chromosome_identifier, start_coord_with_offset, stop_coord_with_offset]
                bed_line_list.append(bed_line)

        return bed_line_list

    def createAndSubmitIgvSbatchScript(self, lookup_file_path, output_dir, igv_mem="10000",
                                       igv_path="/opt/apps/igv/2.4.7/igv.jar"):  # TODO: USE THE IGV IN THE LAB /OPT /opt/apps/labs/mblab/software/IGV_Linux_2.8.13/lib/igv.jar # working igv: /opt/apps/igv/2.4.7/igv.jar
        """
            create sbatch script -- there is an of batchscripts in templates (though there may not be one of this exactly if i forget to put it there. feel free to write an issues report)
            :params lookup_file_path: a list of batchscript files
            :params output_dir: the directory into which the individual batchscripts will write (what actually matters is what is in the batchscript. this only checks to make sure it exists)
            :params igv_mem: the amount of memory to give the igv program
            :params igv_path: path to the igv.jar
            :creates: a sbatch file for this particular list of batchfiles
            :submits: to the scheduler on the htcf cluster
        """
        # error check input
        if not os.path.isfile(lookup_file_path):
            self.logger.critical('lookup file %s does not exist' % lookup_file_path)
            raise FileNotFoundError('LookupFileDoesNotExist')
        if not os.path.isdir(output_dir):
            self.logger.critical('output dir %s does not exist' % output_dir)
            raise NotADirectoryError('OutputDirectoryDoesNotExist')

        # write sbatch job. see https://htcfdocs.readthedocs.io/en/latest/runningjobs/
        line_count_cmd = 'cat %s | wc -l' % lookup_file_path
        line_count = int(subprocess.getoutput(line_count_cmd))
        sbatch_array_line = "--array=1-{}%{}".format(line_count, min(line_count,
                                                                     20))  # this has to be 1 at a time since the WT may be the same for multiple samples
        job = '#!/bin/bash\n\n' \
              '#SBATCH -N 1\n' \
              '#SBATCH --mem=20G\n' \
              '#SBATCH %s\n' \
              '#SBATCH -o %s/igv_snapshot_%s_%s.out\n' \
              '#SBATCH -J igv_snapshot\n\n' \
              'module load java\n\n' \
              'read igv_batchfile < <(sed -n ${SLURM_ARRAY_TASK_ID}p %s )\n\n' \
              'xvfb-run --auto-servernum --server-num=1 java -Xmx%sm -jar %s -b ${igv_batchfile}\n' \
              % (sbatch_array_line, self.sbatch_log, self.year_month_day, utils.hourMinuteSecond(), lookup_file_path,
                 igv_mem, igv_path)

        # write to file
        # check that a job_scripts directory for igv for today has been created, make one if no
        job_script_dir_path = os.path.join(self.job_scripts, 'igv_%s' % self.year_month_day)
        if not os.path.isdir(os.path.join(job_script_dir_path)):
            utils.mkdirp(job_script_dir_path)

        igv_job_script_path = os.path.join(job_script_dir_path,
                                           'igv_%s_%s.sbatch' % (self.year_month_day, utils.hourMinuteSecond()))
        with open(igv_job_script_path, 'w') as igv_job_script_file:
            igv_job_script_file.write(job)

        # error check igv_job_script_path
        if not os.path.isfile(igv_job_script_path):
            self.logger.critical('Failed to write job %s to file %s' % (job, igv_job_script_path))
            raise FileNotFoundError('IgvJobScriptNotSuccessfullyWritten')
        else:
            print('Submitting igv batchscript %s' % igv_job_script_path)
            batchscript_submit_cmd = 'sbatch %s' % igv_job_script_path
            utils.executeSubProcess(batchscript_submit_cmd)

    def qortsPlots(self):
        """
            create the qorts multiplots
        """
        try:  # TODO: clean this up
            if not os.path.isfile(self.query_path):
                raise FileExistsError('QueryPathNotValid')
            else:
                if not hasattr(self, 'standardized_database_df'):
                    query_df = utils.readInDataframe(self.query_path)
                    self.standardized_database_df = DatabaseObject.standardizeDatabaseDataframe(query_df)
        except FileExistsError or AttributeError:
            print('Query Path not valid')

        try:
            if not os.path.isdir(self.qorts_output):
                raise NotADirectoryError('QortsOutputDoesNotExist')
        except NotADirectoryError or AttributeError:
            print('qorts_output is either not passed as an argument, or does not exist.')

        try:
            if not hasattr(self, 'organism'):
                raise AttributeError('NoOrganismPassed')
        except AttributeError:
            print('The QualityAssessmentObject must have an attribute organism')

        try:
            if not os.path.isfile(self.annotation_file):
                raise FileExistsError('AnnotationFilePathNotValid')
        except (FileExistsError, AttributeError):
            organism_object = OrganismData(organism=self.organism)

        try:
            if not hasattr(self, 'bam_file_list'):
                raise AttributeError('NoBamFileList')
        except AttributeError:
            print('You must supply a list of bamfiles')

        # sort list of countfilenames pre and post 2015
        self.standardized_database_df['LIBRARYDATE'] = pd.to_datetime(self.standardized_database_df['LIBRARYDATE'])
        pre_2015_df = self.standardized_database_df[self.standardized_database_df['LIBRARYDATE'] <= '2015-01-01']
        pre_2015_sample_list = list(pre_2015_df['FASTQFILENAME'])
        post_2015_df = self.standardized_database_df[self.standardized_database_df['LIBRARYDATE'] > '2015-01-01']
        post_2015_sample_list = list(post_2015_df['FASTQFILENAME'])

        cmd_list = []

        for fastqfile_simplename in post_2015_sample_list:
            bamfilename_path = [bam_file for bam_file in self.bam_file_list if fastqfile_simplename in bam_file][0]
            try:
                if bamfilename_path not in self.bam_file_list:
                    raise FileNotFoundError('QuerySampleNotInAlignCountDirectory')
            except FileNotFoundError:
                return ('A sample in the query could not be located in the directory with the alignment files.\n'
                        'Make sure the query provided corresponds to the align_count_path directory')

            output_subdir = os.path.join(self.qorts_output, fastqfile_simplename + '_qorts')
            utils.mkdirp(output_subdir)

            try:
                if not os.path.exists(bamfilename_path):
                    raise FileExistsError('BamfilePathNotValid')
            except FileExistsError:
                print('path from align_count_path to bamfile does not exist')
            if not hasattr(self, 'annotation_file'):
                annotation_file = organism_object.annotation_file
            else:
                annotation_file = self.annotation_file

            try:
                if not os.path.exists(bamfilename_path):
                    raise FileExistsError('BamfilePathNotValid')
            except FileExistsError:
                print('path from align_count_path to bamfile does not exist')
            qorts_cmd = 'java -Xmx1G -jar /opt/apps/labs/mblab/software/QoRTs/hartleys-QoRTs-099881f/scripts/QoRTs.jar QC --singleEnded --stranded --keepMultiMapped --generatePlots %s %s %s\n' % (
                bamfilename_path, annotation_file, output_subdir)
            cmd_list.append(qorts_cmd)

        for fastqfile_simplename in pre_2015_sample_list:
            bamfilename_path = [bam_file for bam_file in self.bam_file_list if fastqfile_simplename in bam_file][0]
            try:
                if bamfilename_path not in self.bam_file_list:
                    raise FileNotFoundError('QuerySampleNotInAlignCountDirectory')
            except FileNotFoundError:
                print('A sample in the query could not be located in the directory with the alignment files.\n'
                      'Make sure the query provided corresponds to the align_count_path directory')

            output_subdir = os.path.join(self.qorts_output, fastqfile_simplename + '_qorts')
            utils.mkdirp(output_subdir)

            try:
                if not os.path.exists(bamfilename_path):
                    raise FileExistsError('BamfilePathNotValid')
            except FileExistsError:
                print('path from align_count_path to bamfile does not exist')
            if not hasattr(self, 'annotation_file'):
                if hasattr(organism_object, 'annotation_file_no_strand'):
                    annotation_file = organism_object.annotation_file_no_strand
                else:
                    annotation_file = organism_object.annotation_file
            else:
                annotation_file = self.annotation_file

            qorts_cmd = 'java -Xmx16G -jar /opt/apps/labs/mblab/software/QoRTs/hartleys-QoRTs-099881f/scripts/QoRTs.jar QC --singleEnded --keepMultiMapped --generatePlots %s %s %s\n' % (
                bamfilename_path, annotation_file, output_subdir)
            cmd_list.append(qorts_cmd)

        sbatch_script = os.path.join(self.job_scripts,
                                     'qorts_job_script_' + self.year_month_day + '_' + utils.hourMinuteSecond() + '.sbatch')
        print('writing job script to %s' % sbatch_script)
        with open(sbatch_script, 'w') as sbatch_file:
            sbatch_file.write('\n')
            sbatch_file.write('\n')

            for cmd in cmd_list:
                sbatch_file.write(cmd)
                sbatch_file.write('\n')
