import os
import sys
import subprocess
import pandas as pd
from rnaseq_tools import utils
from rnaseq_tools.QualityAssessmentObject import QualityAssessmentObject

# turn off SettingWithCopyWarning in pandas
pd.options.mode.chained_assignment = None

# TODO: CALCULATE COVERAGE ONCE, STORE AS BED FILE, USE BED RATHER THAN QUANTIFYING BAM EVERYTIME
# TODO: EXTRACT FILEPATHS FOR EG LOG2CPM MUCH MORE CLEARLY, ERROR CHECK (put this in QualityAssessObject, eg)
class CryptoQualityAssessmentObject(QualityAssessmentObject):

    def __init__(self, expected_attributes=None, **kwargs):
        # add expected attributes to super._attributes
        self._add_expected_attributes = []
        # This is a method of adding expected attributes to StandardData from StandardData children
        if isinstance(expected_attributes, list):
            self._add_expected_attributes.extend(expected_attributes)
        # initialize Standard data with the extended _attributes
        # recall that this will check for and/or create the directory structure found at
        super(CryptoQualityAssessmentObject, self).__init__(self._add_expected_attributes, **kwargs)
        # overwrite super.self_type with object type of child (this object)
        self.self_type = 'CryptoQualityAssessmentObject'
        # create logger
        self.logger = utils.createStandardObjectChildLogger(self, __name__)
        # for ordering columns below. genotype1_coverage and genotype2_coverage added if coverage_check is passed
        self.column_order = ['FASTQFILENAME', 'LIBRARY_SIZE', 'EFFECTIVE_LIBRARY_SIZE', 'EFFECTIVE_UNIQUE_ALIGNMENT',
                             'EFFECTIVE_UNIQUE_ALIGNMENT_PERCENT',
                             'MULTI_MAP_PERCENT', 'PROTEIN_CODING_TOTAL', 'PROTEIN_CODING_TOTAL_PERCENT',
                             'PROTEIN_CODING_COUNTED',
                             'PROTEIN_CODING_COUNTED_PERCENT', 'AMBIGUOUS_FEATURE_PERCENT', 'NO_FEATURE_PERCENT',
                             'INTERGENIC_COVERAGE', 'NOT_ALIGNED_TOTAL_PERCENT', 'GENOTYPE1_COVERAGE', 'GENOTYPE1_LOG2CPM',
                             'GENOTYPE2_COVERAGE', 'GENOTYPE2_LOG2CPM', 'OVEREXPRESSION_FOW', 'NAT_COVERAGE', 'NAT_LOG2CPM', 'G418_COVERAGE',
                             'G418_LOG2CPM', 'NO_MAP_PERCENT', 'HOMOPOLY_FILTER_PERCENT', 'READ_LENGTH_FILTER_PERCENT',
                             'TOO_LOW_AQUAL_PERCENT', 'rRNA_PERCENT', 'nctrRNA_PERCENT']

        print('Quantifying noncoding rRNA (rRNA, tRNA and ncRNA)')
        # extract rRNA, tRNA and ncRNA quantification for crypto from bam files -- this takes a long time
        ncRNA_df = self.quantifyNonCodingRna(self.qual_assess_df)
        # merge this into the self.qual_assess_df
        self.qual_assess_df = pd.merge(self.qual_assess_df, ncRNA_df, on='FASTQFILENAME')
        print('Quantifying intergenic coverage')
        self.qual_assess_df = self.calculateIntergenicCoverage(self.qual_assess_df)
        # if coverage_check_flag true, check coverage of perturbed genes
        try:
            if self.coverage_check_flag:
                coverage_df = self.perturbedCheck()
                self.qual_assess_df = pd.merge(self.qual_assess_df, coverage_df, how='left', on='FASTQFILENAME')
        except AttributeError:
            self.logger.info('query_df or coverage_check_flag not present -- no coverage check')
        # format the self.qual_assess_df dataframe
        self.qual_assess_df = self.formatQualAssessDataFrame(self.qual_assess_df)

    def formatQualAssessDataFrame(self, qual_assess_df):
        """
            Format/calculate column values for qual_assess_df
            :param qual_assess_df: a complete qual_assess_df (after all steps in cryptoQualityAssessmentObject have been run)
            :returns: a calculated/re-formatted qual_assess_df
        """
        # EFFECTIVE_LIBRARY_SIZE is LIBRARY_SIZE - (total_rRNA + unique_tRNA_ncRNA)
        qual_assess_df['EFFECTIVE_LIBRARY_SIZE'] = qual_assess_df['LIBRARY_SIZE'].astype('float') - (
                    qual_assess_df['TOTAL_rRNA'] + qual_assess_df['UNIQUE_tRNA_ncRNA'])

        # EFFECTIVE_UNIQUE_ALIGNMENT is number of unique reads minus those unique reads mapping to rRNA and nc + t RNA
        qual_assess_df['EFFECTIVE_UNIQUE_ALIGNMENT'] = qual_assess_df['UNIQUE_ALIGNMENT'] - (
                    qual_assess_df['UNIQUE_rRNA'] +
                    qual_assess_df['UNIQUE_tRNA_ncRNA'] +
                    qual_assess_df['TOO_LOW_AQUAL'])

        # present the following categories as fraction of library_size
        # percent of library made up of rRNA (recall total_rRNA is unique + primary alignment since reads multimap in two spots, both seemingly rRNA)
        qual_assess_df['rRNA_PERCENT'] = qual_assess_df['TOTAL_rRNA'] / qual_assess_df['LIBRARY_SIZE'].astype('float')
        qual_assess_df['nctRNA_PERCENT'] = qual_assess_df['UNIQUE_tRNA_ncRNA'] / qual_assess_df['LIBRARY_SIZE'].astype(
            float)
        qual_assess_df['nctrRNA_PERCENT'] = (qual_assess_df['TOTAL_rRNA'] + qual_assess_df['UNIQUE_tRNA_ncRNA']) / \
                                            qual_assess_df['LIBRARY_SIZE'].astype(float)

        # PROTEIN_CODING_TOTAL (formerly with_feature) is the number of reads mapping to a protein coding gene by htseq plus the unique ambiguous reads mapping to exon portiosn of overlapping protein coding reads
        qual_assess_df['PROTEIN_CODING_TOTAL'] = qual_assess_df['PROTEIN_CODING_COUNTED'] + qual_assess_df[
            'AMBIGUOUS_UNIQUE_PROTEIN_CODING_READS']

        # protein_coding_total as percent of effective unique alignment
        qual_assess_df['PROTEIN_CODING_TOTAL_PERCENT'] = qual_assess_df['PROTEIN_CODING_TOTAL'] / qual_assess_df[
            'EFFECTIVE_UNIQUE_ALIGNMENT'].astype('float')
        # protein_coding_counted as percent of unique alignment
        qual_assess_df['PROTEIN_CODING_COUNTED_PERCENT'] = qual_assess_df['PROTEIN_CODING_COUNTED'] / qual_assess_df[
            'EFFECTIVE_UNIQUE_ALIGNMENT'].astype('float')

        # present the following as fraction of (total) unique_alignment  ### TODO: SHOULD THIS BE OF EFFECTIVE_UNIQUE_ALIGNMENT? PROBABLY YES: NOTE THE GRAPHS IN CURRENT OLD CRYPTO SHEETS ARE OF UNIQUE_ALIGNMENT
        qual_assess_df['NO_FEATURE_PERCENT'] = qual_assess_df['NO_FEATURE'] / qual_assess_df[
            'EFFECTIVE_UNIQUE_ALIGNMENT'].astype(float)
        qual_assess_df['AMBIGUOUS_FEATURE_PERCENT'] = qual_assess_df['AMBIGUOUS_FEATURE'] / qual_assess_df[
            'EFFECTIVE_UNIQUE_ALIGNMENT'].astype(float)
        qual_assess_df['TOO_LOW_AQUAL_PERCENT'] = qual_assess_df['TOO_LOW_AQUAL'] / qual_assess_df[
            'EFFECTIVE_UNIQUE_ALIGNMENT'].astype(float)

        # present EFFECTIVE_UNIQUE_ALIGNMENT as percent of library size (make sure this is the last step
        qual_assess_df['EFFECTIVE_UNIQUE_ALIGNMENT_PERCENT'] = qual_assess_df['EFFECTIVE_UNIQUE_ALIGNMENT'] / \
                                                               qual_assess_df['LIBRARY_SIZE'].astype(float)
        # below is a messy way of ensuring that all expected columns are present, even if there is not a value (eg, if a sample has no over expression)
        for column in self.column_order:
            try:
                x = qual_assess_df[column]
            except KeyError:
                qual_assess_df[column] = None

        return qual_assess_df[self.column_order]

    def quantifyNonCodingRna(self, qual_assess_df):
        """

        """
        num_reads_to_ncRNA_dict = {}
        # set threshold to determine strandedness. note that there is an email from holly to yiming mentioning 10.25.2015. That is the best record we have, if htis message remains
        strandedness_date_threshold = pd.to_datetime('10.01.2015')
        kn99_tRNA_ncRNA_annotations = os.path.join(self.genome_files, 'KN99', 'ncRNA_tRNA_no_rRNA.gff')
        if hasattr(self, 'query_df'):
            for index, row in qual_assess_df.iterrows():
                try:
                    # extract genotype1
                    genotype = [self.extractInfoFromQuerySheet(row['FASTQFILENAME'], 'genotype1'), None]
                    #genotype = [list(self.query_df[self.query_df['fastqFileName'].str.contains(row['FASTQFILENAME'] + '.fastq.gz')]['genotype1'])[0]]
                except ValueError:
                    self.logger.info('genotype cannot be extracted with the fastq filename in this row. Note: if there are null entries in the column fastqFileNames, this is the cause. those need to be remedied or removed in order for this to work: %s' %row)
                try:
                    # extract genotype2 or set it to None
                    genotype[1] = self.extractInfoFromQuerySheet(row['FASTQFILENAME'], 'genotype2')
                except KeyError:
                    self.logger.debug("sample: %s does not have genotype2" %row['FASTQFILENAME'])
                # test if organism is KN99. Proceed if so
                if genotype[0].startswith('CNAG'):
                    # extract fastq_filename without any preceeding path or file extension
                    fastq_simple_name = utils.pathBaseName(row['FASTQFILENAME'])
                    print('...evaluating ncRNA in %s' % fastq_simple_name)
                    # use this to extract bam_path
                    try:
                        # bam_file_list is inherited
                        bam_path = [bam_file for bam_file in self.bam_file_list if fastq_simple_name in bam_file][0]
                    except IndexError:
                        self.logger.info('%s not in bam_file_list' % fastq_simple_name)
                        continue
                    try:
                        if not os.path.isfile(bam_path):
                            raise FileNotFoundError
                    except FileNotFoundError:
                        self.logger.error('bam file not found %s' % bam_path)
                        print('bam file not found: %s' % bam_path)
                    library_date = list(self.query_df[self.query_df['fastqFileName'].str.contains(row['FASTQFILENAME'] + '.fastq.gz')]['libraryDate'])[0]
                    row_date_time = pd.to_datetime(library_date)
                    strandedness = 'no' if row_date_time < strandedness_date_threshold else 'reverse'
                    total_rRNA, unique_rRNA = self.totalrRNA(bam_path, 'CP022322.1:272773-283180', strandedness)
                    unique_tRNA_ncRNA = self.totaltRNAncRNA(bam_path, kn99_tRNA_ncRNA_annotations, strandedness)
                    num_reads_to_ncRNA_dict.setdefault(fastq_simple_name,
                                                       {'total_rRNA': total_rRNA, 'unique_rRNA': unique_rRNA,
                                                        'total_tRNA_ncRNA': unique_tRNA_ncRNA})

        # create dataframe from num_reads_to_ncRNA_dict
        ncRNA_df = pd.DataFrame.from_dict(num_reads_to_ncRNA_dict, orient='index').reset_index()
        # format column headers
        ncRNA_df.columns = ['FASTQFILENAME', 'TOTAL_rRNA', 'UNIQUE_rRNA', 'UNIQUE_tRNA_ncRNA']

        return ncRNA_df

    def parseGeneCount(self, htseq_counts_path):
        """
            NOTE: SPECIFICALLY SET UP FOR CRYPTO
            count the gene counts that mapped either to genes (see COUNT_VARS at top of script for other features)
            :param htseq_counts_path: a path to a  _read_count.tsv file (htseq-counts output)
            :returns: a dictionary with the keys FEATURE_ALIGN_NOT_UNIQUE, TOO_LOW_AQUAL, AMBIGUOUS_FEATURE, NO_FEATURE, NOT_ALIGNED_TOTAL
        """

        sample_name = utils.pathBaseName(htseq_counts_path).replace('_read_count','')
        try:
            genotype = [self.extractInfoFromQuerySheet(sample_name, 'genotype1'), None]
            perturbation = [self.extractInfoFromQuerySheet(sample_name, 'perturbation1'), None]
        except KeyError:
            self.logger.info('Not in query sheet: %s' %htseq_counts_path)
            sys.exit('Count file passed to one of the quality assessment objects was not in the query sheet. These * should be * filtered out in the qual_assess_1 script')
        try:
            # extract genotype2 or set it to None
            genotype[1] = self.extractInfoFromQuerySheet(sample_name, 'genotype2')
            perturbation[1] = self.extractInfoFromQuerySheet(sample_name, 'perturbation2')
        except KeyError:
            self.logger.debug("%s has no genotype2 and/or perturbation2 -- may need to check script if this is expected" %sample_name)
        else:
            library_metadata_dict = {}
            # TODO: error checking on keys
            htseq_file = open(htseq_counts_path, 'r')
            htseq_file_reversed = reversed(htseq_file.readlines())

            crypto_protein_coding_count = 0
            line = next(htseq_file_reversed)
            try:
                while True:
                    line_strip_split = line.strip().split('\t')
                    if line.startswith('CKF44'):
                        # split the line, take the entry in the second column, which is the gene count, and add to crypto_protein_coding_effective_count
                        gene_count = int(line_strip_split[1])
                        crypto_protein_coding_count += gene_count
                    if not (line.startswith('CNAG') or line.startswith('CKF44')):
                        # strip newchar, split on tab
                        line = line.strip().split('\t')
                        # extract the category of metadata count (eg __alignment_not_unique --> ALIGNMENT_NOT_UNIQUE)
                        htseq_count_metadata_category = line_strip_split[0][
                                                        2:].upper()  # drop the __ in front of the category
                        # enter to htseq_count_dict
                        library_metadata_dict.setdefault(htseq_count_metadata_category, int(line[1]))
                        # iterate
                    line = next(htseq_file_reversed)
            except StopIteration:
                pass

            # error check gene count
            try:
                if crypto_protein_coding_count == 0:
                    raise ValueError('NoGeneCountsDetected')
            except ValueError:
                self.logger.info('no lines start with CKF44 -- check organism: %s' %htseq_file)
                print('No lines starting with CKF44 have gene counts')

            # rename some key/value pairs
            library_metadata_dict['NOT_ALIGNED_TOTAL'] = library_metadata_dict.pop('NOT_ALIGNED')
            library_metadata_dict['FEATURE_ALIGN_NOT_UNIQUE'] = library_metadata_dict.pop('ALIGNMENT_NOT_UNIQUE')
            library_metadata_dict['AMBIGUOUS_FEATURE'] = library_metadata_dict.pop('AMBIGUOUS')

            # add PROTEIN_CODING_COUNTED
            library_metadata_dict['PROTEIN_CODING_COUNTED'] = crypto_protein_coding_count
            # add log2cpm data -- note, this will look in the run_####_samples directory of subdir count
            log2cpm_path = os.path.join(utils.dirPath(utils.dirPath(htseq_counts_path)), '%s_log2_cpm.csv' %self.organism)
            try:
                if not os.path.isfile(log2cpm_path):
                    raise FileNotFoundError('log2cpm_pathDNE: %s' %log2cpm_path)
            except FileNotFoundError:
                msg = ' Output of log2cpm.R, which requires output of %s_raw_counts.py, ' \
                      'must be in run_####_samples directory containing subdir count. ' \
                      'This doesn\'t exist in %s' %(self.organism, sample_name)
                print(msg)
                self.logger.critical(msg)

            library_metadata_dict['NAT_LOG2CPM'] = self.extractLog2cpm('CNAG_NAT', sample_name, log2cpm_path)
            library_metadata_dict['G418_LOG2CPM'] = self.extractLog2cpm('CNAG_G418', sample_name, log2cpm_path)
            print("...extracting genotype log2cpm -- TESTING TESTING TESTING")
            if genotype[0] != 'CNAG_00000':
                library_metadata_dict['GENOTYPE1_LOG2CPM'] = self.extractLog2cpm(genotype[0].replace("CNAG", "CKF44"), sample_name, log2cpm_path)
            if genotype[1] is not None:
                library_metadata_dict['GENOTYPE2_LOG2CPM'] = self.extractLog2cpm(genotype[1].replace("CNAG", "CKF44"), sample_name, log2cpm_path)
            if perturbation[0] == "over":
                sample_medium = self.extractInfoFromQuerySheet(sample_name, 'treatment')
                sample_temperature = self.extractInfoFromQuerySheet(sample_name, 'temperature')
                sample_atmosphere = self.extractInfoFromQuerySheet(sample_name, 'atmosphere')
                sample_timepoint = self.extractInfoFromQuerySheet(sample_name, 'timePoint')
                perturbed_gene = genotype[0].replace('CNAG', 'CKF44') #TODO: THIS NEEDS TO BE CHANGED -- WHAT IF DOUBLE OVER?
                # THIS NEEDS TO BE UPDATED WITH NEW MEDIAN_LOG2CPM BY WILDTYPE REPLICATE GROUPS WHEN TREATMENT COLUMNS ARE STABLE AGAIN
                library_metadata_dict['OVEREXPRESSION_FOW'] = 0 #self.foldOverWildtype(perturbed_gene, sample_name, log2cpm_path, [sample_medium, sample_temperature, sample_atmosphere], sample_timepoint)

            htseq_file.close()
            return library_metadata_dict

    def addMarkerCoverageColumns(self, log2cpm_path):
        """

        """
        self.qual_assess_df['NAT_LOG2CPM'] = None
        self.qual_assess_df['G418_LOG2CPM'] = None
        for index, row in self.qual_assess_df.iterrows():
            # extract sample info
            sample_name = str(row['FASTQFILENAME'])
            genotype = self.extractInfoFromQuerySheet(sample_name, 'genotype1')
            sample_treatment = self.extractInfoFromQuerySheet(sample_name, 'treatment')
            sample_timepoint = self.extractInfoFromQuerySheet(sample_name, 'timepoint')
            # add NAT and G418 log2cpm
            self.qual_assess_df.loc[index, 'NAT_LOG2CPM'] = self.extractLog2cpm('CNAG_NAT', sample_name, log2cpm_path)
            self.qual_assess_df.loc[index, 'G418_LOG2CPM'] = self.extractLog2cpm('CNAG_G418', sample_name, log2cpm_path)

    def uniqueAmbiguousProteinCodingCount(self, fastq_simplename):
        """
            intersect bed file with protein coding coords with an alignment file with the htseq annotations added (in this case, grep out __ambiguous)
            :params fastq_simplename: fastq filename minus any path and extention eg /path/to/my_reads_R1.fastq.gz would be my_reads_R1
            :returns: the number of reads (lines) aligning to protein coding coordinates
        """
        # TODO: MAKE THE TRY/CATCH STATEMENT BELOW A FUNCTION AND REPLACE ALL INSTANCES WITH FUNCTION CALL

        # extract bam file from bam_file_list using fastq simple name
        try:
            bam_file = [bam_file for bam_file in self.bam_file_list if fastq_simplename in bam_file][0]
        except IndexError:
            self.logger.debug('bam file not found for %s' % fastq_simplename)  # TODO: improve this logging

        # grep out lines that are ambiguous and map to CKF44 reads (and only CKF44 -- no CNAG)
        extract_ambiguous_protein_coding_reads_cmd = 'samtools view %s| grep ambiguous | grep CKF44| grep -v CNAG| wc -l' % bam_file
        ambiguous_protein_coding_count = int(subprocess.getoutput(extract_ambiguous_protein_coding_reads_cmd))

        return ambiguous_protein_coding_count

    def calculateIntergenicCoverage(self, qual_assess_df):
        """
            calculate coverage of regions in genome between features
            Requires that the number of bases in the intergenic regions of the genome is present in OrganismData_config.ini
            in genome_files/<organism>
        """
        # create new column in qual_assess_df
        qual_assess_df['INTERGENIC_COVERAGE'] = None

        try:
            for index, row in qual_assess_df.iterrows():
                print('...Assessing intergenic coverage for %s' % row['FASTQFILENAME'])
                genotype = [self.extractInfoFromQuerySheet(row['FASTQFILENAME'], 'genotype1')]
                # set total_intergenic_bases and intergenic_region_bed by organism
                # common error message for setting total_intergenic_bases. Requires %(attribute_name, organism) both strings eg %('kn99_total_intergenic_bases', 'KN99')
                error_msg = 'No attribute %s. Check OrganismData_config.ini in %s'
                # this is a method of testing if KN99
                if genotype[0].startswith('CNAG'):
                    try:
                        intergenic_region_bed_path = os.path.join(self.genome_files, 'KN99',
                                                                  self.intergenic_region_bed)
                        total_intergenic_bases = int(self.total_intergenic_bases)
                    except KeyError:
                        kn99_error_msg = error_msg % ('kn99_total_intergenic_bases', 'KN99')
                        self.logger.critical(kn99_error_msg)
                        print(kn99_error_msg)

                # error check intergenic_region_bed_path
                try:
                    if not os.path.isfile(intergenic_region_bed_path):
                        raise FileNotFoundError('IntergenicRegionBedDoesNotExist')
                except FileNotFoundError:
                    intergenic_region_bed_error_msg = 'Intergenic region bed file does not exist at: %s' % intergenic_region_bed_path
                    self.logger.critical(intergenic_region_bed_error_msg)
                    print(intergenic_region_bed_error_msg)

                # extract intergenic_bases_covered from the bam file
                try:
                    bam_file = [bam_file for bam_file in self.bam_file_list if str(row['FASTQFILENAME']) in bam_file][0]
                except IndexError:
                    self.logger.debug(
                        'bam file not found for %s' % str(row['FASTQFILENAME']))  # TODO: improve this logging
                intergenic_bases_covered_cmd = 'samtools depth -aa -Q 10 -b %s %s | cut -f3 | grep -v 0 | wc -l' % (
                intergenic_region_bed_path, bam_file)
                num_intergenic_bases_covered = int(subprocess.getoutput(intergenic_bases_covered_cmd))
                qual_assess_df.loc[index, 'INTERGENIC_COVERAGE'] = num_intergenic_bases_covered / float(
                    total_intergenic_bases)

        except NameError:
            self.logger.critical(
                'Cannot calculate INTERGENIC COVERAGE -- total_intergenic_bases or intergenic bed file not found as attribute for organism. Check genome_files/subdirs and each OrganismData_config.ini')
        # return qual_assess_df with INTERGENIC_COVERAGE added
        return qual_assess_df

    def calculateExonicCoverage(self, bam_file, output_directory):
        """
            calculate coverage of exon regions. deposits file in output directory named utils.pathBaseName(bam_file)+'_exonic_coverage.csv'
            :param bam_file: path to bam file
            :param output_directory: path to output directory
        """
        # create df with two columns -- fastqFileName and EXONIC_COVERAGE
        exonic_df = pd.DataFrame()
        exonic_df['fastqFileName'] = [(bam_file)]
        exonic_df['EXONIC_COVERAGE'] = None

        try:
            for index, row in exonic_df.iterrows():
                bam_file = row['fastqFileName']
                print('...Assessing exonic coverage for %s' % bam_file)
                error_msg = 'No attribute %s. Check OrganismData_config.ini in %s'
                # get exon info from config file
                try:
                    exon_region_bed_path = os.path.join(self.genome_files, 'KN99', self.exon_region_bed)
                    total_exon_bases = int(self.total_exon_bases)
                except KeyError:
                    kn99_error_msg = error_msg % ('kn99_total_intergenic_bases', 'KN99')
                    self.logger.critical(kn99_error_msg)
                    print(kn99_error_msg)

                # error check intergenic_region_bed_path
                try:
                    if not os.path.isfile(exon_region_bed_path):
                        raise FileNotFoundError('IntergenicRegionBedDoesNotExist')
                except FileNotFoundError:
                    exonic_region_bed_path_error_msg = 'Intergenic region bed file does not exist at: %s' % exon_region_bed_path
                    self.logger.critical(exonic_region_bed_path_error_msg)
                    print(exonic_region_bed_path_error_msg)

                # extract exonic bases covered by at least one read
                exonic_bases_covered_cmd = 'samtools depth -aa -Q 10 -b %s %s | cut -f3 | grep -v 0 | wc -l' % (exon_region_bed_path, bam_file)
                num_exonic_bases_covered = int(subprocess.getoutput(exonic_bases_covered_cmd))

                # add to the df
                exonic_df.loc[index, 'EXONIC_COVERAGE'] = num_exonic_bases_covered / float(total_exon_bases)

        except NameError:
            self.logger.critical(
                'Cannot calculate INTERGENIC COVERAGE -- total_intergenic_bases or intergenic bed file not found as attribute for organism. Check genome_files/subdirs and each OrganismData_config.ini')
        # return qual_assess_df with INTERGENIC_COVERAGE added

        exonic_df['fastqFileName'] = exonic_df['fastqFileName'].apply(lambda x: utils.pathBaseName(x))
        # write
        output_path = os.path.join(output_directory, utils.pathBaseName(bam_file) + '_exonic_coverage.csv')
        exonic_df.to_csv(output_path, index=False)

    def perturbedCheck(self):
        """
           calculate gene coverage of genes in the 'genotype' column of the query_df that do not have suffix _over
           split genotype into two columns, genotype1, genotype2 to address double KO
        """
        # error check if bam_file_list is set
        try:
            if not hasattr(self, 'bam_file_list'):
                raise AttributeError('NoBamFileList')
        except AttributeError:
            print('QualityAssessmentObject does not have attribute bam file list')

        # log which run numbers are present
        self.logger.info('\nThe run numbers in the sheet are: %s' % self.query_df['runNumber'].unique())

        # create genotype_df from query_df[['fastqFileName' and 'genotype']]
        # get intersection of set('fastqFileName', 'genotype1', 'perturbation1', 'genotype2', 'perturbation2') set(self.query_df.columns)
        genotype_columns = list({'fastqFileName', 'genotype1', 'perturbation1', 'genotype2', 'perturbation2'}.intersection(set(self.query_df.columns)))
        # create genotype_df from intersection
        genotype_df = self.query_df[genotype_columns]
        # reduce fastqFileName to only simple basename (eg /path/to/some_fastq_R1_001.fastq.gz --> some_fastq_R1_001
        genotype_df['fastqFileName'] = genotype_df['fastqFileName'].apply(lambda x: utils.pathBaseName(x))

        # extract nat and g418 num bases from organism data
        nat_bases_in_cds = int(self.nat_cds_length)
        g418_bases_in_cds = int(self.g418_cds_length)

        # set feature over which to take percentage of reads (CDS in this case)
        feature = 'CDS'

        # create columns genotype1_coverage, genotype2_coverage, overexpression_fow (fold over wildtype), NAT_coverage, G418_coverage
        genotype_df['genotype1_coverage'] = None
        genotype_df['genotype2_coverage'] = None
        genotype_df['overexpression_fow'] = None
        genotype_df['NAT_coverage'] = None
        genotype_df['G418_coverage'] = None

        # loop over rows, calculating coverage for each genotype (testing wither genotype2 is none and perturbation is _over
        for index, row in genotype_df.iterrows():
            # simple name is like this: run_673_s_4_withindex_sequence_TGAGGTT (no containing directories, no extention)
            fastq_simple_name = utils.pathBaseName(row['fastqFileName'])
            # get bam files which correspond to query_df fastqFileNames
            try:
                bam_file = [bam_file for bam_file in self.bam_file_list if fastq_simple_name in bam_file][0]
            except IndexError:
                self.logger.info('bam file not found for %s' % fastq_simple_name)  # TODO: improve this logging
                continue

            # calculate marker coverages
            print('...calculating NAT coverage for %s' % fastq_simple_name)
            genotype_df.loc[index, 'NAT_coverage'] = self.calculatePercentFeatureCoverage(feature, 'CNAG_NAT',
                                                                                          self.annotation_file,
                                                                                          bam_file, nat_bases_in_cds)
            print('...calculating G418 coverage for %s' % fastq_simple_name)
            genotype_df.loc[index, 'G418_coverage'] = self.calculatePercentFeatureCoverage(feature, 'CNAG_G418',
                                                                                           self.annotation_file,
                                                                                           bam_file, g418_bases_in_cds)

            # if deletion, calculate coverage. Currently only set to check genotype1. assumes both are deletions if perturbation1 == 'deletion'
            if row['perturbation1'] == "deletion" or row['perturbation1'] == "geneSwap":
                try:
                    # extract genotype1 TODO: JUST CASE EVERYTHING TO UPPER EARLIER
                    genotype = [self.extractInfoFromQuerySheet(row['fastqFileName'], 'genotype1'), None]
                    #genotype = [list(self.query_df[self.query_df['fastqFileName'].str.contains(row['FASTQFILENAME'] + '.fastq.gz')]['genotype1'])[0]]
                except ValueError:
                    self.logger.info('genotype cannot be extracted with the fastq filename in this row. Note: if there are null entries in the column fastqFileNames, this is the cause. those need to be remedied or removed in order for this to work: %s' %row)
                try:
                    # extract genotype2 or set it to None
                    genotype[1] = self.extractInfoFromQuerySheet(row['fastqFileName'], 'genotype2')
                except KeyError:
                    self.logger.debug("sample: %s does not have genotype2" %row['fastqFileName'])
                # determine which genome to use -- if CNAG, use KN99
                if not genotype[0].startswith('CNAG'):
                    raise ValueError('%sNotRecognizedCryptoGenotype' % genotype[0])
                # replace CNAG with CKF44 (in past version of pipeline, KN99 genes were labelled with H99 names. NCBI required change to CKF. Numbering/order is same -- just need to switch CNAG to CKF44)
                genotype[0] = genotype[0].replace('CNAG', 'CKF44')
                if genotype[1] not in [None, 'nan'] and genotype[1].startswith('CNAG'):
                    genotype[1] = genotype[1].replace('CNAG', 'CKF44')
                print('...checking coverage of %s in %s' % (genotype, fastq_simple_name))
                genotype_df.loc[index, 'genotype1_coverage'] = self.calculatePercentFeatureCoverage(feature, genotype[0],
                                                                                                    self.annotation_file,
                                                                                                    bam_file)
                # do the same for genotype2 if it exists
                if genotype[1] not in [None, 'nan']:
                    genotype_df.loc[index, 'genotype2_coverage'] = self.calculatePercentFeatureCoverage(feature, genotype[1],
                                                                                                        self.annotation_file,
                                                                                                        bam_file)
        # return genotype check
        genotype_df.columns = [column_name.upper() for column_name in genotype_df.columns]
        return genotype_df[['FASTQFILENAME', 'GENOTYPE1_COVERAGE', 'GENOTYPE2_COVERAGE', 'NAT_COVERAGE', 'G418_COVERAGE']]
