import configparser
from rnaseq_tools import utils
from rnaseq_tools.CryptoQualityAssessmentObject import CryptoQualityAssessmentObject
import numpy as np

class CryptoQualAssessAuditObject(CryptoQualityAssessmentObject):

    def __init__(self, expected_attributes=None, **kwargs):
        # add expected attributes to super._attributes
        self._add_expected_attributes = []
        # This is a method of adding expected attributes to StandardData from StandardData children
        if isinstance(expected_attributes, list):
            self._add_expected_attributes.extend(expected_attributes)
        # initialize Standard data with the extended _attributes
        # recall that this will check for and/or create the directory structure found at
        super(CryptoQualAssessAuditObject, self).__init__(self._add_expected_attributes, **kwargs)
        # overwrite super.self_type with object type of child (this object)
        self.self_type = 'CryptoQualAssessAuditObject'
        # create logger
        self.logger = utils.createStandardObjectChildLogger(self, __name__)

        # extract threshold/status from config file #TODO: move to constructor
        qual_assess_config = configparser.ConfigParser()
        qual_assess_config.read(self.config_file)
        qual_assess_1_dict = qual_assess_config['KN99QualityAssessOne']

        # extract thresholds #TODO: CLEAN UP TO DICTIONARY, AUTOMATICALLY EXTRACT
        self.protein_coding_total_threshold = int(qual_assess_1_dict['PROTEIN_CODING_TOTAL_THRESHOLD'])
        self.not_aligned_total_percent_threshold = float(qual_assess_1_dict['NOT_ALIGNED_TOTAL_PERCENT_THRESHOLD'])
        self.perturbed_coverage_threshold = float(qual_assess_1_dict['PERTURBED_COVERAGE_THRESHOLD'])
        self.nat_expected_coverage_threshold = float(qual_assess_1_dict['NAT_EXPECTED_COVERAGE_THRESHOLD'])
        self.nat_expected_log2cpm_threshold = float(qual_assess_1_dict['NAT_EXPECTED_LOG2CPM_THRESHOLD'])
        self.nat_unexpected_coverage_threshold = float(qual_assess_1_dict['NAT_UNEXPECTED_COVERAGE_THRESHOLD'])
        self.nat_unexpected_log2cpm_threshold = float(qual_assess_1_dict['NAT_UNEXPECTED_LOG2CPM_THRESHOLD'])
        self.g418_log2cpm_threshold = float(qual_assess_1_dict['G418_LOG2CPM_THRESHOLD'])
        self.overexpression_fow_threshold = float(qual_assess_1_dict['OVEREXPRESSION_FOW_THRESHOLD'])

        # extract status
        self.protein_coding_total_bit_status = int(qual_assess_1_dict['PROTEIN_CODING_TOTAL_STATUS'])
        self.not_aligned_total_percent_bit_status = int(qual_assess_1_dict['NOT_ALIGNED_TOTAL_PERCENT_STATUS'])
        self.perturbed_coverage_bit_status = int(qual_assess_1_dict['PERTURBED_COVERAGE_STATUS'])
        self.nat_expected_marker_status = int(qual_assess_1_dict['NAT_EXPECTED_MARKER_STATUS'])
        self.nat_unexpected_marker_status = int(qual_assess_1_dict['NAT_UNEXPECTED_MARKER_STATUS'])
        self.g418_expected_marker_status = int(qual_assess_1_dict['G418_EXPECTED_MARKER_STATUS'])
        self.g418_unexpected_marker_status = int(qual_assess_1_dict['G418_UNEXPECTED_MARKER_STATUS'])
        self.overexpression_fow_status = int(qual_assess_1_dict['OVEREXPRESSION_FOW_STATUS'])
        self.no_metadata_marker_status = int(qual_assess_1_dict['NO_METADATA_MARKER_STATUS'])

        self.auditQualAssessDataframe()

    def auditQualAssessDataframe(self):
        """
            use rnaseq_pipeline/config/quality_assess_config.ini entries to add status, auto_audit columns
            :params qual_assess_df: a quality_assess_df created by qual_assess_1 or a path to one
            :params bam_file_list: the bam_file_list used in qual_assess_1 -- TODO: make this a list of log2cpm files
            :returns: qual_assess_df with added status and auto_audit columns
        """
        # instantiate a list to hold the status (sum of the status' values when a threshold is failed)
        status_column_list = []

        # loop over the rows of qual_assess_df
        for index, row in self.qual_assess_df.iterrows():
            # extract sample information
            fastq_simple_name = str(row['FASTQFILENAME'])
            # extract genotype -- if double KO, there will be two genotypes listed. this is a point of weakness in generalizing genetic perturbations
            try:
                genotype1 = self.extractInfoFromQuerySheet(fastq_simple_name, 'genotype1')
                genotype2 = self.extractInfoFromQuerySheet(fastq_simple_name, 'genotype2')
                if genotype2 in ["na", "NA", "nan", "NaN", "Nan", None]:
                    raise KeyError("genotype2 not present")
            except KeyError:
                genotype2 = None
            genotype = [genotype1, genotype2]
            # marker_1 and marker_2 will be either NAT or G418, depending on what is expected, or NA
            marker_1 = self.extractInfoFromQuerySheet(fastq_simple_name, 'marker1')
            marker_2 = self.extractInfoFromQuerySheet(fastq_simple_name, 'marker2')

            # extract quality_assessment_metrics
            library_protein_coding_total = int(row['PROTEIN_CODING_TOTAL'])
            not_aligned_total_percent = float(row['NOT_ALIGNED_TOTAL_PERCENT'])
            if row['GENOTYPE1_COVERAGE'] is not None:
                library_genotype1_coverage = float(row['GENOTYPE1_COVERAGE'])
            else:
                library_genotype1_coverage = -1
            if row['GENOTYPE2_COVERAGE'] is not None:
                library_genotype2_coverage = float(row['GENOTYPE2_COVERAGE'])
            else:
                library_genotype2_coverage = -1
            if row['OVEREXPRESSION_FOW'] is not None:
                overexpression_fow = float(row['OVEREXPRESSION_FOW'])
            else:
                overexpression_fow = 1000
            if row['NAT_COVERAGE'] is not None:
                nat_coverage = float(row['NAT_COVERAGE'])
            else:
                nat_coverage = -1
            nat_log2cpm = float(row['NAT_LOG2CPM'])
            g418_log2cpm = float(row['G418_LOG2CPM']) # g418 evaluated on log2cpm alone

            # set status_total to 0
            status_total = 0

            # check values against thresholds, add status to flag
            if library_protein_coding_total < self.protein_coding_total_threshold:
                status_total += self.protein_coding_total_bit_status
            if not_aligned_total_percent > self.not_aligned_total_percent_threshold:
                status_total += self.not_aligned_total_percent_bit_status
            # if the overexpression_flag is set, eval based on expression
            if overexpression_fow < self.overexpression_fow_threshold:
                status_total += self.overexpression_fow_status
            # else, evaluate KO based on coverage
            else:
                if library_genotype1_coverage > self.perturbed_coverage_threshold or library_genotype2_coverage > self.perturbed_coverage_threshold:
                    status_total += self.perturbed_coverage_bit_status

            # test wildtypes for marker coverage and expression
            if genotype[0] == 'CNAG_00000':
                if nat_coverage > self.nat_unexpected_coverage_threshold and nat_log2cpm > self.nat_unexpected_log2cpm_threshold:
                    status_total += self.nat_unexpected_marker_status
                if g418_log2cpm > self.g418_log2cpm_threshold:
                    status_total += self.g418_unexpected_marker_status
            # if a perturbed sample, first check that marker information is present and flag it if it is not
            else:
                if marker_1 == 'nan' or marker_1 is None or (
                        genotype[1] != None and (marker_2 == 'nan' or marker_2 == 'none')):
                    status_total += self.no_metadata_marker_status
                # if perturbed, and marker information is present, test the markers
                else:
                    if marker_1 == 'NAT':
                        if nat_coverage < self.nat_expected_coverage_threshold or nat_log2cpm < self.nat_expected_log2cpm_threshold:
                            status_total += self.nat_expected_marker_status
                        if g418_log2cpm > self.g418_log2cpm_threshold and genotype[1] != None:
                            status_total += self.g418_unexpected_marker_status
                    elif marker_1 == 'G418':
                        if g418_log2cpm < self.g418_log2cpm_threshold:
                            status_total += self.g418_expected_marker_status
                        # test coverage if single KO
                        if (nat_coverage > self.nat_unexpected_coverage_threshold and nat_log2cpm > self.nat_unexpected_log2cpm_threshold) and genotype[1] == None:
                            status_total += self.nat_unexpected_marker_status
                    # TODO: THIS NEEDS TO BE FIXED FOR THE INSTANCE IN WHICH A SINGLE MARKER IS NOTED WITH BOTH MARKERS B/C OF MANUAL ENTRY (SEE STRAINS 2274 AND 1351)
                    # TEST THIS WITH THE STRAINS ABOVE
                    # this is a double KO,  or there is supposed to be a double ko according to the filled marker_2
                    if genotype[1] != None or marker_2 != 'nan':  # note: unentered 2nd markers for double KO should be caught in the if statement above
                        if marker_2 == 'NAT':
                            if marker_1 == 'NAT':
                                self.logger.critical('%s has two NAT markers in the metadata' % fastq_simple_name)
                            if nat_coverage < self.nat_expected_coverage_threshold or nat_log2cpm < self.nat_expected_log2cpm_threshold:
                                status_total += self.nat_expected_marker_status
                            if g418_log2cpm > self.g418_log2cpm_threshold and genotype[1] == None:
                                status_total += self.g418_unexpected_marker_status
                        elif marker_2 == 'G418':
                            if marker_1 == 'G418':
                                self.logger.critical('%s has two G418 markers in the metadata' % fastq_simple_name)
                            if g418_log2cpm < self.g418_log2cpm_threshold and genotype[1] == None:
                                status_total += self.g418_expected_marker_status
                            # if nat_coverage > unexpected_nat_coverage_threshold or nat_log2cpm > unexpected_nat_log2cpm_threshold:
                            #     status_total += nat_unexpected_marker_status

            status_column_list.append(status_total)

        self.qual_assess_df['STATUS'] = status_column_list
        self.qual_assess_df['AUTO_AUDIT'] = np.where(self.qual_assess_df.STATUS > 0, 1, 0)
        self.qual_assess_df['STATUS_DECOMP'] = None
        # add status decomposition
        for index, row in self.qual_assess_df.iterrows():
            status = int(row['STATUS'])
            self.qual_assess_df.loc[index, 'STATUS_DECOMP'] = str(utils.decomposeStatus2Bit(status))

        return self.qual_assess_df
