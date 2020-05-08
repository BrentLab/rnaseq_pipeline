from rnaseq_tools.StandardDataObject import StandardData
from rnaseq_tools.DatabaseObject import DatabaseObject


class QualityAssessmentObject(StandardData):

    def __init__(self, expected_attributes=None, **kwargs):
        # add expected attributes to super._attributes
        self._add_expected_attributes = ['quality_assessment_filename']
        # TODO: This is a messy and repetitive way of adding expected attributes from children of OrganismData to add to StandardData
        if isinstance(expected_attributes, list):
            self._add_expected_attributes.extend(expected_attributes)
        # initialize Standard data with the extended _attributes
        # recall that this will check for and/or create the directory structure found at
        super(QualityAssessmentObject, self).__init__(self._add_expected_attributes, **kwargs)
        # overwrite super.self_type with object type of child (this object)
        self.self_type = 'QualityAssessmentObject'

        # set optional kwarg arguments
        try:
            self.standardized_database_df = kwargs['standardized_database_df']
        except KeyError:
            pass
        try:
            self.genotype_list = kwargs['genotype_list']
        except KeyError:
            self.genotype_list = []
        try:
            self.ko_gene_list = kwargs['ko_gene_list']
        except KeyError:
            self.ko_gene_list = [] # eg ['CNAG_01020' , ['CNAG_39392','CNAG_48382'], 'CNAG_23421'] where the center item is a double ko
        try:
            self.overexpress_gene_list = kwargs['overexpress_gene_list']
        except KeyError:
            self.overexpress_gene_list = [] # expecting no nested lists in this

    def setCryptoGenotypeList(self):
        """ # TODO: Move this to DatabaseObject
            from database_df extract crypto genotype_list (no wildtype CNAG_00000)
        """
        try:
            self.genotype_list = self.standardized_database_df.GENOTYPE.unique()
        except AttributeError:
            print('QualityAssessObject has no standardized_database_df attribute')

    def parseCryptoGenotypeList(self):
        """
            parse a list of genotypes extracted from database_df into ko_gene_list (list of knockout genes) and
            overexpress_gene_list (overexpression genes)
        """
        if len(self.genotype_list) == 0:
            raise IndexError('GenotypeListEmpty')
        for genotype in self.genotype_list:
            if not genotype == 'CNAG_00000':
                if '_over' in genotype:
                    self.overexpress_gene_list.append(genotype.split('_over')[0])
                else:
                    if '.' in genotype:
                        self.ko_gene_list.append(genotype.split('.'))
                    else:
                        self.ko_gene_list.append(genotype)


    def cryptoPerturbationGenotypeCheck(self):
        # check if necessary paths are entered as attributes
        if not hasattr(self, 'query_sheet_path'):
            raise AttributeError('NoQuerySheetPath')
        if not hasattr(self, 'log2_cpm_path'):
            raise AttributeError('NoLog2CpmPath')
        if not hasattr(self, 'experiment_columns'):
            raise AttributeError('NoExperimentColumns')

        # if standardized_database_df is not yet set, do so
        # if align_count_path is present, this is for quality_assess_1 and the prefix should be set as such in column COUNTFILENAME
        if hasattr(self, 'align_count_path') and not hasattr(self, 'standardized_database_df'):
            self.standardized_database_df = DatabaseObject.standardizeDatabaseDataframe(self.query_sheet_path,
                                                                                     self.align_count_path)
        # else, no prefix in COUNTFILENAME (see DatabaseObject.StandardizeDatabaseDataframe)
        elif not hasattr(self, 'standardized_database_df'):
            self.standardized_database_df = DatabaseObject.standardizeDatabaseDataframe(self.query_sheet_path)

        for gene in self.ko_gene_list:
            pass
            # extract log2cpm expression and treatment_timepoint
            # extract treatment_timepoint
            # get quantile from wildtype dataframe for gene in given treatment_timepoint
            # if gene expression is greater than 20% treatment_timempoint, flag, create diagnostic dataframe and browser shot

        for gene in self.overexpress_gene_list:
            pass
            # extract log2cpm expression and treatment_timepoint
            # extract treatment_timepoint
            # get quantile from wildtype dataframe for gene in given treatment_timepoint
            # if gene expression is less than 99% wildtype, flag and create diagnostic dataframe and browser shot

    def updateStatusColumn(self):
        pass
