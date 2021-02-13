import os
import re
import pandas as pd
from glob import glob
import numpy as np
from itertools import combinations, product
import yaml
import sys
import subprocess
import configparser
import time
import logging
import logging.config
import math
import requests


def getRunNumber(fastq_path):
    """
        extract run number from -f input. this *should be* to a file called sequence/run_####_samples
        :param fastq_path
        :returns: the run number as string (leading zero will be retained)
    """
    try:
        regex = r"run_(\d*)"
        run_num = re.search(regex, fastq_path).group(1)
    except AttributeError:
        raise AttributeError('NoRunNumber')  # raise this, handle in the program that uses utils.getRunNumber
    else:
        return run_num


def decomposeStatus2Bit(status):
    """
        Decompose a value in the column of the output of quality_assess_2 to bit.
        eg 18 = 2 + 16 or the powers of 2 [1.0,4.0]. See templates/qc_config.yaml
        :param status: an int from the status column of the .csv output of quality_assess_2
        :returns: a list of powers of 2 representing the bit
    """
    if status == 0:
        return None
    decomp = []
    for i in np.arange(np.floor(np.log2(status)), -1, -1):
        if status - 2 ** i >= 0:
            decomp.append(i)
            status -= 2 ** i
    return decomp


def makeCombinations(lst):
    """
        Make all possible replicate combinations of the inputted list
        :param lst: a list of items
        :returns: a list of combinations of the inputted list
    """
    if len(lst) < 2:
        return [list]
    combos = []
    for i in range(len(lst), 1, -1):
        for s in combinations(lst, i):
            combos.append(sorted(s))
    return combos


# def makeListProduct(lst):  #  TODO: IF COMMETING THIS OUT DOESN'T BERAK ANYTHING, REMOVE
#     """
#     Make all possible value combination, with each value coming from a each list. Support up to 10 lists.
#     :param lst: a list of items
#     :returns: all value combinations of the inputted list
#     """
#     if len(lst) == 0:
#         return None
#     elif len(lst) == 1:
#         return list(product(lst[0]))
#     elif len(lst) == 2:
#         return list(product(lst[0], lst[1]))
#     elif len(lst) == 3:
#         return list(product(lst[0], lst[1], lst[2]))
#     elif len(lst) == 4:
#         return list(product(lst[0], lst[1], lst[2], lst[3]))
#     elif len(lst) == 5:
#         return list(product(lst[0], lst[1], lst[2], lst[3], lst[4]))
#     else:
#         sys.exit('ERROR: The number of contrasts is greater than 10. This length beyond what I can handle')


def loadConfig(json_file):
    """
        Load configuration file (JSON) for QC thresholding and scoring
        :param json_file: the config for the rnaseq_pipeline, in particular, is found in templates/qc_config.yaml
        :returns: the read-in json object
    """
    with open(json_file) as json_data:
        d = yaml.safe_load(json_data)
    return d


def mkdirp(path_to_directory):
    """
        Function to create a directory. Equivalent to mkdir -p in bash.
        Will create directory if the path does not already exist.
        :param path_to_directory: path to a directory to be created if DNE already
    """
    if not os.path.exists(path_to_directory):
        try:
            os.makedirs(path_to_directory)
        except OSError:
            print("Directory {} cannot be created. See utils.mkdirp function call "
                  "in the script from which this error occurred".format(path_to_directory))


def addForwardSlash(path):
    """
        add forward slash to path
        :param path: a path, typically to the lowest subdir prior to adding a file
        :returns: the path with a trailing slash added
    """
    if not path.endswith('/'):
        path = path + '/'
    return path


def removeForwardSlash(path):
    """
        removes forward slash from path
        :param path: filepath
        :returns: path without final forward slash
    """
    if path.endswith('/'):
        path = path[:-1]
    return path


def checkCSV(file_path):
    """
        test whether a given file is a .csv or something else
        :param file_path: a file that ends with .csv, .xlsx, .tsv, etc.
        :returns: True if .csv, false otherwise
    """
    if not os.path.isfile(file_path):
        raise FileNotFoundError('PathToColumnarDataDNE')
    else:
        # test whether a given file is a .csv or .xlsx
        if re.search('\.csv', file_path):
            return True
        else:
            return False


def checkTSV(file_path):
    """
        test whether a given file is a .tsv
        :param file_path: a file path
        :returns: True if .tsv, false otherwise
        :raises: FileNotFoundError
    """
    if not os.path.isfile(file_path):
        raise FileNotFoundError('PathToColumnarDataDNE')
    else:
        if re.search('\.tsv', file_path):
            return True
        else:
            return False


def checkExcel(file_path):
    """
        test whether a given file is a .xlsx
        :param file_path: a file path
        :returns: True if .tsv, false otherwise
    """
    if not os.path.isfile(file_path):
        raise FileNotFoundError('PathToColumnarDataDNE')
    else:
        # test whether a given file is a .xlsx
        if re.search('\.xlsx', file_path):
            return True
        else:
            return False


def readInDataframe(path_to_csv_tsv_or_excel):
    """
        read in .csv, .tsv or .xlsx
        :param path_to_csv_tsv_or_excel: path to a .csv, .tsv .xlsx
        :returns: a pandas dataframe
    """
    if not (path_to_csv_tsv_or_excel.endswith('csv') or
            path_to_csv_tsv_or_excel.endswith('tsv') or
            path_to_csv_tsv_or_excel.endswith('xlsx')):
        raise ValueError('UnrecognizedFileExtension')
    try:
        if checkCSV(path_to_csv_tsv_or_excel):
            return pd.read_csv(path_to_csv_tsv_or_excel)
        elif checkTSV(path_to_csv_tsv_or_excel):
            return pd.read_csv(path_to_csv_tsv_or_excel, sep='\t')
        elif checkExcel(path_to_csv_tsv_or_excel):
            return pd.read_excel(path_to_csv_tsv_or_excel)
    except FileNotFoundError:
        raise FileNotFoundError('PathToColumnarDataDNE')


def fileBaseName(file_name):
    """
        strip everything after first . in path (may be too extreme -- make sure there are no 'necessary' . in filepaths. However, please try to avoid putting a . in a filepath other than the extension(s)
        :param file_name: the base filename of a file (not a path!)
        :returns: the filename stripped of all extensions
    """
    # https://stackoverflow.com/a/46811091
    if '.' in file_name:
        separator_index = file_name.index('.')
        base_name = file_name[:separator_index]
        return base_name
    else:
        return file_name


def pathBaseName(path):
    """
        This function will work for fileBaseName and also where you want the basename of a file from a path of any length
        strips directory structure and extensions. eg /path/to/file.fastq.gz --> file
        :param path: a filepath that includes a directory structure and file (see description)
        :returns: the basename of the file stripped of the directory structure and extensions
    """
    # This gets the basename of a given path, and then strips all file extensions (even if multiple). see fileBaseName
    file_name = os.path.basename(path)
    return fileBaseName(file_name)


def dirName(path):
    """
        get the directory name of the directory one level up from the end of the path provided.
        eg /path/to/file.fastq.gz returns 'to'
        :param path: a path to a directory or file
        :returns: the name of the directory one level up from the final level of path. see description for example
    """
    path = removeForwardSlash(path)
    if path == '.':
        raise Exception('Must pass a path, not \\., in order to extract directory name')
    directory_path_split = os.path.split(path)
    # search second group for a period (this is not very specific -- a non file ending period would also return true)
    if re.search('\\.', directory_path_split[1]):
        dir_name = os.path.basename(directory_path_split[0])
    else:
        dir_name = os.path.basename(path)
    return dir_name


def dirPath(path):
    """
        return the path (first half of result of os.path.split() to the directory one up from file
        :param path: a path to a directory or file
        :returns: the name of the directory one level up from the final level of path. see description for example
        :raises: NotADirectoryError
    """
    path_to_directory_one_up = os.path.split(path)[0]

    if not os.path.isdir(path_to_directory_one_up):
        raise NotADirectoryError

    return path_to_directory_one_up


def extractTopmostFiles(path_to_directory):
    """
        extract a list of all files in a given directory
        :param path_to_directory: path to the directory
        :returns: a list of files in topmost level of directory (not recursive)
    """
    return glob(os.path.join(path_to_directory, '*'))


def countsPerMillion(raw_count_path, output_FULL_path):
    """ TODO: re-write this with python subprocess to control input/output of R script
        submit raw_count_path to log2_cpm.R (in tools/)
        :param raw_count_path: path to output of raw_count.py
        :param output_FULL_path: the full path (including the file and extension) of the output of log2_cpm.R.
        eg <experiment_name>_log2_cpm.csv
    """
    cmd = 'log2_cpm.R -r %s -o %s' % (raw_count_path, output_FULL_path)
    executeSubProcess(cmd)


def executeSubProcess(cmd):
    """ TODO: re-do this with package subprocess and store rather than write. write to logger
        executes command, sys.exit with message if the subprocess fails
        :param cmd: the full command to be run
    """
    exit_status = subprocess.call(cmd, shell=True)
    if exit_status == 1:
        raise IOError('%s failed to execute. check the code' % cmd)


def configure(object_instance, config_file, config_header, prefix=''):
    """
        reads and sets the attributes in a config_file.ini in type output by configparser.
        :param object_instance: an object to be configured
        :param config_file: a .ini
        :param config_header: the [header] in the .ini file to read (config format created by configparser).
        :param prefix:
        The function will loop through these key, value pairs and set attributes accordingly
    """
    # read config file
    config = configparser.ConfigParser()
    config.read(config_file)
    # set attributes for StandardData
    for key, value in config[config_header].items():
        setattr(object_instance, key, os.path.join(prefix,
                                                   value))  # by default, values are read in as strings. Currently, all filepaths, so this is good


def submitSbatch(sbatch_path, email=None):
    """
        submit an sbatch script (see sbatch script examples in templates)
        :param sbatch_path: path to an sbatch file
        :param  email: user email, default None
    """
    print('...submitting sbatch job')
    # Submit sbatch job
    if email is None:
        cmd = "sbatch {}".format(sbatch_path)
        executeSubProcess(cmd)
    else:
        cmd = "sbatch --mail-type=END,FAIL --mail-user={0} {1}".format(email, sbatch_path)
        executeSubProcess(cmd)


def yearMonthDay():
    """
        :returns: the year-month-day as 20200403
    """
    return time.strftime("%Y%m%d")


def hourMinuteSecond():
    """
        :returns: hour-minute-second as 135301
    """
    hour_minute_second = str(time.strftime("%H%M%S"))
    return hour_minute_second


def softLinkAndSetAttr(object_instance, list_of_target_directories, source_path, intended_target_dir):
    """
        creates soft links and stores the path as attributes of object_instance
        :param object_instance: an instance of a given object (StandardData in this case)
        :param list_of_target_directories: a list of directories that exists in origin_dir_path
        :param source_path: path to where the directories exist (actual data)
        :param intended_target_dir: the place where you'd like the directories soft linked
        :raises: FileNotFoundError
    """
    if not os.path.exists(source_path):
        raise FileNotFoundError('SourceDoesNotExist')
    if not os.path.exists(intended_target_dir):
        raise FileNotFoundError('TargetDirectoryOfSoftLinkDoesNotExist')

    # loop through directories in list_of_dirs
    for target_directory_name in list_of_target_directories:
        # store the path (either it exists or it will after the if not block) to the directory in intended_dir_path)
        target_path = os.path.join(intended_target_dir, target_directory_name)
        # check if directory exists in the intended_dir_path
        if not os.path.exists(target_path):
            # if it does not, soft link ln -s origin_dir_path/directory to intended_dir_path/directory
            cmd = 'ln -s {}/{} {}'.format(source_path, target_directory_name, target_path)
            executeSubProcess(cmd)
        # set attribute named directory (from for loop above) that points towards variable path which stores intended_dir_path/directory
        setattr(object_instance, target_directory_name, target_path)


def setAttributes(sd_object, input_dict):
    """
        a makeshift method to mimic a constructor.
        (a check for entering 'query' instead of 'query_sheet_path'
        :param sd_object: a instance of a StandardData object in which you wish to set attributes
        :param input_dict: kwargs entered upon instantiation
    """
    for key, value in input_dict.items():
        setattr(sd_object, key, value)


def userInputCorrectPath(message, my_object, attribute, **kwargs):
    """
        asks user to correct input path
        usage: object = userInputCorrectPath(...)
        :param message: message to print to user
        :param my_object: object whose attribute you wish to correct
        :param attribute: attribute to be corrected
        :param kwargs: arbitrary list of keyword arguments. none currently handled. intended use is for logger
        :returns: my_object with attribute set to new value
    """
    # method to prompt user to correct path
    print(message)

    # args and the if block below are for testing purposes
    if kwargs:
        new_attr_value = kwargs[0][0]
    else:
        new_attr_value = input()

    return setattr(my_object, attribute, new_attr_value)


def userInputCorrectAttributeName(message, my_object, old_attribute_name, *args):
    """
        ask user to correct attribute name in object
        usage: object = userInputCorrectAttributeName(...)
        :param message: message to print to user
        :param my_object: object whose attribute you wish to correct the name of
        :param old_attribute_name: name of the attribute you wish to change
        :param kwargs: arbitrary list of keyword arguments. none currently handled. intended use is for logger
        :returns: object with renamed attribute
    """
    print(message + " ")
    # this if statement and the args parameter are for testing purposes
    if args:
        new_attribute_name = args[0][0]
    else:
        new_attribute_name = input()

    setattr(object, new_attribute_name, getattr(object, old_attribute_name))
    delattr(object, old_attribute_name)
    # to use this, user will need to object = userInputCorrectAttributeName(...)
    return my_object


def createLogger(log_file_path, logger_name, logger_level, logging_conf=None):
    """
    create logger in filemode append and format name-levelname-message with package/module __name__ (best practice from logger tutorial)
    :param log_file_path: name of the file in which to log.
    :param logger_name: __name__ is recommended as a best practice in logger.config eg you can call this like so: createLogger(<your_filename>, __name__)
                     (__name__ is a special variable in python)
    :param logger_level: at what level to publish to log file see logger docs
    :param logging_conf: path to logging configuration file
    :returns: an instance of the configured logger
    """
    try:
        log_directory_path = dirPath(log_file_path)
    except NotADirectoryError:
        print('log directory %s not found and therefore can\'t create log_file here.' % log_directory_path)
    # a config file is passed, load it
    if logging_conf:
        logging.config.fileConfig(logging_conf)  # should include at least what is below
    # if it is not, configure as follows
    else:
        # create log for the year-month-day
        logging.basicConfig(
            filename='%s' % log_file_path,
            filemode='a',
            format='%(name)s-%(levelname)s-%(asctime)s-%(message)s',
            datefmt='%I:%M:%S %p',  # set 'datefmt' to hour-minute-second AM/PM
            level=logger_level
        )
    # return an instance of the configured logger
    return logging.getLogger(logger_name)


def createStandardObjectChildLogger(StandardDataObjectChild, name):
    """
        create logger for StandardDataObjectChild
        :param StandardDataObjectChild: a child of StandardData. If this is not called, the child's logger will be named "StandardDataObject"
        :param name: probably passed with __name__. The name of the logger to display in the log file.
        :raises: NotADirectoryError if logger_directory_path does not exist
    """
    logger_directory_path = dirPath(StandardDataObjectChild.log_file_path)
    if os.path.isdir(logger_directory_path):
        try:
            logger = createLogger(StandardDataObjectChild.log_file_path, name, StandardDataObjectChild.logger_level)
        except AttributeError:
            print('StandardDataObject constructor not yet run -- logger_level attribute not set.\n'
                  'Must run StandardDataObject constructor first (see another StandardDataObject constructor)')
    else:
        raise NotADirectoryError('LogDirectoryDoesNotExist')

    return logger


def createStdOutLogger(**kwargs):
    """
        create a logger that writes to stdout
        :param kwargs: optional key word arguments, right now only set to handle name
        :credit: https://stackoverflow.com/a/14058475/9708266
    """
    # create logger
    try:
        logger = logging.getLogger(kwargs['name'])
    except KeyError:
        logger = logging.getLogger()
    # configure logger
    logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%I:%M:%S %p')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    return logger


def getFileListFromDirectory(dir_path,
                             list_of_file_suffixes_to_extract):  # TODO: currently set up to work specifically for align_counts. re-write usage and make more flexible
    """
    write fastq filepaths in a list stored as a .txt. Used in slurm job script
    :param dir_path: path to a diretory with files you wish to extract
    :param list_of_file_suffixes_to_extract: a list of suffixes
    :returns: list of files from directory with a suffix matching one in list_of_file_suffixes_to_extract
    """
    try:
        if not isinstance(list_of_file_suffixes_to_extract, list):
            raise TypeError('NotAList')
    except TypeError:
        print(
            'You must pass a list, even if it is just one item, to getFileListFromDirectory for list_of_file_suffixes_to_extract')
    else:
        file_paths = []
        for suffix in list_of_file_suffixes_to_extract:
            file_paths += glob(dir_path + "/*." + suffix)
        return file_paths


def extractFiles(containing_directory, file_pattern):
    """
        extract a list of files, expected to be non-empty, from containing_directory. NOTE: this is recursive, so be
        specific in your pattern or risk returning everyting in a given containing_directory that matches
        :params containing_directory: a path to the directory expected to contain the files with a given pattern
        :params file_pattern: a pattern to search the containing directory for
        :raises: UserWarning if file_pattern starts with astrisk, NotADirectoryError if containig_directory dne, FileNotFoundError if list is length zero
        :returns: a list of files matching the file_pattern in the containing_directory
    """
    try:
        if file_pattern.startswith('*'):
            raise UserWarning('DoNotStartPatternWithAstrisk')
    except UserWarning:
        print(
            'When using utils.extractFiles, the argument file_pattern should not start with an astrisk. That is added in the method.')

    try:
        if os.path.isdir(containing_directory):
            file_list = glob('%s/**/*%s' % (containing_directory, file_pattern), recursive=True)
        else:
            raise NotADirectoryError('AlignDirectoryDoesNotExist')
        if len(file_list) == 0:
            raise FileNotFoundError('NoFileFoundMatchingPatternInContainingDirectory')
    except NotADirectoryError:
        print('Directory %s does not exist' % containing_directory)
    except FileNotFoundError:
        print(
            'No Files found matching the pattern %s in containing directory %s' % (file_pattern, containing_directory))
    else:
        return file_list


def logit2probability(log_odds):
    """
        convert log odds to probability
        params log_odds: a log odds score
        :returns: probability converted from input
        credit: https://sebastiansauer.github.io/convert_logit2prob/
    """
    # convert log odds to odds
    odds = math.exp(log_odds)
    # return probability
    return odds / float(1 + odds)


def twoParameterGlmTemplate(intercept, coefficient_1, coefficient_2):
    """
        factory function to create a two parameter glm formula.
        usage: to create the formula, first call twoParameterGlmTemplate:
            my_glm = twoParameterGlmTemplate(1,2,3)
        then, you can use my_glm to return x if x = 1 + 2*x_1 + 3*x_2 where x_1 and x_2 are inputs to my_glm
            for example, my_glm(1,1) == 6
        :params intercept: intercept of the model you wish to create
        :params coefficient_1: first coefficient of the model you wish to create
        :params coefficient_2: second coefficient of the model you wish to create
        :returns: a two parameter glm with the intercept and two coefficients
    """

    # define a specific instance of a two parameter glm
    def twoParameterGlm(x_1, x_2):
        if not testNumeric(x_1):
            raise ValueError('VariablesNonNumeric %s' % x_1)
        if not testNumeric(x_2):
            raise ValueError('VariablesNonNumeric %s' % x_2)
        x_1 = float(x_1)
        x_2 = float(x_2)
        return intercept + (coefficient_1 * x_1) + (coefficient_2 * x_2)

    return twoParameterGlm


def testNumeric(var):
    """
    test if var is numeric
    :params var: any variable
    :returns: true if var is int or float, false otherwise
    """
    if isinstance(var, int):
        return True
    elif isinstance(var, float):
        return True
    else:
        return False


# def writeSamtoolsIndexSbatchScript(self): #TODO: update for nextflow pipeline (different bam file suffix, should already have .bai)
#     """
#        write sbatch script to index bam files
#        This needs to be done b/c indexing can only take place via srun/sbatch (on compute node. samtools not available on login node as of 3/2020)
#        script will be place in rnaseq_pipeline/job_scripts
#    """
#     if not hasattr(self, 'samtools_index_sbatch_script'):
#         self.samtools_index_sbatch_script = os.path.join(self.job_scripts, utils.dirName(
#             self.scratch_alignment_source) + '_igv_index.sbatch')
#
#     job = '#!/bin/bash\n' \
#           '#SBATCH -N 1\n' \
#           '#SBATCH --mem=10G\n' \
#           '#SBATCH -o {0}/index_bams_%A.out\n' \
#           '#SBATCH -e {0}/index_bams_%A.err\n' \
#           '#SBATCH -J index_bams\n'.format(self.sbatch_log)
#     if hasattr(self, 'email'):
#         job += '#SBATCH --mail-type=END,FAIL\n' \
#                '#SBATCH --mail-user=%s\n' % self.email
#     job += '\nml samtools\n'
#
#     for alignment_file in self.bam_file_to_index_list:
#         job += '\nsamtools index -b %s\n' % alignment_file
#
#     with open(self.samtools_index_sbatch_script, 'w') as file:
#     file.write('%s' % job)

def convertFastqFilename(fastq_filename, filetype):
    """
        given a fastq filename, return a bam filename with a given suffix
        :param fastq_filename: a fastqFilename or path with a fastqfilename
        :param filetype: one of ['bam', 'count', 'novoalign_log']
        :returns: the file with .fastq.gz removed and the filetype suffix added name
                  (eg sample1_sorted_alignment_with_annote.bam for filetype 'bam')
    """
    suffix_dict = {'bam': '_sorted_aligned_reads_with_annote.bam',
                   'count': '_read_count.tsv',
                   'novoalign_log': '_novoalign.log'}
    fastq_basename = pathBaseName(fastq_filename)
    return fastq_basename + suffix_dict[filetype]


def extractInfoFromQuerySheet(query_df, fastq_filename, extract_column):
    """
        extract information from query sheet given sample_name from qual_assess_df (which is the basename, no ext, of the fastq.gz)
        :param query_df: a query dataframe from which info will be extracted
        :param fastq_filename: fastq filename used to identify a sample
        :param extract_column: column from which to extract a value from the query_df based on sample_name
        :returns: value extracted from query_df based on sample name and extract column
    """
    # remove rows containing NaN in fastqFileName
    tmp_query_df = query_df[~query_df.fastqFileName.isnull()]
    try:
        extract_value = list(tmp_query_df[tmp_query_df['fastqFileName'].str.contains(fastq_filename)][extract_column])[
            0]
    except AttributeError:
        print('You must pass a query df')

    return str(extract_value)


def extractInfoFromQuerySheet(query_df, sample_name, extract_column):
    """ TODO: THIS WAS REMOVED FROM ORGANISMDATA OR QUALITY ASSESSMENT...DON'T REMEMBER...BUT IT NEEDS TO BE UPDATED IN THE REST OF THE CODEBASE TO USE THIS IN UTILS NOW (AND HANDLE ERRORS IN CALLING SCRIPT)
        extract information from query sheet given sample_name from qual_assess_df (which is the basename, no ext, of the fastq.gz)
        :param sample_name: name of sample -- basename of fastq.gz, no containing directory, no extension
        :param extract_column: column from which to extract a value from the query_df based on sample_name
        :returns: value extracted from query_df based on sample name and extract column
    """
    try:
        extract_value = \
            list(query_df[query_df['fastqFileName'].str.contains(sample_name + '.fastq.gz')][extract_column])[
                0]
    except AttributeError:
        raise AttributeError  # TODO: CLEAN THIS -- I JUST REMOVED THE LOGGER STATEMENT WHILE DEBUGGING IGVSHOTS
    else:
        return str(extract_value)


def extractGenotypeList(query_df_row, genotype_columns=["genotype1", "genotype2"], convert_CNAG_to_CKF44=False):
    """
        for use in a loop over rows of a metadata_df
    """
    genotype_list = query_df_row[genotype_columns].values
    if np.isnan(genotype_list[1]):
        genotype_list[1] = None

    if convert_CNAG_to_CKF44:
        for i in range(len(genotype_list)):
            try:
                genotype_list[i] = genotype_list[i].replace("CNAG", "CKF44")
            except AttributeError:  # this is to catch the case in which genotype[1] is none. TODO: should have more error checking in this function
                pass

    return list(genotype_list)


def extractRunNumber(query_df_row):
    """
        for use in a loop over rows of a metadata_df
    """
    # TODO SEE STANDARD DATA -- ALREADY IMPLEMENTED
    raise NotImplementedError


def postData(url, data, logger, **kwargs):
    """
        post data to url
        :params url: url to the site (this should be complete, eg to update counts http://13.59.167.2/api/Counts) NOTE: pass this url, no id, etc, and use the kwargs option for a put
        :params data: the body of the request
        :throws: HTTPError
    """
    r = requests.post(url, data=data)
    r.raise_for_status()
