#!/usr/bin/env python

# name: verify_metadata_accuracy
# purpose: check metadata tables for expected entries
# input: an individual sheet, a directory of sheets, or the top most directory
# output: Warnings and prompts to the user to either ignore or change entries
# written by: chase mateusiak(chase.mateusiak@gmail.com)
# date included in rnaseq_pipe: 1/21/2020
# See the end of the script for a description of the environment used to write/test

#import queryDB
import pandas as pd
import sys
import argparse
import re

def main(argv):
    args = parseArgs(argv)

    if not len(args.uniqueKeys) == 0:
        print(args.key, " ", args.sheet_path)
        uniqueKeys(args.key, args.sheet_path)
    print('\n::verify_metadata_complete::')

def checkCSV(file):
    # test whether a given file is a .csv or .xlsx
    if re.search('\.csv', file):
        return True
    else:
        return False

def getKeys(datadir_keys, concat_dict):
    # create list of shared columns between successive pairs of keys in datadir_keys i.e. the columns which are shared between the sheets in fastqFiles and Library
    # Args: a list of subdirectories in datadir and a dictionary of concatenated sheets from within each one of those directories
    # PLEASE NOTE: DICTIONARIES ARE ORDERED IN PYTHON 3.6+. KEY/DIRECTORIES NEED TO BE ADDED TO CONCAT DICT IN ORDER YOU WISH TO MERGE
    # Returns: a list of key columns
    key_cols = []
    num_keys = len(datadir_keys)
    for i in range(num_keys-1):
        key_cols.append(concat_dict[datadir_keys[i]].keys().intersection(concat_dict[datadir_keys[i+1]].keys()))
    return key_cols

def parseArgs(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument('-u', '--uniqueKeys', required=False,
                        help='Enter True to access uniqueKeys function. Must also pass -k and -s')
    parser.add_argument('-k', '--key', nargs='+', default=[],
                        help = 'list of keys to check for uniqueness')
    parser.add_argument('-s', '--sheet_path', type=str, required=False,
                        help='path to data, either .csv or .xlsx')

    return parser.parse_args(argv[1:])

def uniqueKeys(key, df_path):
    # test whether the key columns in a given sheet are unique
    # Args: the keys (passed as list) and the path to the dataframe
    # Return: none
    # to std_out: info on uniqueness of key
    if checkCSV(df_path):
        sheet = pd.read_csv(df_path)
    else:
        sheet = pd.read_excel(df_path)
    # make a tuple of the key columns and store them as pandas series
    key_tuples = sheet[key].apply(tuple, axis=1)
    num_keys = key_tuples.size
    num_unique_keys = key_tuples.unique().size
    print("\nThe number of unique keys is {}. The number of rows is {}. If these are equal, the keys are unique.".format(num_keys, num_unique_keys))
    if not num_keys == num_unique_keys:
        print("\nThe following indicies are not unique:\n\t{}".format(sheet[sheet[key].duplicated()]))
        print("\nDo you want to continue? Enter y or n: ")
        user_response = input()
        if user_response == 'n':
            quit()

def coerceAllCols(sheet):
    # converts columns specific in functions below to floats and datetime
    sheet = floatColCoerce(sheet)

    sheet = datetimeColCoerce(sheet)

    # forcing all to lowercase is not working in the search
    #sheet = strColCoerce(sheet)

    return sheet

def colCoerce(cols, sheet, dtype):
    # general function to coerce column to certain data type
    # Args: a list of columns to coerce (can be a subset of columns of a data frame. must be a list), a dataframe,
    #       and a datatype (see pandas .astype documentation for which datatypes are possible. didn't work for datetime, for example)
    # Return: the dataframe with the specified columns coerced
    coerce_cols = cols

    for col in sheet.columns[sheet.columns.isin(coerce_cols)]:
        sheet[col] = sheet[col].astype(dtype)
    return sheet

def floatColCoerce(sheet):
    # coerce columns specified below to int
    int_columns = ['librarySampleNumber', 'runNumber', 'tapestationConc', 'readsObtained', 'rnaSampleNumber',
                   'inductionDelay', 'replicate', 'timePoint']

    sheet = colCoerce(int_columns, sheet, 'float')

    return sheet

def datetimeColCoerce(sheet):
    # coerce columns specified below to datetime
    datetime_columns = ['libraryDate', 's2cDNADate', 's1cDNADate', 'rnaDate', 'harvestDate']

    for col in datetime_columns:
        sheet[col] = pd.to_datetime(sheet[col])

    return sheet

def strColCoerce(sheet):
    # coerce any column NOT in int_col or datetime_col to a string and make it lower case

    int_columns = ['librarySampleNumber', 'tapestationConc', 'readsObtained', 'rnaSampleNumber','inductionDelay', 'replicate', 'timePoint']
    datetime_columns = ['libraryDate', 's2cDNADate', 's1cDNADate', 'rnaDate', 'harvestDate']

    typed_cols = int_columns + datetime_columns

    for col in sheet.columns:
        if col not in typed_cols:
            sheet[col] = sheet[col].astype(str).apply(lambda x: x.lower())

    return sheet

if __name__ == '__main__':
	main(sys.argv)