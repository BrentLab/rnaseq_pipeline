#!/usr/bin/env Rscript

# Title     : log2_counts_per_million
# Objective : convert raw counts to log2 counts per million using edgeR cpm().
#             From the edgeR documentation: If log-values are computed, then a small count, given by prior.count
#             but scaled to be proportional to the library size, is added to y to avoid taking the log of zero.
# Created by: chasem@wustl.edu chase.mateusiak@gmail.com
# Created on: 3/17/20
# Updated: 07/18/2020 to remove nctr RNA from KN99 sample counts

suppressMessages(library(optparse))
suppressMessages(library(edgeR))
suppressMessages(library(tidyverse))


# TODO: ERROR WHEN THERE IS ONLY ONE SMAPLE IN RAW COUNTS. COLUMN NAMES NOT WRITTEN IN THIS CASE
main = function(args){

  # parse cmd line arguments
  parsed = args
  # rename inputs
  print('...parsing cmd line input')
  path_to_raw_counts = parsed$raw_counts
  organism = parsed$organism

  sprintf('...Reading in raw counts %s', path_to_raw_counts)
  # read in raw counts data.frame
  raw_counts = read.csv(path_to_raw_counts, row.names = 'gene_id', check.names = FALSE) # without check.names = FALSE, R will insert X in front of colnames that start with a number https://stackoverflow.com/a/58951644/9708266

  # for KN99, remove nonconding, transfer, and ribosomal RNA (in the genome, these are currently annotated by taking
  # the top blast hits with some more filtering (see genome_files supplemental file) of the nc,t,r - RNA annotations in H99 against KN99
  if (organism == 'KN99'){
      print('...filtering out nctrRNA genes from KN99 counts')
      # create boolean vector, with TRUE for index rows containing a nctr RNA annotation (CNAG_12345 for example)
      # NOTE: currently, the drug markers are identified as CNAG_NAT and CNAG_G418
      # NOTE: the ! negates the filter, so using this diretly will only return CKF44 and the drug markers CNAG_NAT and CNAG_G418
      nctr_rna_filter = !grepl('CNAG_[[:digit:]]+', rownames(raw_counts))
      raw_counts = raw_counts[nctr_rna_filter, ]
  }

  print('...Creating log2_cpm from the raw counts')
  # convert to edgeR DGEList object
  dgelist = DGEList(raw_counts)
  # cpm returns the log2 of counts per million
  log2_cpm = cpm(dgelist, log=TRUE)
  # convert index row to column 1
  log2_cpm = rownames_to_column(as.data.frame(log2_cpm))
  # rename the converted index row to gene_id
  colnames(log2_cpm)[1] = 'gene_id'

  # write to output path -- NOTE: in cmd line input, the -o is the FULL output path (including filename and extension)
  output_full_path = str_replace(path_to_raw_counts, 'raw_count.csv', 'log2_cpm.csv')
  sprintf("Writing log2_cpm matrix to: %s", output_full_path)
  write_csv(log2_cpm, output_full_path)

} # end main()

parseArguments = function() {
  option_list = list(
    make_option(c('-r', '--raw_counts'),
                help='[REQUIRED]raw count matrix produced by raw_counts.py'),
    make_option(c('-g', '--organism'),
                 help='[REQUIRED]Either KN99 or None (if your sample is KN99, you must enter KN99)'))
  args = parse_args(OptionParser(option_list=option_list))
  return(args)
} # end parseArguments()

main(parseArguments()) # call main method
