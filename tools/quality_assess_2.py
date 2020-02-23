#!/usr/bin/env python
import sys
import argparse
import re
import pandas as pd
import numpy as np
from utils import *
import os

def main(argv):
	parsed = parse_args(argv)
	## validate args
	output_dir = parsed.output
	experiment_dir = os.path.dirname(parsed.experiment_directory)
	filename = experiment_dir + '_quality_summary.xlsx'
	output_name = os.path.join(output_dir, filename)

	if os.path.exists(output_name):
		sys.exit('WARNING: %s already exists, rename the file to proceed.' % output_name)
	if not os.path.exists(parsed.count_matrix):
		sys.exit('ERROR: %s does not exist.' % parsed.count_matrix)
	if not os.path.exists(os.path.dirname(output_name)):
		sys.exit('ERROR: %s does not exist.' % os.path.dirname(output_name))

	## load QC config data
	## TODO: complexity.thresh <- mean(alignment.sum$COMPLEXITY[indx]) - 2*sd(alignment.sum$COMPLEXITY[indx]);
	global QC_dict
	QC_dict = load_config(parsed.qc_configure)
	## get conditions
	conditions = None if parsed.condition_descriptors is None else \
				[c.strip() for c in parsed.condition_descriptors.split(',')]

	## do QA
	print('... Preparing QA dataframe')
	if parsed.resistance_cassettes is None:
		resistance_cassettes = None
		resistance_cassettes_columns = []
	else:
		resistance_cassettes = [rc.strip() for rc in parsed.resistance_cassettes.split(',')]
		resistance_cassettes_columns = [rc+'_FOM' for rc in resistance_cassettes]
	df_columns = ['GENOTYPE','REPLICATE','FASTQFILENAME'] \
				+ conditions \
				+ ['STATUS', 'AUTO_AUDIT', 'MANUAL_AUDIT', 'USER', 'NOTE'] \
				+ ['TOTAL','ALIGN_PCT','MUT_FOW'] \
				+ resistance_cassettes_columns \
				+ ['COV_MED_REP'+''.join(np.array(combo, dtype=str)) for combo in make_combinations(range(1,parsed.max_replicates+1))]
	df, rep_max = initialize_dataframe(parsed.query_sheet, df_columns, conditions)
	if rep_max != parsed.max_replicates:
		print('The max number of replicates is {}. Please re-launch this script with -r {}.'.format(rep_max, rep_max))
	if rep_max > 7:
		print('Calculating the power set of 7 replicates for the CoV will result in many columns. Are you sure you want to proceed? Enter y or n')
		user_response = input()
		if user_response == 'n':
			quit()
	expr, sample_dict = load_expression_data(df, parsed.count_matrix, parsed.gene_list, conditions)
	print('... Assessing reads mapping')
	df = assess_mapping_quality(df, parsed.experiment_directory)
	print('... Assessing efficiency of gene mutation')
	if parsed.descriptors_specific_fow:
		df = assess_efficient_mutation(df, expr, sample_dict, parsed.wildtype, conditions)
	else:
		df = assess_efficient_mutation(df, expr, sample_dict, parsed.wildtype)
	print('... Assessing insertion of resistance cassette')
	df = assess_resistance_cassettes(df, expr, resistance_cassettes, parsed.wildtype)
	print('... Assessing concordance among replicates')
	df = assess_replicate_concordance(df, expr, sample_dict, conditions)
	print('... Auto auditing')
	df = update_auto_audit(df, parsed.auto_audit_threshold)
	save_dataframe(output_name, df, df_columns, conditions, len(conditions))

def parse_args(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('-qs', '--query_sheet', required=True,
						help='[REQUIRED] The output of queryDB that was used to create this experiment')
	parser.add_argument('-e', '--experiment_directory', required=True,
						help='[REQUIRED] the path to the experiment directory created by create_experiment')
	parser.add_argument('-r', '--max_replicates', required=True, type=int,
						help='[REQUIRED] Maximal number of replicate in experiment design.')
	parser.add_argument('-o', '--output', required=True,
						help='[REQUIRED] directory in which to deposit the sample quality summary.')
	parser.add_argument('-c', '--count_matrix', required=True,
						help='[REQUIRED] Normalized count matrix. If not given, the filepath will be guessed based on analysis group number.')
	parser.add_argument('-l', '--gene_list',
						help='Use a custom gene list other than the list in gene annotation file.')
	parser.add_argument('-w', '--wildtype',
						help='Wildtype genotype, e.g. CNAG_00000 for crypto, BY4741 for yeast.')
	parser.add_argument('-m', '--resistance_cassettes',
						help='Resistance cassettes inserted to replace the deleted genes. Use "," as delimiter if multiple cassettes exist.')
	parser.add_argument('--condition_descriptors', default='TREATMENT,TIMEPOINT',
						help='Experimental conditions that describe the sample are used to identify subgroups within each genotype. Use delimiter "," if multiple descriptors are used.')
	parser.add_argument('--descriptors_specific_fow', action='store_true',
						help = 'Set this flag to find the wildtype samples that match the condition descriptors of the mutant sample when calcualting the fold change over wildtype (FOW).')
	parser.add_argument('--qc_configure', default='/opt/apps/labs/mblab/software/rnaseq_pipeline/1.0/templates/qc_config.yaml',
						help='Configuration file for quality assessment.')
	parser.add_argument('--auto_audit_threshold', type=int, default=0,
						help='Threshold for automatical sample audit.')
	return parser.parse_args(argv[1:])


def initialize_dataframe(query, df_cols, conditions):
	"""
	Define the QC dataframe.
	"""
	df1 = pd.DataFrame(columns = df_cols)
	# read in queryDB dataframe describing experiment
	if checkCSV(query):
		df2 = pd.read_csv(query, dtype=np.str)
	else:
		df2 = pd.read_excel(query, dtype=np.str)
	# force all column headers to upper case
	df2.columns = df2.columns.str.upper()

	# cast replicate to int
	df2 = df2.astype({'REPLICATE': 'float'})
	df2 = df2.astype({'REPLICATE': 'int32'})

	if conditions:
		conditions = [x.upper() for x in conditions]

	df2 = df2[['GENOTYPE','REPLICATE','FASTQFILENAME'] + conditions]
	df2 = df2.reset_index().drop(['index'], axis=1)
	df2 = pd.concat([df2, pd.Series([0]*df2.shape[0], name='STATUS')], axis=1)
	df2 = pd.concat([df2, pd.Series([np.nan]*df2.shape[0], name='AUTO_AUDIT')], axis=1)

	# create sample_summary dataframe
	sample_summary_df = df1.append(df2)
	## re-index replicates in case there are technical replicates of a biological replicate
	# add GENOTYPE to conditions
	conditions.append('GENOTYPE')
	sample_summary_df['REPLICATE'] =  sample_summary_df.groupby(conditions).cumcount() +1
	rep_max = sample_summary_df.groupby(conditions).cumcount().max()

	return sample_summary_df, rep_max


def load_expression_data(df, cnt_mtx, gene_list, conditions):
	"""
	Load count matrix, and make a sample dictionary.
	"""
	## load count matrix
	count = pd.read_csv(cnt_mtx)
	count = count.rename(columns={'Unnamed: 0':'gene'})
	## find the intersected gene list
	if gene_list is not None:
		gids = pd.read_csv(gene_list, names=['gene'])
		if len(np.setdiff1d(gids, count['gene'])) > 0:
			print('WARNING: The custom gene list contains genes that are not in count matrix. Proceeding using the intersection.')
		gids = np.intersect1d(gids, count['gene'])
		count = count.loc[count['gene'].isin(gids)]
	## make sample dict with (genotype, condiition1, condition2, ...) as the key
	sample_dict = {}
	for i,row in df.iterrows():
		genotype = row['GENOTYPE']
		key = tuple([genotype]) if len(conditions) == 0 else \
				tuple([genotype] + [row[c] for c in conditions])
		#sample = str(row['FASTQFILENAME'])
		sample = 'JRtimeCourse/' + fileBaseName(row['FASTQFILENAME']) +'_read_count.tsv'
		if sample in count.columns.values:
			if key not in sample_dict.keys():
				sample_dict[key] = {}
			sample_dict[key][row['REPLICATE']] = sample
	return count, sample_dict


def assess_mapping_quality(df, exp_dir, aligner_tool='novoalign'):
	"""
	Assess percentage of uniquely mapped reads over all reads.
	"""
	for i,row in df.iterrows():
		# fileBaseName from utils
		sample = fileBaseName(str(row['FASTQFILENAME'])) + '_' + aligner_tool + '.log'
		filepath = os.path.join(exp_dir, sample)
		## read alignment log
		reader = open(filepath, 'r')
		lines = reader.readlines()
		for line in lines:
			reg_total = re.search(r'Read Sequences:( +)(.\d+)', line)
			reg_uniq = re.search(r'Unique Alignment:( +)(.\d+)', line)
			if reg_total:
				total_reads = int(reg_total.group(2))
			if reg_uniq:
				uniq_mapped_reads = int(reg_uniq.group(2))
		reader.close()
		align_pct = uniq_mapped_reads/float(total_reads)
		## set mapping quality
		row['TOTAL'] = total_reads # read Sequences
		row['ALIGN_PCT'] = align_pct # Unique Alignment
		if total_reads < QC_dict['TOTAL_READS']['threshold']:
			row['STATUS'] += QC_dict['TOTAL_READS']['status']
		if align_pct < QC_dict['ALIGN_PCT']['threshold']:
			row['STATUS'] += QC_dict['ALIGN_PCT']['status']
		df.iloc[i] = row
	return df


def assess_efficient_mutation(df, expr, sample_dict, wt, conditions=None):
	"""
	Assess the completeness of gene deletion or efficiency of gene 
	overexpression by caluclating the expression of the perturbed gene
	in mutant sample over mean expression of the same gene in wildtype.
	"""
	if wt is None:
		return df
	descr_match = True if conditions is not None else False
	## get wildtype samples if not matching descriptors
	if not descr_match:
		wt_samples = [] 
		for key in sample_dict.keys(): 
			if key[0] == wt:
				wt_samples += sample_dict[key].values()
		## calculate mean expression level of each gene
		wt_expr = pd.Series(pd.DataFrame.mean(expr[wt_samples], 
							axis=1), name='mean_fpkm')
		wt_expr = pd.concat([expr, wt_expr], axis=1)
	## calculate efficiency of gene deletion, ignoring overexpression(*_over)
	for i,row in df[df['GENOTYPE'] != wt].iterrows():
		sample = str(row['FASTQFILENAME'])
		## check for each mutant gene (there could be multiple mutant genes, delimited by '.')
		mut_fow_list = []
		for mut_gene in row['GENOTYPE'].split('.'):
			## get wildtype samples if not matching descriptors
			if descr_match:
				mut_descr = [row[c] for c in conditions]
				wt_samples = [] 
				for key in sample_dict.keys(): 
					descr_matched = all([key[j+1] == mut_descr[j] for j in range(len(mut_descr))])
					if key[0] == wt and descr_matched:
						wt_samples += sample_dict[key].values()
				if len(wt_samples) == 0:
					print('\tSample %s has no WT sample that matches its condition descriptors. Skipping this sample' % sample)
					continue
				## calculate mean expression level of each gene
				wt_expr = pd.Series(pd.DataFrame.mean(expr[wt_samples], 
									axis=1), name='mean_fpkm')
				wt_expr = pd.concat([expr, wt_expr], axis=1)
			## get mutant gene expression in mutatnt sample 
			mut_gene2 = mut_gene.strip("_over")
			if mut_gene2 not in expr['gene'].tolist():
				print('\t%s not in gene list. Skipping this genotype' % mut_gene2)
				continue
			wt_mean = float(wt_expr[wt_expr['gene'] == mut_gene2]['mean_fpkm'])
			if wt_mean == 0:
				print('\t%s has 0 mean expression in WT samples' % mut_gene2)
				mut_fow = np.inf
			else:
				mut_fow = float(expr[expr['gene'] == mut_gene2][sample])/wt_mean
			
			if mut_gene.endswith('_over'):
				## check overexpression
				if (mut_fow < QC_dict['MUT_FOW']['OVEREXPRESSION']['threshold']) and (row['STATUS'] < QC_dict['MUT_FOW']['OVEREXPRESSION']['status']):
					row['STATUS'] += QC_dict['MUT_FOW']['OVEREXPRESSION']['status']
			else:
				## check deletion
				if (mut_fow > QC_dict['MUT_FOW']['DELETION']['threshold']) and (row['STATUS'] < QC_dict['MUT_FOW']['DELETION']['status']):
					row['STATUS'] += QC_dict['MUT_FOW']['DELETION']['status']
			mut_fow_list.append(str(mut_fow))
		row['MUT_FOW'] = ','.join(mut_fow_list)
		df.iloc[i] = row
	return df


def assess_resistance_cassettes(df, expr, resi_cass, wt):
	"""
	Assess drug resistance marker gene expression, making sure the proper
	marker gene is swapped in place of the perturbed gene.
	"""
	if resi_cass is None:
		return df
	## get the median of resistance cassettes
	rc_med_dict = {}
	mut_samples = [s for s in expr.columns.values if (not s.startswith(wt)) and (s != 'gene')]
	for rc in resi_cass:
		## exclude wildtypes and markers expressed < 150 normalized counts
		rc_fpkm = expr.loc[expr['gene'] == rc, mut_samples]
		rc_fpkm = rc_fpkm.loc[:, (np.sum(rc_fpkm, axis=0) > 150)]
		rc_med_dict[rc] = rc_fom = np.nan if rc_fpkm.empty else np.median(rc_fpkm)
	## calcualte FOM (fold change over mutant) of the resistance cassette
	for i,row in df.iterrows():
		genotype = row['GENOTYPE']
		sample = genotype +'-'+ str(row['FASTQFILENAME'])
		## update FOM
		for rc in rc_med_dict.keys():
			row[rc+'_FOM'] = np.nan if np.isnan(rc_med_dict[rc]) else float(expr.loc[expr['gene'] == rc, sample])/rc_med_dict[rc]
		## flag those two problems: 
		## 1. the resistance cassette is exprssed in WT
		## 2. more than one resistance cassette is expressed in single mutant
		## TODO: add criteria for multi-mutants 
		fom_check = [row[rc+'_FOM'] > QC_dict['MARKER_FOM']['threshold']*rc_med_dict[rc] for rc in rc_med_dict.keys()]
		if genotype == wt and sum(fom_check) > 0:
			row['STATUS'] += QC_dict['MARKER_FOM']['status']
		if genotype != wt and len(genotype.split('.')) > 1 and sum(fom_check) > 1:
			row['STATUS'] += QC_dict['MARKER_FOM']['status']
	return df


def assess_replicate_concordance(df, expr, sample_dict, conditions):
	"""
	Assess the concordance among the replicates of each genotype by calculating
	the COV of each combination of replicates. Then find the maximal number of
	concordant replicates.
	"""

	## calcualte COV medians for replicate combinations
	for key in sorted(sample_dict.keys()):
		sample_ids = [s for s in sample_dict[key].items()]
		cov_meds_dict = {}
		rep_combos = make_combinations(sample_dict[key].keys())
		for rep_combo in rep_combos:
			## sort as integers
			#rep_combo = np.array(sorted(np.array(rep_combo,dtype=int)), dtype=str)
			rep_num = len(rep_combo)
			sample_combo = [sample_dict[key][rep] for rep in rep_combo]
			## calculate COV median
			cov_median = calculate_cov_median(expr[sample_combo])
			rep_combo_col = 'COV_MED_REP'+''.join(np.array(rep_combo, dtype=str))
			fastq_file_names = [os.path.basename(fileBaseName(tup[1])[:-11]) +'.fastq.gz' for tup in sample_ids]
			df.loc[df['FASTQFILENAME'].isin(fastq_file_names), rep_combo_col] = cov_median
			## store COV median at the respective rep number
			if rep_num not in cov_meds_dict.keys():
				cov_meds_dict[rep_num] = {'rep_combos': [], 'cov_meds': []}
			cov_meds_dict[rep_num]['rep_combos'].append(rep_combo)
			cov_meds_dict[rep_num]['cov_meds'].append(cov_median)
		## find the maximal number of replicates that pass concordance threshold
		for rep_num in sorted(cov_meds_dict.keys())[::-1]:
			rep_combo = cov_meds_dict[rep_num]['rep_combos']
			cov_meds = cov_meds_dict[rep_num]['cov_meds']
			if sum([c < QC_dict['COV_MED']['threshold'] for c in cov_meds]) > 0:
				best_combo = rep_combo[np.argmin(cov_meds)]
				break 
		max_rep_combo = cov_meds_dict[max(cov_meds_dict.keys())]['rep_combos'][0]
		outlier_reps = set(max_rep_combo) - set(best_combo)
		## update status
		for rep in outlier_reps:
			outlier_indx = set(df.index[(df['GENOTYPE'] == key[0]) & \
							(df['REPLICATE'] == rep)])
			for ci in range(len(conditions)):
				outlier_indx = outlier_indx & \
							set(df.index[df[conditions[ci]] == key[ci+1]])
			df.loc[list(outlier_indx), 'STATUS'] += QC_dict['COV_MED']['status']
	return df


def calculate_cov_median(x):
	"""
	Calculate the median of COVs (coefficient of variation) among replicates
	"""
	covs = np.std(x, axis=1) / np.mean(x, axis=1)
	return np.nanmedian(covs)


def update_auto_audit(df, threshold):
	"""
	Automatically flag sample with status over threshold 
	"""
	df.loc[df['STATUS'] > threshold, 'AUTO_AUDIT'] = 1
	return df


def save_dataframe(filepath, df, df_cols, conditions, fp_ext=0):
	"""
	Save dataframe of quality assessment
	"""
	df = df.sort_values(['GENOTYPE'] + conditions + ['REPLICATE'])
	#if not filepath.endswith('.xlsx'):
	#	filepath += '.xlsx'
	df.to_excel(filepath, columns=df_cols, index=False, freeze_panes=(1,3+fp_ext))
	

if __name__ == '__main__':
	main(sys.argv)
