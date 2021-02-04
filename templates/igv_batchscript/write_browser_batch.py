#!/usr/bin/env python

# this script will create a batch script for IGV that you can use on your local computer to take a series of browser shots

# first argument is path to metadata .xlsx
# second argument is path to align dir
# third argument the suffix of the bam path (in this case, _sorted_aligned_reads.bam) appended to filename in metadata
# fourth argument is locus in IGV format, eg chr15:start_pos-stop_pos (its a good idea to add some margin on either side. i used 500)
# fifth argument is path to .batch output
# sixth argument is output dir

# ./write_browser_batch.py ../../data/metadata.xlsx $HOME/projects/laura/dlk/align _sorted_aligned_reads.bam chr15:102,498,000-102,518,000 make_browser_shots.batch $PWD

import pandas as pd
import sys
import os

metadata_path = sys.argv[1]

align_dir_path = sys.argv[2]

bam_suffix = sys.argv[3]

locus = sys.argv[4]

batchshot_path = sys.argv[5]

output_dir = sys.argv[6]

metadata_df = pd.read_excel(metadata_path)

batchshot_dict_keys = ["new", "snapshotDirectory", "genome", "maxPanelHeight", "preference", "load", "goto", "sort", "collapse", "snapshot"]

with open(batchshot_path, 'w') as file:
    for index, row in metadata_df.iterrows():
        batchshot_dict = dict(zip(batchshot_dict_keys, [None]*len(batchshot_dict_keys)))
        filename = row.filename
        bam_path = os.path.join(align_dir_path, filename+bam_suffix)
        batchshot_dict["genome"] = "Mouse"
        batchshot_dict["snapshotDirectory"] = output_dir
        batchshot_dict["load"] = bam_path
        batchshot_dict["goto"] = locus
        batchshot_dict["snapshot"] = "_".join([row.genotype, row.perturbation, str(row.time), str(index), row.filename])
        batchshot_dict["maxPanelHeight"] = "300"
        for key, value in batchshot_dict.items():
            append_value = "" if value is None else value
            file.write(key+" "+append_value)
            file.write("\n")
        file.write("\n\n")
