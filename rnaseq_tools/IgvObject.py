"""
tools to create browser shots

makeIgvSnapshotDict create a dict with the following format
igv_snapshot_dict = {'wildtype.fastq.gz': {'gene': [gene_1, gene_2, gene_3...],
                     'bam': file.bam, 'bed': file.bed},
                     'perturbed_sample.fastq.gz': {'gene': [gene_1, gene_2],
                     'bam': file.bam, 'bed': file.bed},
                     'perturbed_sample.fastq.gz': {'gene': [gene_1],
                     'bam': file.bam, 'bed': file.bed}
                  }
and ensures that all alignment files (.bam) and their index companions (.bam.bai) are in the experiment directory, which
is either supplied or created in rnaseq_tmp

"""

from rnaseq_tools import utils
from rnaseq_tools.OrganismData import OrganismData
from rnaseq_tools import annotation_tools
import sys
import os
import time
import subprocess


# igv_snapshot_dict = {'wildtype.fastq.gz': {'gene': [gene_1, gene_2, gene_3...], 'bam':
#                  file.bam, 'bed': file_1.bed},
#                  'perturbed_sample.fastq.gz': {'gene': [gene_1, gene_2], 'bam':
#                   file.bam, 'bed': file_1.bed},
#                  'perturbed_sample.fastq.gz': {'gene': [gene_1], 'bam': file.bam,
#                  'bed': file_1.bed}
#                  }

class IgvObject(OrganismData):
    def __init__(self, **kwargs):
        # additional attributes to add to the _attributes in StandardData
        self._igv_attributes = ['sample_list', 'igv_genome', 'output_dir', 'wildtype', 'experiment_dir']
        # initialize Standard data with the extended _attributes
        super(IgvObject, self).__init__(self._igv_attributes, **kwargs)


    def makeIgvSnapshotDict(self):
        """
        from list of samples, create dictionary in format {'sample': {'gene': [gene_1, gene_2,...], 'bed': /path/to/bed, 'bam': /path/to/bam}
        # ASSUMPTION: THE SAMPLE LIST IS EITHER A LIST OF .FASTQ.GZ OR *_read_count.tsv
        """
        # module load samtools
        os.system('ml samtools')
        # raise attribute error if sample list is not found
        if not hasattr(self, 'sample_list'):
            raise AttributeError('Your instance of IgvObject does not have a list of samlpes. '
                                 'It must have a list of samples (list of fastq or count file names) to create the igv_dict.')
        # raise attribute error if no (standardized) query_df. this should be created by passing query_sheet_path to the constructor. See StandardData
        if not hasattr(self, 'query_df'):
            raise AttributeError('No query_df (standardized query_df object in StandardData object). This is necessary '
                                 'and created if you input a query_sheet_path to the constructor of any StandardData '
                                 'object instance. However, you can input the full database (queryDB.py with flag -pf). '
                                 'It will just make the searches for a given sample somewhat longer (though likely not much).')
        # if an experiment_dir is not passed in constructor of IgvObject, create one in rnaseq_tmp/<timenow>_igv_files and store the path in self.experiment_dir
        if not hasattr(self, 'experiment_dir'):
            print('creating a directory in rnaseq_tmp to store alignment files so that the scheduler has access to them.'
                  'This assumes that the count and alignment files have already been processed and moved to '
                  '/lts/mblab/Crypto/rnaseq_data/align_expr')
            timestr = time.strftime("%Y%m%d_%H%M%S")
            setattr(self, 'experiment_dir', os.path.join(self.rnaseq_tmp, '{}_igv_files'.format(timestr)))
            utils.mkdirp(self.experiment_dir)
        igv_snapshot_dict = {}
        genotype_list = []
        wildtype_sample_list = []
        for sample in self.sample_list:
            # strip .fastq.gz if it is there and add _read_count.tsv -- remember that self.query_df which has fastqFileName --> COUNTFILENAME and .fastq.gz converted to _read_count.tsv extensions. All column headings CAPITAL
            if sample.endswith('.fastq.gz'):
                sample = utils.pathBaseName(sample) + '_read_count.tsv'
            genotype = self.extractValueFromStandardRow('COUNTFILENAME', sample, 'GENOTYPE')
            # extract run_number just in case needed to find bam file in align_expr
            run_number = self.extractValueFromStandardRow('COUNTFILENAME', sample, 'RUNNUMBER', check_leading_zero=True)
            # create bamfile name
            bamfile = sample.replace('_read_count.tsv', '_sorted_aligned_reads.bam')
            # if it is not in the exp dir, then add it
            if not os.path.exists(os.path.join(self.experiment_dir, bamfile)):
                prefix = utils.addForwardSlash(self.lts_align_expr)
                bamfile_align_expr_path = '{}run_{}/{}'.format(prefix, run_number, bamfile)
                cmd = 'rsync -aHv {} {}'.format(bamfile_align_expr_path, utils.addForwardSlash(self.experiment_dir))
                utils.executeSubProcess(cmd)
            # store full path to bam file in experiment dir
            bamfile_fullpath = os.path.join(self.experiment_dir, bamfile)
            if not os.path.exists(bamfile_fullpath + '.bai'): # test if indexed bam exists
                self.indexBam(bamfile_fullpath)
            # split on period if this is a double perturbation. Regardless of whether a '.' is present,
            # genotype will be cast to a list eg ['CNAG_00000'] or ['CNAG_05420', 'CNAG_01438']
            genotype = genotype.split('.')
            # if genotype ends with _over, remove _over
            for index in range(len(genotype)):
                genotype[index] = genotype[index].replace('_over', '')
            # if the object has an attribute wildtype, and genotype is not wildtype, add to igv_snapshot_dict
            if hasattr(self, 'wildtype') and genotype[0] == self.wildtype:
                # add the genotypes to genotype_list (not as a list of list, but as a list of genotypes)
                genotype_list.extend(genotype)
                # add to igv_snapshot_dict
                igv_snapshot_dict.setdefault(sample, {}).setdefault('gene', []).extend(genotype) # TODO: clean this up into a single line, make function to avoid repeated code below w/wildtype
                igv_snapshot_dict[sample]['bam'] = bamfile_fullpath
                igv_snapshot_dict[sample]['bed'] = None
            # if genotype is equal to wildtype, then store the sample as the wildtype (only one, check if this is right)
            else:
                igv_snapshot_dict.setdefault(sample, {'gene': None, 'bam': None, 'bed': None})
                wildtype_sample_list.append([sample, bamfile_fullpath]) # wildtype_sample_list will be a list of lists
        # if the wildtype genotype was found, create entry in the following form
        # {sample_read_counts.tsv: {'gene': [perturbed_gene_1, perturbed_gene_2, ...], 'bam': wt.bam, 'bed': created_bed.bed}
        for wt_sample in wildtype_sample_list:
            igv_snapshot_dict[wt_sample[0]]['gene'] = genotype_list
            igv_snapshot_dict[wt_sample[0]]['bam'] = wt_sample[1]
            igv_snapshot_dict[wt_sample[0]]['bed'] = None
        # set attribute pointing toward igv_snapshot_dict
        setattr(self, 'igv_snapshot_dict', igv_snapshot_dict)

    def indexBam(self, bamfile_fullpath):
        """
        Index the bam files. The .bai file (indexed bam) will be deposited in the same directory as the bam file itself.
        for igv, as long as the .bai is in the same directory as the .bam file, it will work.
        :param bamfile_fullpath: full path (absolute) to bamfile
        """
        print('\nindexing {}'.format(bamfile_fullpath))
        # subprocess.call will wait until subprocess is complete
        exit_status = subprocess.call('samtools index -b {}'.format(bamfile_fullpath), shell = True)
        if exit_status == 0:
            print('\nindexing complete. .bai is deposited in {}'.format(self.experiment_dir))
        else:
            sys.exit('failed to index {}. Cannot continue'.format(bamfile_fullpath))

    def createBedFile(self, flanking_region = 500, file_format='png'):
        """
        Create bed files to describe IGV region of interest
        :param flanking_region: how far up/down stream from the gene to capture in the snapshot
        :param file_format: what format to use for image file
        """
        ## get gene dictionary with chromsome, gene coordinates, strand
        if self.annotation_file.endswith('gtf'):
            self.annotation_dict = annotation_tools.parseGtf(self.annotation_file)
        elif self.annotation_file.endswith('gff') or self.annotation_file.endswith('gff3'):
            self.annotation_dict = annotation_tools.parseGff3(self.annotation_file)
        else:
            sys.exit("ERROR: The gene annotation format cannot be recognized.") # TODO: clean up preceeding blocks -- move parseGFF to OrganismData
        ## create gene body region bed file
        for sample in self.igv_snapshot_dict.keys():
            igv_bed_filepath = os.path.join(self.experiment_dir, utils.pathBaseName(sample) + '.bed')
            self.igv_snapshot_dict[sample]['bed'] = igv_bed_filepath
            writer = open(igv_bed_filepath, 'w')
            for gene in self.igv_snapshot_dict[sample]['gene']:
                d = self.annotation_dict[gene]
                writer.write('%s\t%d\t%d\t[%s]%s.%s\n' % \
                             (d['chrm'], d['coords'][0] - flanking_region,
                              d['coords'][1] + flanking_region, sample, gene, file_format))


    def writeIgvJobScript(ineffmut_dict, igv_genome, igv_output_dir, fig_format='png', email=None,
                          job_script='job_scripts/igv_snapshot.sbatch'):
        """
        Write sbatch job script to make IGV snapshot
        """
        num_samples = len(ineffmut_dict.keys())
        job = '#!/bin/bash\n#SBATCH -N 1\n#SBATCH --mem=5G\n'
        job += '#SBATCH -D ./\n#SBATCH -o log/igv_snapshot_%A.out\n#SBATCH ' \
               '-e log/igv_snapshot_%A.err\n#SBATCH -J igv_snapshot\n'
        if email is not None:
            job += '#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=%s\n' % email
        job += '\nml java\n'
        for sample in ineffmut_dict.keys():
            bam_file = ineffmut_dict[sample]['bam']
            bed_file = ineffmut_dict[sample]['bed']
            # this is a call to another script in the rnaseq_pipeline/tools
            job += 'python -u tools/make_IGV_snapshots.py %s ' \
                   '-bin /opt/apps/igv/2.4.7/igv.jar -nf4 -r %s -g %s -fig_format %s -o %s\n' \
                   % (bam_file, bed_file, igv_genome, fig_format, igv_output_dir)
        # write job to script
        writer = open(job_script, 'w')
        writer.write('%s' % job)
        writer.close()

