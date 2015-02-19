#!/usr/bin/env python

__author__ = 'olga'

import argparse
from glob import glob
import os
import sys

import pandas as pd

class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = parser = argparse.ArgumentParser(
            description='Combine sailfish output files')
        parser.add_argument('-g', '--glob-command', required=False,
                            default='*sailfish', type=str, action='store',
                            help='Where to find sailfish output directories. '
                                 'Default is folders in the current directory '
                                 'whose names end with "sailfish"')
        parser.add_argument('-o', '--out-dir', required=False,
                            default='combined_output', type=str,
                            action='store',
                            help='Where to output the combined matrices. Does '
                                 'not need to exist already. '
                                 'Default is to create a folder called '
                                 '"combined_input"')
        parser.add_argument('-n', '--n-progress', required=False, type=int,
                            default=10, action='store',
                            help="Number of files to show per iterative "
                                 "progress, e.g. 10/58 files completed, "
                                 "20/58 files completed. Can increase this if"
                                 "you have thousands of files, or decrease to"
                                 "1 if you only have a few")
        if inOpts is None:
            self.args = vars(self.parser.parse_args())
        else:
            self.args = vars(self.parser.parse_args(inOpts))

    def do_usage_and_die(self, str):
        '''
        If a critical error is encountered, where it is suspected that the
        program is not being called with consistent parameters or data, this
        method will write out an error string (str), then terminate execution
        of the program.
        '''
        import sys

        print >> sys.stderr, str
        self.parser.print_usage()
        return 2

# Class: Usage
class Usage(Exception):
    '''
    Used to signal a Usage error, evoking a usage statement and eventual
    exit when raised
    '''

    def __init__(self, msg):
        self.msg = msg


class CombineSailfish(object):
    def __init__(self, glob_command, out_dir, n_progress):
        """Any CamelCase here is directly copied from the STAR inputs for
        complete compatibility
        """
        # Make the directory if it's not there already
        if n_progress < 1:
            raise ValueError('"n_progress" must be 1 or greater')

        out_dir = os.path.abspath(os.path.expanduser(out_dir))
        try:
            os.mkdir(out_dir)
        except OSError:
            pass

        tpm_dfs = []
        columns = ['transcript', 'length', 'tpm', 'rpkm', 'kpkm',
                   'EstimatedNumKmers', 'EstimatedNumReads']

        glob_command = '{}/quant_bias_corrected.sf'.format(glob_command)
        filenames = glob(glob_command)
        n_files = len(filenames)

        sys.stdout.write("Reading {} sailfish's quant_bias_corrected.sf "
                         "files ...\n".format(n_files))

        for i, filename in enumerate(filenames):
            # Read "tabluar" data, separated by tabs.
            # Arguments:
            # skiprows=5      Skip the first 5 rows
            # names=columns   Use the column names in the list "columns"
            # index_col=0     The first column is the row names (the row names
            #                 are called the "index" in pandas terms)
            df = pd.read_table(filename, skiprows=5, names=columns,
                               index_col=0)

            # Get the "series" (aka single column) of TPM
            tpm = df.tpm

            # To get the sample ID, split by the folder identifier, "/",
            # and take the second-to-last item (via "[-2]"), which has the
            # sample id, then split on the period, and take the first item
            # via "[0]"
            sample_id = filename.split('/')[-2].split('.')[0]

            # Change the name of the
            tpm.name = sample_id
            tpm_dfs.append(tpm)

            if (i+1) % n_progress == 0:
                sys.stdout.write("\t{}/{} files read\n".format(i+1, n_files))
        sys.stdout.write("\tDone.\n")
        tpm = pd.concat(tpm_dfs, axis=1)

        sys.stdout.write("Separating out spike-ins from regular genes ...\n")
        # Get nonstandard genes, i.e. everything that's not an ensembl ID
        spikein_columns = tpm.columns.map(lambda x: not x.startswith('ENST'))
        tpm_spikein = tpm.ix[:, spikein_columns]
        sys.stdout.write("\tDone.\n")

        sys.stdout.write("Summing TPM expression of all transcripts in a "
                         "gene ...\n")
        # Sum expression of all transcripts of a gene
        tpm_transcripts = tpm.ix[:, ~spikein_columns]
        ensembl_ids = tpm_transcripts.columns.map(
            lambda x: x.split('|')[1].split('.')[0])
        tpm_transcripts.columns = ensembl_ids
        tpm_genes = tpm_transcripts.groupby(level=0, axis=1).sum()
        sys.stdout.write("\tDone.\n")
        import pdb;pdb.set_trace()

        # Save the output files
        filename_to_df = {'tpm_spikein.csv': tpm_spikein,
                          'tpm_genes.csv': tpm_genes}

        sys.stdout.write("Writing output files ...\n")
        for filename, df in filename_to_df.items():
            full_filename = '{}/{}'.format(out_dir, filename)
            df.to_csv(full_filename)
            sys.stdout.write("\tWrote {}\n".format(full_filename))
        sys.stdout.write("Done, son.\n")


if __name__ == '__main__':
    try:
        cl = CommandLine()



        CombineSailfish(cl.args['glob_command'], cl.args['out_dir'],
                        cl.args['n_progress'])
    except Usage, err:
        cl.do_usage_and_die()