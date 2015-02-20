#!/usr/bin/env python

__author__ = 'olga'

import argparse
from glob import iglob
import os
import re
import sys

import numpy as np
import pandas as pd

class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = parser = argparse.ArgumentParser(
            description='Combine sailfish output files')
        parser.add_argument('-g', '--glob-command', required=False,
                            default='./miso/*/*/summary/*.miso_summary',
                            type=str, action='store',
                            help='Where to find sailfish output directories. '
                                 'The default assumes you have your MISO '
                                 'summary files are located as so: '
                                 '<current_directory>/miso/'
                                 '<sample_id>/<splice_type>/summary/'
                                 '<splice_type>.miso_summary')
        parser.add_argument('-o', '--out-dir', required=False,
                            default='./combined_output', type=str,
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
        parser.add_argument('--ci-max', required=False, type=float,
                            default=0.5, action='store',
                            help="Used for filtering. Maximum size of the "
                                 "confidence interval for the "
                                 "'percent-spliced-in' aka PSI score of a "
                                 "splicing event")
        parser.add_argument('--per-isoform-counts-min', required=False,
                            type=int, default=10, action='store',
                            help="Used for filtering. Minimum number of reads"
                                 "that MISO deemed unique to one isoform, for "
                                 "that splicing event to pass filtering.")
        parser.add_argument('--downsampled', required=False,
                            action='store_true',
                            help='If given, then assumed that the samples '
                                 'were downsampled at certain numbers of '
                                 'reads, and assumes your summary files are '
                                 'located in,'
                                 '<current_directory>/miso/<sample_id>_prob<p>'
                                 '_iter<i>/<splice_type>/summary/')
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


class CombineMiso(object):
    def __init__(self, glob_command, out_dir='./combined_outputs',
                 n_progress=10, downsampled=False, ci_max=0.5,
                 per_isoform_counts_min=10):
        """Combine MISO output files and write to disk

        Parameters
        ----------
        glob_command : str
            Where to find sailfish output directories
        out_dir : str
            Where to output the combined matrices. Will be created if it
            doesn't exist
        n_progress : int
            Integer step size to show progress. E.g. for 10/58 completed

        """
        out_dir = out_dir.rstrip('/')
        out_dir = os.path.abspath(os.path.expanduser(out_dir))

        if n_progress < 1:
            raise ValueError('"n_progress" must be 1 or greater')

        # Make the directory if it's not there already
        out_dir = os.path.abspath(os.path.expanduser(out_dir))
        try:
            os.mkdir(out_dir)
        except OSError:
            pass
        dfs = []

        filenames = iglob(glob_command)
        n_files = sum(1 for i in filenames)
        sys.stdout.write("Reading {} MISO summary files ...\n".format(n_files))

        # re-initialize iterator
        filenames = iglob(glob_command)

        n_files_true = 0
        for i, filename in enumerate(filenames):
            # Check that more than just the header is there
            if os.path.getsize(filename) > 113:
                n_files_true += 1
                df = self.read_miso_summary(filename)

                splice_type = os.path.basename(filename).split('.')[0]
                sample_id = filename.split('/')[-4]

                if downsampled:
                    fragments = sample_id.split('_')
                    real_id = '_'.join(fragments[:-2])
                    probability = float(fragments[-2].lstrip('prob'))
                    iteration = int(fragments[-1].lstrip('iter'))
                    sys.stdout.write('\t{}\t{}\t{}\t{}\n'.format(i, real_id,
                                                                 probability,
                                                                 iteration))
                else:
                    real_id = filename
                    sys.stdout.write('\t{}\t{}\n'.format(i, real_id))

                df['sample_id'] = sample_id
                df['splice_type'] = splice_type

                if downsampled:
                    df['probability'] = probability
                    df['iteration'] = iteration

                dfs.append(df.reset_index())
                if (i + 1) % n_progress == 0:
                    sys.stdout.write(
                        "\t{}/{} files attempted to read\n".format(i + 1, n_files))
            else:
                sys.stdout.write("\tOnly found header and an empty table for "
                                 "{}\n".format(filename))
        sys.stdout.write("\tDone.\n")

        sys.stdout.write("Merging all {} MISO summaries into a gigantic "
                         "one ...\n".format(n_files_true))
        summary = pd.concat(dfs)
        sys.stdout.write("\tDone.\n")

        sys.stdout.write("Writing raw MISO summary files ...\n")
        csv = '{}/miso_summary_raw.csv'.format(out_dir)
        summary.to_csv(csv)
        sys.stdout.write("\tWrote {}\n".format(csv))

        sys.stdout.write("Filtering MISO summaries with ci_max={}, "
                         "per_isoform_counts={} ...\n".format(
            ci_max, per_isoform_counts_min))
        summary = self.filter_miso_summary(summary, ci_max,
                                           per_isoform_counts_min)
        sys.stdout.write("\tDone.\n")

        if downsampled:
            sys.stdout.write("Sorting downsampled files by probability "
                             "and iteration...\n")
            summary = summary.sort(columns=['probability', 'iteration'])
            sys.stdout.write("\tDone.\n")

            summary.index = np.arange(summary.shape[0])

            def remove_inconsistent(x, thresh=0.8):
                """Remove iterations with fewer events than the threshold
                fraction
                """
                size = x.groupby('iteration').size()
                return x.groupby('iteration').filter(
                    lambda y: len(y) > (thresh * size.mean()))

            sys.stdout.write("Removing iterations that had too few "
                             "events ...\n")
            summary = summary.groupby(['splice_type', 'probability']).apply(
                remove_inconsistent)
            sys.stdout.write("\tDone.\n")

        sys.stdout.write("Writing filtered MISO summary files ...\n")
        csv = '{}/miso_summary_filtered.csv'.format(out_dir)
        summary.to_csv(csv)
        sys.stdout.write("\tWrote {}\n".format(csv))

        if not downsampled:
            sys.stdout.write("Creating ((event_name, splice_type), samples) "
                             "PSI matrix ...\n")
            psi = summary.pivot_table(rows=('event_name', 'splice_type'),
                                      cols='sample_id',
                                      values='miso_posterior_mean')
            psi.to_csv('{}/psi.csv'.format(out_dir))
            sys.stdout.write("\tWrote {}\n".format(csv))

    @staticmethod
    def max_csv(x):
        '''Take the maximum of integers separated by commas

        Parameters
        ----------
        x : str
            Integers separated by columns

        Returns
        -------
        xmax : int
            Maximum of integers

        >>> max_csv("100,200")
        200
        '''
        return max(map(int, x.split(',')))

    @staticmethod
    def min_csv(x):
        '''Take the minimum of integers separated by commas

        Parameters
        ----------
        x : str
            Integers separated by columns

        Returns
        -------
        xmin : int
            Minimum of integers

        >>> min_csv("100,200")
        100
        '''
        return min(map(int, x.split(',')))

    def read_miso_summary(self, filename):
        '''Read a miso summary file and add helpful columns

        Reads a MISO summary file as a pandas dataframe, and adds these columns:
        1. a copy-paste-able genome location at the end, based on the minimum
           mRNA_starts and maximum mRNA_ends. (df.genome_location)
        2. The difference between df.ci_high and df.ci_low (df.ci_diff)
        3. The left and right halves of the confidence interval, e.g. the right
           half is df.ci_high - df.miso_posterior_mean. (df.ci_left_half and
           df.ci_right_half)
        4. The max of the two left and right confidence interval halves
           (df.ci_halves_max)

        Parameters
        ----------
        filename : str
            Full path location of the miso output file

        Returns
        -------
        summary : pandas.DataFrame
            A (n_events, n_columns) dataframe of
        '''
        df = pd.read_table(filename)
        genome_location = pd.DataFrame(
            ['%s:%d-%d' % (chrom, self.min_csv(starts), self.max_csv(stops))
             for chrom, starts, stops in zip(df.chrom,
                                             df.mRNA_starts,
                                             df.mRNA_ends)],
            columns=['genome_location'], index=df.index)
        ci_diff = pd.DataFrame(df.ci_high - df.ci_low, columns=['ci_diff'],
                               index=df.index)
        ci_halves = pd.DataFrame(
            {'ci_left_half': (df.ci_high - df.miso_posterior_mean),
             'ci_right_half': (df.miso_posterior_mean - df.ci_low)},
            index=df.index)
        ci_halves_max = pd.DataFrame(ci_halves.max(axis=1),
                                     columns=['ci_halves_max'])
        return pd.concat([df, genome_location, ci_diff, ci_halves,
                          ci_halves_max], axis=1)

    @staticmethod
    def counts_pair_to_ints(x):
        """Convert a string of isoform and counts to tuples of python integers

        >>> counts_pair_to_ints("(0,0):552")
        ((0,0), 552)
        """
        x = x.split(':')

        isoforms = tuple(map(int, x[0].strip('()').split(',')))
        counts = int(x[1].rstrip(','))
        return isoforms, counts

    def counts_col_to_dict(self, counts):
        """Turn a column of MISO-formatted isoform counts into a dict

        For example, given a column like this:
            (0,0):552,(1,0):449,(1,1):224
        - (0,0): 552 reads were thrown out, not used for either isoform
        - (1,0): 449 reads were unique to isoform0 ("True" or "1" in the first
            position, and the first position is the 0th position)
        - (1,1): 224 reads could have been attributed to either isoform

        If there was an entry with (0,1) then those reads were unique to
        isoform1, which has "True" or "1" in the second position.

        Parameters
        ----------
        counts : pandas.Series
            A (n_events,) column of a pandas dataframe which has the
            per-isoform counts

        Returns
        -------
        d : dict of dicts
            A dictionary of events to a dictionary of isoform counts
        """
        return dict(
            map(self.counts_pair_to_ints, re.findall('\(\d,\d\):\d+,?',
                                                     counts)))


    def filter_miso_summary(self, summary, ci_max=0.5,
                            per_isoform_counts_min=10):
        """Filter a MISO summary on confidence intervals and read depth

        This filters on the maximum confidence interval size and number of
        "junction reads" (reads that are specific to one isoform)

        Parameters
        ----------
        summary : pandas.DataFrame
            A "tall" dataframe of all samples and splice types
        ci_max : float, optional
            Maximum confidence interval size of the percent spliced in value
        per_isoform_counts_min : int, optional

        Returns
        -------
        filtered_summary : pandas.DataFrame

        """
        original_events = summary.shape[0]
        summary = summary.ix[summary.ci_diff <= ci_max]
        after_ci_events = summary.shape[0]
        isoform_counts = pd.DataFrame.from_dict(
            dict(zip(summary.index,
                     summary.counts.map(self.counts_col_to_dict).values)),
            orient='index')

        # Get counts that support only one specific isoform "junction reads"
        specific_isoform_counts = isoform_counts.ix[:, [(0, 1), (1, 0)]].sum(
            axis=1)

        # Filter on at least 10 "junction reads"
        summary = summary.ix[
            specific_isoform_counts >= per_isoform_counts_min]
        after_counts_events = summary.shape[0]

        # Set the index as just the range now that we've filtered everything
        summary.index = np.arange(after_counts_events)

        sys.stdout.write(
            ' {} events removed with poor confidence (ci >{:.2f})\n'
            .format(after_ci_events - original_events, ci_max))
        sys.stdout.write(
            ' {} events removed with too few reads are unique'
            ' to individual isoforms (n < {})\n'.format(
                after_counts_events - after_ci_events,
                per_isoform_counts_min))
        return summary


if __name__ == '__main__':
    try:
        cl = CommandLine()

        CombineMiso(cl.args['glob_command'], cl.args['out_dir'],
                        cl.args['n_progress'], cl.args['ci_max'])
    except Usage, err:
        cl.do_usage_and_die()