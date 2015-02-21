#!/usr/bin/env python

__author__ = 'olga'

import argparse
from glob import iglob
import os
import string
import sys

import numpy as np
import pandas as pd

class CommandLine(object):
    def __init__(self, inOpts=None):
        self.parser = parser = argparse.ArgumentParser(
            description='Combine sailfish output files')
        parser.add_argument('-g', '--glob-command', required=False,
                            default='./*Log.final.out', type=str, action='store',
                            help='Where to find sailfish output directories. '
                                 'Default is folders in the current directory '
                                 'whose names end with "Log.final.out"')
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


class CombineSTARLogFinalOut(object):
    def __init__(self, glob_command, out_dir, n_progress):
        """
        Given a glob command describing where all the Log.final.out files are from
        STAR, return a pd.DataFrame with each sample (id) as its own column.

        @param glob_command: A string that will be passed to glob
        @param ids_function: A function (could be an anonymous function like a
        lambda) that specifies how to get the sample ID from the filename. Could
        also be a list of IDs, but they must be in the exact order as in the
        directories, which is why a function can be easier.

        Example:
        >>> glob_command = '/Users/olga/workspace-git/single_cell/analysis/mapping_stats/*.Log.final.out'
        >>> mapping_stats = self.log_final_out(glob_command,
        ...             lambda x: '_'.join(x.split('/')[-1].split('_')[:2]))
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

        series = []

        filenames = iglob(glob_command)
        n_files = sum(1 for i in filenames)
        sys.stdout.write("Reading {} of STAR's *Log.final.out "
                         "files ...\n".format(n_files))

        # re-initialize iterator
        filenames = iglob(glob_command)

        for i, filename in enumerate(filenames):
            s = pd.read_table(filename, header=None, index_col=0, squeeze=True)
            s.index = s.index.map(
                lambda x: x.rstrip(' |').rstrip(':').rstrip().lstrip())
            converted = [self.maybe_convert_to_float(x.strip('%'))
                         if type(x) != float else x for x in s]
            sample_id = os.path.basename(filename).split('.')[0]
            series.append(pd.Series(converted, index=s.index, name=sample_id))

            if (i + 1) % n_progress == 0:
                sys.stdout.write("\t{}/{} files read\n".format(i + 1, n_files))
        sys.stdout.write("\tDone.\n")

        sys.stdout.write("Merging STAR outputs into a single dataframe...\n")
        mapping_stats = pd.concat(series, axis=1)
        sys.stdout.write("\tDone.\n")


        sys.stdout.write("Adding percentages of splicing events ...\n")
        # Turn all the number of splicing events into percentages for
        # statistical testing
        number_splicing_event_names = ['Number of splices: Annotated (sjdb)',
                                       'Number of splices: GT/AG',
                                       'Number of splices: GC/AG',
                                       'Number of splices: AT/AC',
                                       'Number of splices: Non-canonical']
        percent_splicing_event_names = [x.replace('Number of', '%')
                                        for x in number_splicing_event_names]

        import pdb; pdb.set_trace()
        total_splicing_events = mapping_stats.ix['Number of splices: Total',
                                :].replace(0, np.nan).values.astype(float)

        pieces = []
        for num_events in zip(number_splicing_event_names):
            pieces.append(100.0 * mapping_stats.ix[num_events,
                                  :].values \
                          / total_splicing_events)
        pieces = [np.reshape(piece, len(mapping_stats.columns)) for piece in
                  pieces]
        percent_splicing = pd.DataFrame(pieces,
                                        index=percent_splicing_event_names,
                                        columns=mapping_stats.columns)
        df = pd.concat((mapping_stats, percent_splicing)).T.sort_index()
        sys.stdout.write("\tDone.\n")

        csv = '{}/mapping_stats.csv'.format(out_dir)

        sys.stdout.write("Writing mapping stats ...\n")
        df.to_csv(csv)
        sys.stdout.write("\tWrote {}".format(csv))

    @staticmethod
    def maybe_convert_to_float(x):
        try:
            return float(x)
        except ValueError:
            return x

    @staticmethod
    def maybe_convert_to_int(x):
        try:
            return int(x)
        except ValueError:
            return x

    @staticmethod
    def merge_mapping_stats(s1, s2):
        '''
        Given two series of mapping stats data created within log_final_out,
        merge their mapping stats in an intuitive way: add the raw values,
        and take weighted averages of the percentages.
        '''
        S = pd.Series(index=s1.index, dtype=object)

        # Since this Log.final.out file is a known format, hardcode the field we
        # know to exist, and have to take a weighted average of uniquely or
        # multimapping reads
        weighted_avg_input_reads = ['Mapping speed, Million of reads per hour',
                                    'Average input read length']
        weighted_avg_uniquely_mapped = ['Uniquely mapped reads %',
                                        'Average mapped length',
                                        'Mismatch rate per base, %',
                                        'Deletion rate per base',
                                        'Deletion average length',
                                        'Insertion rate per base',
                                        'Insertion average length']
        weighted_avg_multi_mapped = ['% of reads mapped to multiple loci',
                                     '% of reads mapped to too many loci']
        weighted_avg_un_mapped = ['% of reads unmapped: too many mismatches',
                                  '% of reads unmapped: too short',
                                  '% of reads unmapped: other']

        s1_n_input_reads = int(s1['Number of input reads'])
        s2_n_input_reads = int(s2['Number of input reads'])
        n_input_reads = float(s1_n_input_reads + s2_n_input_reads)

        s1_n_uniquely_mapped = int(s1['Uniquely mapped reads number'])
        s2_n_uniquely_mapped = int(s2['Uniquely mapped reads number'])
        n_uniquely_mapped = float(s1_n_uniquely_mapped + s2_n_uniquely_mapped)

        s1_n_un_mapped = int(
            s1['Number of input reads']) - s1_n_uniquely_mapped
        s2_n_un_mapped = int(
            s2['Number of input reads']) - s2_n_uniquely_mapped
        n_un_mapped = float(s1_n_un_mapped + s2_n_un_mapped)
        for ind in S.index:
            # print 'ind: "%s"' % ind
            stat1 = s1[ind]
            stat2 = s2[ind]

            try:
                stat1_float = np.float(stat1)
            except ValueError:
                S[ind] = stat1
                continue
            if np.isnan(stat1_float):
                S[ind] = str(stat1)
            elif ind in weighted_avg_input_reads:
                stat1 = float(stat1)
                stat2 = float(stat2)
                S[ind] = (
                             s1_n_input_reads * stat1 + s2_n_input_reads * stat2) / n_input_reads
            elif ind in weighted_avg_uniquely_mapped:
                stat1 = float(stat1)
                stat2 = float(stat2)
                S[ind] = (
                             s1_n_uniquely_mapped * stat1 + s2_n_uniquely_mapped * stat2) / n_uniquely_mapped
            elif ind in weighted_avg_multi_mapped:
                stat1 = float(stat1)
                stat2 = float(stat2)
                S[ind] = (
                             s1_n_un_mapped * stat1 + s2_n_un_mapped * stat2) / n_un_mapped
            elif ind in weighted_avg_un_mapped:
                stat1 = float(stat1)
                stat2 = float(stat2)
                S[ind] = (
                             s1_n_un_mapped * stat1 + s2_n_un_mapped * stat2) / n_un_mapped
            else:
                stat1 = int(stat1)
                stat2 = int(stat2)
                S[ind] = stat1 + stat2
                # print "S['Number of splices: Total']", S['Number of splices: Total']
        return S


    def fix_duplicate_columns(self, mapping_stats):
        '''
        Given a dataframe of mapping stats with "duplicate" column names
        (actually duplicate column names aren't allowed in pandas,
        but this assumes you have columns like M1_01 and M1_01a, and the M1_01a
        column was created from the second *.Log.final.out file from RNA-STAR),
        this detects the duplicate columns and sends them to merge_mapping_stats,
        which combines the duplicate colmns into single column in a reasonable way
        '''
        duplicate_columns = mapping_stats.filter(regex='[a-z]$')
        mapping_stats.drop((x[:-1] for x in duplicate_columns), axis=1)
        merged_columns = {}
        for col in duplicate_columns:
            original_col = col[:-1]
            merged_columns[original_col] = self.merge_mapping_stats(
                mapping_stats.ix[:, col], mapping_stats.ix[:, original_col])

        # Remove duplicately-named columns, e.g. M1_01 and M1_01a
        mapping_stats = mapping_stats.drop((x for x in duplicate_columns),
                                           axis=1)
        mapping_stats = mapping_stats.drop((x[:-1] for x in duplicate_columns),
                                           axis=1)

        merged_df = pd.DataFrame(data=merged_columns,
                                 index=mapping_stats.index)
        # print 'merged_df.index', merged_df.index
        # print " merged_df.ix['Number of splices: Total', :]"
        # print merged_df.ix['Number of splices: Total', :]
        # print 'merged_df.columns', merged_df.columns

        #
        # print 'mapping_stats.index', mapping_stats.index
        # print 'mapping_stats.columns', mapping_stats.columns

        mapping_stats = pd.concat((mapping_stats, merged_df), axis=1)

        # print "mapping_stats.ix['Number of splices: Total','P3_01']",
        # print mapping_stats.ix['Number of splices: Total','P3_01']

        return mapping_stats


    def make_unique(seq, idfun=None):
        '''
        if an object appears more than once in a list, append a letter to it

        Modified from: http://www.peterbe.com/plog/uniqifiers-benchmark
        '''
        if idfun is None:
            def idfun(x): return x
        seen = {}
        result = []
        for item in seq:
            marker = idfun(item)
            # in old Python versions:
            # if seen.has_key(marker)
            # but in new ones:
            if marker in seen:
                seen[marker] += 1
                result.append(item + string.lowercase[seen[marker] - 2])
                continue
            seen[marker] = 1
            result.append(item)
        return result




if __name__ == '__main__':
    try:
        cl = CommandLine()

        CombineSTARLogFinalOut(cl.args['glob_command'], cl.args['out_dir'],
                        cl.args['n_progress'])
    except Usage, err:
        cl.do_usage_and_die()