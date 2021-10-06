"""program "pandas_TMT_IRS_norm.py"
Reads in TMT summary files from PAW processing and
computes normalizations (within each TMT and between TMT runs).

IRS normalization should be computed now. Aug. 2016, -PW
Converted for Python 3. Aug. 2017, -PW
Modified for PAW protein summary files. -PW 10/4/2017

written by Phil Wilmarth, OHSU, June 2016.
updated for use with PAW pipeline, -PW Oct. 2017.

The MIT License (MIT)

Copyright (c) 2017 Phillip A. Wilmarth and OHSU

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Direct questions to:
Technology & Research Collaborations, Oregon Health & Science University,
Ph: 503-494-8200, FAX: 503-494-4729, Email: techmgmt@ohsu.edu.
"""
# modification history:

# 12/4/2018 -PW
# added support for unused channels (also any channel where SL norm is not right)
#    the unused channels gets adjustments to match experiments; otherwise passed through

# 12/4/2018 -PW
# checks for duplicate sample names and suffixes with capital letters per plex

# 12/4/2018 -PW
# added a log file for console output

# 7/2/2019 -PW
# made TXT file table parsing a little more robust (trim header line, etc.)


import os
import sys
import re
import functools
import csv

from pandas import Series, DataFrame
import pandas as pd
import numpy as np

import PAW_lib

# any globals defined here
VERSION = 'v1.0.1'

class TMT_experiment(object):
    """Each TMT experiment has several things to keep track of."""
    def __init__(self, label):
        self.label = label          # TMT experiment name
        self.number = None          # integer plex number
        self.all_channels = []      # all plex channels
        self.plex = 0               # number of TMT channels
        self.channels = []          # channel designations
        self.sample_key = {}        # channel designation to biological sample dictionary
        self.grand_totals = []      # column totals for each channel
        self.average_total = 0.0    # average of the column totals
        self.target_total = 0.0     # study-wide global normalization target value
        self.loading_factors = []   # sample loading correction factors
        return

    def _snoop(self):
        print('\nlabel:', self.label)
        print('plex:', self.plex)
        print('channels:', self.channels)
        print('sample_key:', self.sample_key)
        print('grand_totals:', self.grand_totals)
        print('average_total:', self.average_total)
        print('target_total:', self.target_total)
        print('loading_factors:', self.loading_factors)

class PAW_TMT_results(object):
    """Data container for PAW TMT protein summary files."""
    def __init__(self):
        self.prefix_lines = []      # any lines before table
        self.suffix_lines = []      # any lines after table
        self.pre_headers = []       # the headers in the row before the headers
        self.headers = []           # the headers from the main header row
        self.frame = None           # pandas DataFrame with main table
        self.num_cols = 0           # number of columns in main table
        self.num_rows = 0           # number of rows in main table
        self.num_prot = 0           # actual number of proteins to normalize
        self.results_file = ''      # original data file's name
        self.TMT_exps = []          # list of TMT_experiment objects for each TMT experiment
        return

    def save_prefix(self, lines):
        """Saves contents of lines before table"""
        lines = [x.rstrip() for x in lines]
        for line in lines:
            if line:
                self.prefix_lines.append(line)
            else:
                self.prefix_lines.append(' ')
        return

    def save_suffix(self, lines):
        """Saves contents of lines after table"""
        lines = [x.rstrip() for x in lines]
        for line in lines:
            if line:
                self.suffix_lines.append(line)
            else:
                self.suffix_lines.append(' ')
        return
    
    def label_duplicates(self, samples):
        """Make TMT samples unique per plex."""
        # make list of suffixes
        suffixes = ['-A', '-B', '-C', '-D', '-E', '-F', '-G', '-H',
                    '-I', '-J', '-K', '-L', '-M', '-N', '-O', '-P']   # up to 16 duplicate times
        
        # find set of duplicate values
        values = list(samples.values())
        dupes = set([item for idx, item in enumerate(values) if item in values[idx+1:]])

        # replace duplicated values with incremental suffixed values
        for dupe in dupes:
            replace = 0
            for key in samples:
                if samples[key] == dupe:
                    samples[key] += suffixes[replace]
                    replace += 1
        return
    
    def get_sample_keys(self):
        """Makes TMT channel to sample key dictionary."""
        # get the indexes for each experiment
        for exp in self.TMT_exps:
            print('experiment:', exp.label, exp.number)
##            exp.all_channels = [head for head in self.headers if ('TotInt_' in head) and (exp.label in head)] # original statement
            # the test below is more specific
            exp.all_channels = [head for head in self.headers if (head.startswith('TotInt_')) and (head.endswith(exp.label))]
            exp.sample_key = {head: self.pre_headers[j] for (j, head) in enumerate(self.headers) if head in exp.all_channels}
            self.label_duplicates(exp.sample_key)
            exp.channels = [x for x in exp.all_channels if
                            'UNUSED' not in exp.sample_key[x].upper() and 'NOTUSED' not in exp.sample_key[x].upper()]
            print('sample key:')
            for x in exp.channels:
                print(x, exp.sample_key[x])
            exp.plex = len(exp.channels)
        return
      
    def load_table(self, file_name):
        """Loads a table from a file into a pandas data frame."""
        self.results_file = file_name
        # read file into memory, remove EOL from lines
        self.contents = open(file_name, 'rt').readlines()
        for i, line in enumerate(self.contents):
            self.contents[i] = '\t'.join([x.strip() if x.strip() else ' ' for x in line.split('\t')])
        
        # find table start and table end indexes
        table_start, table_end = 0, 0
        # start is after header line
        for i, row in enumerate(self.contents):
            if row.startswith('ProtGroup\tCounter'): # header line
                table_start = i + 1
                
                # replace spaces with underscores in header lines
                self.pre_headers = [re.sub(r' ', r'_', x) for x in self.contents[i-1].split('\t')]
                self.headers = [re.sub(r' ', r'_', x) for x in row.strip().split('\t')]

                self.num_cols = len(self.headers)
                labels = [x[len('Total_'):] for x in self.headers if x.startswith('Total_')]
                self.TMT_exps = [TMT_experiment(label) for label in labels] # create the TMT exp containers
                for j, exp in enumerate(self.TMT_exps):
                    exp.number = j + 1
                break
            
        # end is after the actual data
        for i, row in enumerate(self.contents):
            if i > table_start:
                try:
                    # lines after table are "short" in original TXT files. When Excel is used as an
                    # editor to add sample keys and flag contaminants, the text export has padded lines.
                    # can test for a ProtGroup number in the first cell.
                    x = float(row.split('\t')[0])
                except ValueError:
                    table_end = i
                    break
                
        # trap case where there are no lines after the data
        if table_end == 0:
            table_end = i+1     # last row was getting dropped -PW 20210812
        
        # save table length, get sample keys, save prefix and suffix lines
        self.num_rows = table_end - table_start
        self.get_sample_keys()
        self.save_prefix(self.contents[:table_start-2])
        self.save_suffix(self.contents[table_end:])
                   
        # parse columns into dictionaries (header: list) for data frame loading
        table_dict = {header: [] for header in self.headers}
        for line in self.contents[table_start:table_end]:
            for i, cell in enumerate(line.split('\t')[:self.num_cols]): 
                table_dict[self.headers[i]].append(cell)
                
        # cast intensity columns to correct np array data types
        for exp in self.TMT_exps:
            for header in exp.all_channels:
                table_dict[header] = np.array(table_dict[header], dtype='float64')
                
        # make a data frame from the table column dictionary and return it
        self.frame = DataFrame(table_dict, columns=self.headers)
        return

    def find_study_intersection(self, write):
        """Finds the intersection of protein identifications."""
        
        total_proteins = len(self.frame['Filter'][self.frame['Filter'] == ' '])
        contams = len(self.frame) - total_proteins
        for obj in write:
            print('proteins in results file:', len(self.frame), file=obj)
            print('total number of non-contaminant proteins:', total_proteins, file=obj)
            print('contamaminant and decoy count:', contams, file=obj)

        # compute plex-average columns for each TMT experiment and for all
        self.exp_frames = []
        self.frame['Average_all'] = 0.0
        for exp in self.TMT_exps:
            exp_frame = self.frame[exp.all_channels].copy()
            exp_frame['Average_' + exp.label] = exp_frame.mean(axis=1, skipna=True)
            self.frame['Average_' + exp.label] = exp_frame.mean(axis=1, skipna=True)
            self.frame['Average_all'] += exp_frame.mean(axis=1, skipna=True) # accumulate running total
            self.exp_frames.append(exp_frame)
        self.frame['Average_all'] /= len(self.exp_frames) # make running total into average

        # get experiment labels for tagging missing data
        self.exp_labels = [x.label for x in self.TMT_exps]
        self.all = 0

        def _make_missing_flag(row):
            """Helper function to make flag strings."""
            flags = []
            for i, intensity in enumerate(row):
                if intensity == 0.0:
                    flags.append(self.exp_labels[i])
            if flags:
                if len(flags) == len(self.exp_labels):
                    flag = 'missing in all'
                    self.all += 1
                else:
                    flag = 'missing in ' + '+'.join(flags)
            else:
                flag = ' '
            return flag
       
        # get the average intensities for each experiment in separate frame       
        ave_cols = [x for x in self.frame.columns if x.startswith('Average_') and not x.endswith('_all')]
        missing_df = self.frame[ave_cols]
        missing_flags = missing_df.apply(_make_missing_flag, axis=1)
        self.frame.insert(7, 'Missing', missing_flags)

        # individual experiment averages of zero are missing in that experiment
        for exp in self.TMT_exps:
            not_quant = len(self.frame.loc[self.frame['Average_'+exp.label] == 0.0])
            for obj in write:
                print('..%s proteins were not quantifiable in %s' % (not_quant, exp.label), file=obj)
        for obj in write:
            print('..proteins with no intensities experiment-wide: %s' % self.all, file=obj)
        return

    def move_rows_to_bottom(self):
        """Sorts contams and missing reporter ion proteins to bottom of table.
        Also order row by decreasing average intensity."""
        self.frame.loc[self.frame['Filter'] != ' ', 'Average_all'] = 0.0
        self.frame.sort_values(by=['Filter', 'Missing', 'Average_all'], ascending=[True, True, False], inplace=True)

        # set the index
        self.frame.index = list(range(len(self.frame)))

        # get number of rows without the contaminants
        self.num_prot = len(self.frame[self.frame['Filter'] == ' ']) - 1     # dataframe indexing is inclusive
        return

    def compute_loading_norm_factors(self, write):
        """Compute reporter ion columns totals and loading correction factors."""
        tmt_exps = [x for x in self.TMT_exps if x.channels]
        # column totals per TMT
        for exp in tmt_exps:
            exp.grand_totals = self.frame.loc[:self.num_prot, exp.all_channels].sum(axis=0)
            total, number = 0.0, 0
            for label, colsum in exp.grand_totals.iteritems(): # clunky way to skip unused channels
                if label in exp.channels:
                    total += colsum
                    number += 1
            exp.average_total = total / number

        # experiment-wide average and adjustment factors            
        target_total = np.mean([x.average_total for x in tmt_exps])
        for exp in tmt_exps:
            exp.target_total = target_total
            exp.loading_factors = exp.target_total / exp.grand_totals

        # print out factors
        for obj in write:
            print(file=obj)
        for exp in tmt_exps:
            for obj in write:
                print('%s SL Norm Factors:' % exp.label, file=obj)
            for ch, fac in zip(exp.all_channels, exp.loading_factors):
                for obj in write:
                    print('  %s: %0.6f' % (ch, fac), file=obj)
        return

    def add_loading_normalized_columns(self):
        """Add new columns with sample loading corrected intensities."""
        all_int_cols = [x.all_channels for x in self.TMT_exps]
        self.all_int_cols = [item for sublist in all_int_cols for item in sublist]   # get list of all reporter ion columns
        self.work_frame = self.frame.loc[:self.num_prot, self.all_int_cols].copy()    # get sub-frame of just the data

        # adjust unused channel column sums to 
        

        # loop over TMT experiments and correct each channel by single factor
        for exp in [x for x in self.TMT_exps if x.channels]:
            exp.SLNorm_cols = []
            for label, factor in exp.loading_factors.iteritems():
                SLNorm_col = 'SLNorm_' + exp.sample_key[label] + '_' + str(exp.number)
                exp.SLNorm_cols.append(SLNorm_col)
                if label in exp.channels:
                    self.work_frame[SLNorm_col] = self.work_frame[label] * factor
                else:
                    self.work_frame[SLNorm_col] = self.work_frame[label] * (exp.target_total / exp.average_total)

            # compute the average of the pool(s) in each TMT experiment
            pool_list = [x for x in exp.SLNorm_cols if 'POOL' in x.upper()]
            self.work_frame['AvePool_%d' % exp.number] = self.work_frame[pool_list].mean(axis=1)
        return

    def compute_IRS_norm_factors(self):
        """Now do the IRS stuff."""                
        # find the two pooled standard channels and make temporary frame for IRS computations
        ave_pool_list = [x for x in self.work_frame.columns if 'AVEPOOL' in x.upper()]
        temp_frame = self.work_frame.loc[:, ave_pool_list].copy()

        def _geo_mean(row):
            """Computes gemoetric mean, skipping any zeros."""
            net_row = row[row != 0]
            try:
                power = 1 / len(net_row)
            except ZeroDivisionError:
                return 0.0
            return net_row.product() ** power        

        # compute the geometric mean of the pools and the correction factor columns        
        temp_frame['GeoMeanPools'] = temp_frame.apply(_geo_mean, axis=1)
        self.IRSFactors = ['IRS_Fac_' + x for x in ave_pool_list]        
        for i, col in enumerate(self.IRSFactors):
            ave_pool = temp_frame[ave_pool_list[i]]
            temp_frame[col] = np.where(ave_pool > 0.0, temp_frame['GeoMeanPools'] / temp_frame[ave_pool_list[i]], 0.0)

#        temp_frame.to_csv(os.path.join(self.here, 'IRS_dump.txt'), quoting=csv.QUOTE_NONE, index=False)

        # add the columns back to the original data frame
        col_list = ['GeoMeanPools'] + self.IRSFactors
        for col in col_list:
            self.work_frame[col] = temp_frame[col].copy()
        return

    def add_IRS_normalized_columns(self):
        """Add the IRS normalized columns."""
        # loop over TMT experiments and correct each protein by IRS factors
        for j, exp in enumerate([x for x in self.TMT_exps if x.channels]):
            exp.IRSNorm_cols = []
            for i, channel in enumerate(exp.all_channels):
                IRSNorm_col = 'IRSNorm_' + exp.sample_key[channel] + '_' + str(exp.number)
                exp.IRSNorm_cols.append(IRSNorm_col)
                self.work_frame[IRSNorm_col] = self.work_frame[exp.SLNorm_cols[i]] * self.work_frame[self.IRSFactors[j]]                
        return

######## starts here ##############
print('=================================================================')
print(' pandas_TMT_IRS_norm.py, %s, written by Phil Wilmarth 2017-8 ' % VERSION)    
print('=================================================================')

print('select a protein results file with TMT intensities and sample key:')
default_location = os.getcwd()
results_file = PAW_lib.get_file(default_location,
                                [('Text files', '*.txt')],
                                'Select a protein results file with sample keys')
if not results_file:
    sys.exit()

# set up a log file in location where results file is located
log_obj = open(os.path.join(os.path.dirname(results_file), 'PAW_IRS_log.txt'), 'a')

print('\n=================================================================', file=log_obj)
print(' pandas_TMT_IRS_norm.py, %s, written by Phil Wilmarth 2017-8 ' % VERSION, file=log_obj)    
print('=================================================================', file=log_obj)
print('select a protein results file with TMT intensities and sample key:', file=log_obj)

write = [None, log_obj]

# load the results file
for obj in write:
    print('processing:', os.path.split(results_file)[1], file=obj)
frame = PAW_TMT_results()               # create a container for results
frame.load_table(results_file)          # load the results file into container

frame.here = os.path.split(results_file)[0] # keeps track of path

## frame.frame.to_csv(os.path.join(frame.here, 'temp_dump.txt'), sep='\t', quoting=csv.QUOTE_NONE, index=False)

frame.find_study_intersection(write)         # flag proteins not seen in all TMT experiments
frame.move_rows_to_bottom()             # move contaminants and excluded proteins to botton of table
frame.compute_loading_norm_factors(write)    # compute the sample loading correction factors
frame.add_loading_normalized_columns()  # add the sample-loading-normalized columns
frame.compute_IRS_norm_factors()        # add the IRS normalized column factors
frame.add_IRS_normalized_columns()      # finally add the IRS normalized columns

# dump table and see what we have
"""Output matches by hand analysis.
Need to get column header key added to suffix
Need to get proper prefix and suffix framing the table
"""
# need to drop the raw instensity columns from the work_frame before merging
frame.work_frame.drop(frame.all_int_cols, axis = 1, inplace = True)
final = frame.frame.merge(frame.work_frame, how='left', left_index=True, right_index=True)
path, ext = os.path.splitext(results_file)
new_path = path + '_IRS_normalized' + ext
final.to_csv(new_path, sep='\t', quoting=csv.QUOTE_NONE, index=False)

for obj in write:
    print('\nnumber of non-contaminant/decoy proteins:', frame.num_prot + 1, file=obj)
    print('number of proteins with full intensity sets:',
          len(final[(final['Missing'] == ' ') & (final['Filter'] == ' ')]), file=obj)

log_obj.close()
# end
    

    
