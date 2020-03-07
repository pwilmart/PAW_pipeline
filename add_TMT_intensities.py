"""program "add_TMT_intensities.py"
Adds TMT intensities to PAW protein results files.
Written by Phil Wilmarth, OHSU, Sept. 2017.

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
# Requires Python 3.6 or greater and uses numpy.

import os
import sys
import re

from pandas import Series, DataFrame
import numpy as np

import PAW_lib

# globals and constants
# SPS MS3 zero replacement values
MS3_INTENSITY = 500.0    # MS3: individual PSM test of minimum trimmed average reporter ion intensity
MS3_MISSING = 150.0      # MS3: replacement for zero total intensities only at the summed protein level
# MS2 reporter ion zero replacement values
MS2_INTENSITY = 2000.0   # MS2: individual PSM test of minimum trimmed average reporter ion intensity
MS2_MISSING = 500.0      # MS2: replacement for zero total intensities only at the summed protein level
# version
VERSION = 'v1.0.1'

# output labels for different TMT sets
TMT6 = ['TotInt_126', 'TotInt_127', 'TotInt_128',
        'TotInt_129', 'TotInt_130N', 'TotInt_131']
TMT10 = ['TotInt_126C', 'TotInt_127N', 'TotInt_127C', 'TotInt_128N', 'TotInt_128C',
       'TotInt_129N', 'TotInt_129C', 'TotInt_130N', 'TotInt_130C', 'TotInt_131N']
TMT11 = ['TotInt_126C', 'TotInt_127N', 'TotInt_127C', 'TotInt_128N', 'TotInt_128C',
       'TotInt_129N', 'TotInt_129C', 'TotInt_130N', 'TotInt_130C', 'TotInt_131N', 'TotInt_131C']
TMT16 = ['TotInt_126C', 'TotInt_127N', 'TotInt_127C', 'TotInt_128N', 'TotInt_128C',
         'TotInt_129N', 'TotInt_129C', 'TotInt_130N', 'TotInt_130C', 'TotInt_131N',
         'TotInt_131C', 'TotInt_132N', 'TotInt_132C', 'TotInt_133N', 'TotInt_133C', 'TotInt_134N']

def base_peptide_sequence(sequence, mask=False):
    """Returns the peptide amino acid residues from SEQUEST peptide strings
    """
    # remove bounding residues (SEQUEST/Comet format: A.BCD.E)
    prefix, peptide, suffix = split_peptide(sequence)

    # remove the 2017 Comet style mod strings
    peptide = re.sub(r'\[[-+]?[0-9]*(.)?[0-9]*\]', '', peptide)        
    # remove modification symbols: '*', '#', '@', '^', '~', '$', '%', '!', '+', 'n', 'c', '[', ']', "(', ')', '{', '}'
    peptide = re.sub(r'[*#@^~$%!+nc\[\]\{\}\(\)]', '', peptide)
    
    # mask I/L if needed:
    if mask:
        return re.sub(r'[IL]', 'j', peptide)
    else:
        return peptide

def split_peptide(sequence):
    """Splits peptide assuming that there might be single preceeding and following residues with periods."""
    if re.match(r'[-A-Z]\..+\.[-A-Z]', sequence):
        return sequence[0], sequence[2:-2], sequence[-1]
    else:
        if min(sequence.count('['), sequence.count(']')) != sequence.count('.'):
            print('   WARNING: possible malformed peptide string:', sequence)
        return '', sequence, ''

class PAWTable(object):
    """General container for PAW tabular results files.
    Assumes that table may have information before (prefix) and after (suffix)
    the table. Assumes that there is one header row with some unique start string.
    The number of table columns is assumed to be the number of tab-delimited
    elements in the header row. The table can have a ragged last column. Any
    optional lines after the table rows are assumed to be shorter than
    the table lines.
    """
    def __init__(self, file_path, header_test_str):
        self.prefix = []    # any lines in file before header line
        self.headers = []    # the table header line
        self.col_map = {}   # maps header to column number
        self.num_cols = 0   # number of elements in header line
        self.table = []     # the actual table lines (not parsed)
        self.num_rows = 0   # the number of lines in the table
        self.suffix = []    # any (optional) lines after the table

        # read in the file (strip out any Excel CSV stuff)
        try:
            contents = open(file_path, 'rt').readlines()
            for i, line in enumerate(contents):
                items = [x.strip() if x.strip() else ' ' for x in line.split('\t')]
                contents[i] = '\t'.join(items)
        except:
            print('...Could not read in file:', file_path)
            raise

        # find the header line
        start, end = 0, 0
        for i, line in enumerate(contents):
            if line.startswith(header_test_str):
                start = i + 1
                self.headers = line.split('\t')
                self.col_map = {v: i for i, v in enumerate(self.headers)}
                self.num_cols = len(self.col_map)
                self.prefix = contents[:i]
                break

        # assume any rows after main table have fewer columns than header
        for i, line in enumerate(contents):
            if i < start:
                continue
            if len(line.split('\t')) < (self.num_cols - 1):   # allows for a ragged last column
                end = i
                break
            

        # if table had no extra rows then end=0
        if end == 0:
            self.table = contents[start:]
        else:
            self.table = contents[start:end]
            self.suffix = contents[end:]

        # add table length
        self.num_rows = len(self.table)

        return

class PAWPeptideSummary(PAWTable):
    """Container for PAW peptide summary files.
    Has methods to make sets of peptides keyed by protein group
    number and return a set of peptides given a group number.
    """
    def __init__(self, file_path, header_test_str, write):
        """See PAWTable class for more details on attributes."""
        PAWTable.__init__(self, file_path, header_test_str)
        self.write = write
        self.peptide_file = file_path
        self.peptide_sets_dict = {}
        self.make_peptide_sets_dict()

    def make_peptide_sets_dict(self):
        """Makes dictionary of peptide sequence sets keyed by protein group number.
        Need to have peptide sequence keys match peptide table rows, i.e. separated by charge."""
        for line in self.table:
            items = line.split('\t')
            key = int(items[self.col_map['ProtGroup']])
            value = items[self.col_map['Sequence']].split('.')[1] + '_' + items[self.col_map['Z']]
            if key in self.peptide_sets_dict:
                self.peptide_sets_dict[key].add(value)
            else:
                self.peptide_sets_dict[key] = set([value])
        for obj in self.write:
            print('length peptide_sets_dict is:', len(self.peptide_sets_dict), file=obj) 

    def get_peptide_set(self, group_number):
        """Returns peptide set give a protein group number."""
        return self.peptide_sets_dict[group_number]

class PeptideSummaryRow(object):
    """Make a container for original line and TMT intensities summed to the peptide level."""
    def __init__(self, line, col_map):
        self.original_line = line
        self.col_map = col_map
        self.intensities = None
        self.new_line = None
        # make the peptide sequence key
        items = line.split('\t')
        self.seq_key = items[col_map['Sequence']].split('.')[1] + '_' + items[col_map['Z']]

    def load_intensities(self, psm_lists_dict_list, tmt_intensity_dict, number_channels):
        """Get key information from line and fetch TMT intensities for each PSM.
        We need col_map to parse line,
        we need the psm_list_dict to get psm lists for each peptide sequence,
        we need the TMT intensity dictionary keyed by psm.
        Computes and adds the summed intensities.
        """
        if len(self.original_line.split('\t')) == (len(self.col_map) - 1):
            self.original_line += '\t'
        self.new_line = self.original_line
        
        for psm_lists_dict in psm_lists_dict_list:
            int_list = []
                
            # loop over all psms in the psm list and collect intensities
            try:
                for psm in psm_lists_dict[self.seq_key]:
                    int_list.append(tmt_intensity_dict[psm])
            except KeyError:
                pass
                
            # make 2D array of reporter ion intensities and sum over psms
            if len(int_list) > 1:
                int_array = np.array(int_list)
                self.intensities = int_array.sum(axis=0)
            elif len(int_list) == 1:
                self.intensities = np.array(int_list[0])
            else:
                self.intensities = np.zeros(number_channels)

            # add peptide reporter ion intensities to peptide summary line
            self.new_line = self.new_line + '\t' + '\t'.join([str(x) for x in self.intensities])

class ProteinSummaryRow(object):
    """Make a container for original line and TMT intensities summed to the protein level."""
    def __init__(self, line, col_map, minimum_intensity=500.0, missing_intensity=50.0):
        self.original_line = line
        self.col_map = col_map
        self.minimum_intensity = minimum_intensity
        self.missing_intensity = missing_intensity
        self.intensities = None
        self.new_line = None
        # make the protein group key
        items = line.split('\t')
        self.group_key = int(float(items[col_map['ProtGroup']]))
                             
    def load_intensities(self, peptide_sets_dict, tmt_intensity_dict,
                         unique_peptide_dict, psm_lists_dict_list, number_channels):
        """Get key information from line and fetch TMT intensities for each PSM.
        We need col_map to parse line,
        we need the psm_list_dict to get psm lists for each peptide sequence,
        we need the TMT intensity dictionary keyed by psm.
        Computes and adds the summed intensities.
        """
        if len(self.original_line.split('\t')) == (len(self.col_map) - 1):
            self.original_line += '\t'
        self.new_line = self.original_line
        # loop over all peptides in the peptide set and collect intensities
        protein_psm_set = set()
        try:
            unique_set = [x for x in peptide_sets_dict[self.group_key] if x in unique_peptide_dict]
        except KeyError:
            print('KEYERROR')
            print(self.group_key)
            unique_set = []
        for psm_lists_dict in psm_lists_dict_list:
            int_list = []
            for pep in unique_set:
                
                # loop over all psms in the psm list and collect intensities (just once per psm!)
                try:
                    for psm in psm_lists_dict[pep]:
                        if psm in protein_psm_set:
                            continue
                        else:
                            protein_psm_set.add(psm)
                        int_list.append(tmt_intensity_dict[psm])
                except KeyError:
                    pass

            # make 2D array of reporter ion intensities and sum over psms
            if len(int_list) > 1:
                int_array = np.array(int_list)
                self.intensities = int_array.sum(axis=0)
                self.trimmed_test(number_channels)
            elif len(int_list) == 1:
                self.intensities = np.array(int_list[0])
                self.trimmed_test(number_channels)
            else:
                self.intensities = np.zeros(number_channels)

            # add peptide reporter ion intensities to protein summary line
            middle = '\t%d\t' % len(int_list)
            self.new_line = self.new_line + middle + '\t'.join([str(x) for x in self.intensities])

    def trimmed_test(self, number_channels):
        """Finds average intensity of channels, excluding the top and bottom values."""
        int_vector = sorted(list(self.intensities))

        # compute the trimmed average
        average = sum(int_vector[1:-1])/float(len(int_vector[1:-1]))
            
        # test threshold
        if average >= self.minimum_intensity:
            self.intensities[self.intensities == 0.0] = self.missing_intensity
        else:
            self.intensities = np.zeros(number_channels)
        return
                        
class PAWProteinSummary(object):
    """Data container for PAW protein and grouped protein summary files."""
    def __init__(self):
        self.protein_file = ''  # protein results file name
        self.peptide_file = ''  # peptide results file name
        self.results_path = ''  # results folder path
        self.project_path = ''  # main project folder (should be one level above results)
        self.write = [None]     # setting up for log file writing
        self.table = []         # the actual table lines
        self.new_table = []     # non-redundant table (if plain report)
        self.prefix = []        # any lines before table
        self.suffix = []        # any lines after table
        self.pre_headers = []   # the row before the headers
        self.headers = []       # original headers
        self.col_map = {}       # column header to index mapping
        self.new_headers = []   # headers with no spaces and label suffixed
        self.frame = DataFrame()# main data table
        self.num_rows = 0       # number of rows in main table
        self.num_cols = 0       # number of columns in main table

        # add some dictionaries for tacking peptide status, TMT scans, etc.
        self.tmt_intensity_dict = {}    # maps TMT intensities to short DTA names (psms)
        self.unique_peptide_dict = {}   # maps uniqueness to peptide sequences
        self.psm_lists_dict_list = []   # maps lists of short DTA names (psms) to peptide sequences
        
        self.average_int = 0.0  # within TMT Total Normalization target
        return
        
    def make_new_headers(self):
        """Replaces spaces with underscores in headers."""
        self.new_headers = [re.sub(r' ', r'_', x) for x in self.headers] # changes spaces to underscores
        return
        
    def load_table(self, file_name):
        """Loads a results file into a various data structures.
        This is kind of general allowing content both before and after the table.
        Some preparations of data are done with regular Python before loading the
        table into a pandas dataframe."""
        self.protein_file = file_name
        self.results_path = os.path.dirname(self.protein_file)
        self.project_path = os.path.dirname(self.results_path)

        # read file contents, remove any Excel CSV stuff
        contents = open(self.protein_file, 'rt').readlines()
        for i, line in enumerate(contents):
            contents[i] = '\t'.join([x.strip() if x.strip() else ' ' for x in line.split('\t')])
                
        # find table start row number and process header line
        table_start, table_end = 0, 0
        for i, row in enumerate(contents):
            if row.startswith('ProtGroup\tCounter\tAccession'): # header line
                table_start = i+1

                # process header line (needs some fixing to get unique headers)
                self.headers = row.rstrip().split('\t')
                self.num_cols = len(self.headers)
                self.pre_headers = contents[i-1].rstrip().split('\t')    # get items from line before main header line
                self.make_new_headers()
                self.col_map = {v: i for i, v in enumerate(self.new_headers)}
                break

        # find table end row number
        for i, row in enumerate(contents[table_start:]):
            try:
                group = float(row.split('\t')[0]) # table rows start with group numbers
            except:
                group = None
            if (len(row.split('\t')) < (self.num_cols-1)) and not group:
                table_end = i + table_start - 1
                break
        if table_end == 0:
            table_end = i + table_start
        
        # save table, prefix, and suffix lines
        self.num_rows = table_end - table_start
        self.table = [x.rstrip() for x in contents[table_start:table_end+1]]
        self.prefix = [x.rstrip() for x in contents[:table_start-2]]
        self.suffix = [x.rstrip() for x in contents[table_end+1:]]

        # get the list of biological samples
        self.sample_list = [x.replace("Total_", "") for x in self.headers if x.startswith("Total_")]

        # if regular protein summary, make the table non-redundant by row
        if self.new_headers[-1] == 'OtherLoci':
            self.make_nr_table()
        else:
            self.new_table = self.table
                   
        # parse columns into dictionaries (header: list) for data frame loading
        table_dict = {header: [] for header in self.new_headers}
        for line in self.new_table:
            cells = line.split('\t')
            if len(cells) < self.num_cols:
                cells.append('')
            for i, cell in enumerate(cells): 
                table_dict[self.new_headers[i]].append(cell)
                
        # cast numerical columns to correct np array data types
        # most columns can remain as text, only convert a few to numbers.
        floats = ['ProtGroup', 'UniqFrac'] + ['Corrected_' + x for x in self.sample_list]
        for key in floats:
            table_dict[key] = np.array(table_dict[key], dtype='float64')
        ints = (['Counter', 'CountsTot', 'UniqueTot'] +
                ['Total_' + x for x in self.sample_list] +
                ['Unique_' + x for x in self.sample_list])
        for key in ints:
            table_dict[key] = np.array(table_dict[key], dtype='int32')
                
        # make a data frame from the table column dictionary and return it
        self.frame = DataFrame(table_dict, columns=self.new_headers)
        self.frame.index = self.frame.Accession

    def make_nr_table(self):
        """Makes redundant protein summaries into 'one-row-per-protein-group' tables.
        This applies to regular protein summaries not the grouped protein summaries."""
        # take care of adding and moving columns first
        # Add column for identicals
        self.new_headers.insert(self.col_map['Filter'], 'Identical')
        for i, line in enumerate(self.table):
            cells = line.split('\t')
            if len(cells) < self.num_cols:  # pad table if needed
                cells.append(' ')
            cells.insert(self.col_map['Filter'], ' ')
            self.table[i] = '\t'.join(cells)

        # update column mapping
        self.col_map = {v: i for i, v in enumerate(self.new_headers)}

        # move 'OtherLoci' column
        self.new_headers.insert(self.col_map['Filter'], self.headers[-1])
        for i, line in enumerate(self.table):
            cells = line.split('\t')
            cells.insert(self.col_map['Filter'], cells[-1])
            self.table[i] = '\t'.join(cells[:-1])

        # update new headers and column mapping
        self.new_headers = self.new_headers[:-1]
        self.col_map = {v: i for i, v in enumerate(self.new_headers)}

        # gather up accessions for each protein group
        idx = self.col_map['ProtGroup']
        groups = {}
        for line in self.table:
            cells = line.split('\t')
            group_key = int(float(cells[idx]))
            if group_key in groups:
                groups[group_key].append(cells[self.col_map['Accession']])
            else:
                groups[group_key] = [cells[self.col_map['Accession']]]

        # copy primary ProtGroup to new table and add "Identicals"
        self.new_table = []
        for line in self.table:
            cells = line.split('\t')
            if float(cells[idx]).is_integer():
                key = int(float(cells[idx]))
                if len(groups[key]) > 1:
                    cells[self.col_map['Accession']] += (' (+%d)' % len(groups[key]))
                    cells[self.col_map['Identical']] = '&'.join(groups[key])
                self.new_table.append('\t'.join(cells))

    def load_TMT_intensities(self, minimum_intensity=500.0, missing_intensity=0.0):
        """Loads PSM TMT intensities."""
        total = set()
        reject = set()
        txt_list = None
        # get the TMT intensities from the filtered SQT/TXT files
        for obj in self.write:
            print('getting intensities from SQT/TXT files', file=obj)
        filtered_folder = os.path.join(self.project_path, 'filtered_files')
        # first: look for filtered SQT file names
        sqt_list = [x for x in os.listdir(filtered_folder) if x.endswith('_filtered.sqt')]
        if sqt_list:
            txt_list = [x.replace('_filtered.sqt', '_filtered.txt') for x in sqt_list]
        # second: look for *.PAW.txt files
        if not txt_list:
            txt_list = [x for x in os.listdir(filtered_folder) if x.endswith('.PAW.txt')]
        # finally: have user select the TXT files
        if not txt_list:
            ext_list = [('TXT files', '*.txt'), ('All files', '*.*')]
            title = 'Select one or more PAW text files'
            txt_list = PAW_lib.get_files(filtered_folder, ext_list, title)
            if not txt_list:
                sys.exit()
            txt_list = [os.path.split(x)[1] for x in txt_list] # need just the filename from path
        for txt_file in txt_list:
            lc_name = txt_file.replace('.txt', '')
            lc_name = os.path.splitext(txt_file)[0]
            lc_total = set()
            lc_reject = set()
            print_all = True
            with open(os.path.join(filtered_folder, txt_file), 'r') as txt_obj:
                for line in txt_obj:
                    line = line.rstrip()
                    if line.startswith('start\tend'):
                        col_map = {v: i for i, v in enumerate(line.split('\t'))}
                        heights = [x for x in line.split('\t') if x.startswith('height_')]
                        continue
                    intensities = []
                    items = line.split('\t')
                    key = lc_name + '.' + items[col_map['start']]
                    total.add(key)
                    lc_total.add(key)
                    for height in heights:
                        intensities.append(float(items[col_map[height]]))
                    # test minimum trimmed intensities against cutoff
                    test = sorted(intensities)[1:-1]
                    try:
                        if (sum(test) / len(test)) < minimum_intensity:
                            intensities = [0.0 for x in intensities]    # zero out low intensity PSMs
                            reject.add(key)
                            lc_reject.add(key)
                        else:
                            intensities = [x if x > 0.0 else missing_intensity for x in intensities] # option to replace zero
                    except ZeroDivisionError:
                        print('Divide by zero exception code:')
                        print('..key:', key)
                        print('..test array:', test)
                        intensities = [0.0 for x in heights]
                        print('..intensities:', intensities)
                    # update dictionary
                    if key in self.tmt_intensity_dict:
                        continue    # scans are redundant, one for each loci
                    else:
                        self.tmt_intensity_dict[key] = np.array(intensities)
                        
            # print LC run reject rate
            if print_all:
                print('\nLC run:', txt_file)
                print('Total: %d, reject: %d, reject rate: %0.2f'
                      % (len(lc_total), len(lc_reject), 100.0 * len(lc_reject)/len(lc_total)))
                            
        self.number_channels = len(heights)
        if print_all:
            print()
        for obj in self.write:
            print('reporter ions total:', len(total), file=obj)
            print('below intensity cutoff:', len(reject), file=obj)
##            print('   intensity cutoff was:', minimum_intensity, file=obj)
            print('length of tmt_intensity_dict:', len(self.tmt_intensity_dict), file=obj)

    def load_peptide_unique_dict(self):
        """Reads the peptide summary file and saves unique peptide sequences in dictionary."""
        # look for the peptide file that matches the protein file
        self.peptide_file = self.protein_file.replace('protein_summary_9', 'peptide_summary_9')
        if (self.peptide_file == self.protein_file) or not os.path.exists(self.peptide_file):
            print('peptide file missing:', self.protein_file)
            self.peptide_file = ''
            sys.exit()

        # count all peptides
        peptide_dict = {}
        nr_peptides = {}

        # open the peptide file and save the unique peptides
        with open(self.peptide_file) as fin:
            process = False
            for i, line in enumerate(fin):
                line = line.rstrip()

                # look for header line
                if line.startswith('ProtGroup\tAccession'):
                    cols = {v: i for i, v in enumerate(line.split('\t'))}
                    process = True
                    continue
                
                # process main table lines
                if process and len(line.split('\t')) > 3:
                    items = line.split('\t')
                    sequence = items[cols['Sequence']]
                    if sequence == 'X.X.X':
                        continue
                    base_pep = base_peptide_sequence(sequence)
                    nr_peptides[base_pep] = True
                    seq = sequence.split('.')[1] + '_' + items[cols['Z']]
                    peptide_dict[seq] = True
                    if items[cols['Unique']] == 'TRUE':
                        self.unique_peptide_dict[seq] = True
                            
        for obj in self.write:
            print('length of peptide_dict:', len(peptide_dict), file=obj)
            print('length of unique_peptide_dict:', len(self.unique_peptide_dict), file=obj)
            print('non-redundant peptide sequence count:', len(nr_peptides), file=obj)

    def load_psm_list_dict(self):
        """Loads the lists of psm DTA-like names for each peptide sequence
        for each biological sample."""
        for sample in self.sample_list:
            psm_lists_dict = {}
            nr_psm_lists_dict = {}
            
            # look for the detailed peptide summary file
            peptide_file = os.path.join(self.results_path, sample + '_peptide_results_9.txt')
            if not os.path.exists(peptide_file):
                print('detailed peptide file issue in:', self.results_path)
                print('peptide_file:', peptide_file)
                sys.exit()

            # read in the detailed peptide file and build the psm lists
            with open(peptide_file, mode='rt') as fin:
                process = False
                for line in fin:
                    line = line.rstrip()

                    # look for header line
                    if line.startswith('ProtGroup\tAccession'):
                        cols = {v: i for i, v in enumerate(line.split('\t'))}
                        process = True
                        continue

                    # process the main lines
                    if process and len(line.split('\t')) > 3:
                        items = line.split('\t')
                        seq = items[cols['Sequence']].split('.')[1] + '_' + items[cols['Z']]
                        dta_short = '.'.join(items[cols['DTA_filename']].split('.')[:-2])
                        if seq == 'X.X.X':
                            continue
                        if seq in psm_lists_dict:
                            psm_lists_dict[seq].append(dta_short)
                        else:
                            psm_lists_dict[seq] = [dta_short]

            # psm lists will have duplicates due to redundancy
            for k in psm_lists_dict:
                nr_psm_lists_dict[k] = set(psm_lists_dict[k])

            # add the dictionary to the list 
            self.psm_lists_dict_list.append(nr_psm_lists_dict)

            # print some summary stats
            both = [len(nr_psm_lists_dict[k]) for k in nr_psm_lists_dict]
            unique = [len(nr_psm_lists_dict[k]) for k in nr_psm_lists_dict if k in self.unique_peptide_dict]
            for obj in self.write:
                print(os.path.split(peptide_file)[1], file=obj)
                print('..length psm_list_dict:', len(nr_psm_lists_dict), file=obj)
                print('..sum of all peptides:', sum(both), file=obj)
                print('..sum of all unique peptides:', sum(unique), file=obj)
            
def make_new_prefix(protein_summary):
    """Pads the prefix for headers are in Row 5 and adds sample labels above TMT reporter ion channels."""
    if len(protein_summary.prefix) < 4:
        for i in range(4-len(protein_summary.prefix)):
            protein_summary.prefix.append(' ')

###############################################################################                        
############################# program starts here ############################# 

print('Select a PAW protein results file')
default_location = r'F:\PSR_Core_Analysis'
if not os.path.exists(default_location):
    default_location = os.getcwd()
results_file = PAW_lib.get_file(default_location,
                                [('Text files', '*.txt')],
                                'Select a protein results file')
if not results_file:
    sys.exit()

# get type of reporter ions (MS2 or MS3) and set test and replacement values
if PAW_lib.get_yesno('Reporter ions', 'Is this SPS MS3 data?'):
    INTENSITY = MS3_INTENSITY
    MISSING = MS3_MISSING
else:
    INTENSITY = MS2_INTENSITY
    MISSING = MS2_MISSING
    
# get location and set up log file
# not sure about platform default encoding between PC and Mac
log_obj = open(os.path.join(os.path.dirname(results_file), 'add_TMT_intensities.log'), mode='at')
write = [None, log_obj]

# print name, version information
for obj in write:
    print(file=obj)
    print('===============================================================', file=obj)
    print(' add_TMT_intensities.py, %s, written by Phil Wilmarth, OHSU' % VERSION, file=obj)
    print('===============================================================', file=obj)

# read in the protein summary file
for obj in write:
    print('results file:', os.path.split(results_file)[1], file=obj)
    print('trimmed average intensity must exceed:', INTENSITY, file=obj)
    print('final zero values replaced with:', MISSING, file=obj)
protein_summary = PAWProteinSummary()
protein_summary.write = write
protein_summary.load_table(results_file)

# compute various dictionaries
protein_summary.load_TMT_intensities(INTENSITY)
protein_summary.load_peptide_unique_dict()
protein_summary.load_psm_list_dict()    # this is a list of dictionaries, one for each sample

# read in the peptide summary file
peptide_summary = PAWPeptideSummary(protein_summary.peptide_file, 'ProtGroup\tAccession', write)

# get summed reporter ions for each peptide
for i, line in enumerate(peptide_summary.table):
    row_obj = PeptideSummaryRow(line, peptide_summary.col_map)
    row_obj.load_intensities(protein_summary.psm_lists_dict_list, protein_summary.tmt_intensity_dict,
                             protein_summary.number_channels)
    peptide_summary.table[i] = row_obj.new_line

# now do the total protein intensities (unique peptides only)
for i, line in enumerate(protein_summary.new_table):
    prot_obj = ProteinSummaryRow(line, protein_summary.col_map, INTENSITY, MISSING)
    prot_obj.load_intensities(peptide_summary.peptide_sets_dict, protein_summary.tmt_intensity_dict,
                              protein_summary.unique_peptide_dict, protein_summary.psm_lists_dict_list,
                              protein_summary.number_channels)
    protein_summary.new_table[i] = prot_obj.new_line

# write out the new protein summary file
if protein_summary.number_channels == 6:
    TMT = TMT6
elif protein_summary.number_channels == 10:
    TMT = TMT10
elif protein_summary.number_channels == 11:
    TMT = TMT11
elif protein_summary.number_channels == 16:
    TMT = TMT16
new_file = protein_summary.protein_file.replace('summary_9', 'summary_TMT_9')
if new_file == protein_summary.protein_file:
    print('error creating new protein file name')
    sys.exit()
tmt_headers = []
for sample in protein_summary.sample_list:
    tmt_headers += ['PSMs_Used_'+sample] + [x+'_'+sample for x in TMT]
    
with open(os.path.join(protein_summary.results_path, new_file), 'w') as fout:
    protein_summary.prefix.insert(0, 'New PAW TMT results file!')
    make_new_prefix(protein_summary)
    for line in protein_summary.prefix:
        if not line:
            line = ' '
        print(line, file=fout)
    print('\t'.join(protein_summary.new_headers + tmt_headers), file=fout)
    for line in protein_summary.new_table:
        print(line, file=fout)
    for line in protein_summary.suffix:
        if not line:
            line = ' '
        print(line, file=fout)

# write out the new peptide summary file
new_file = peptide_summary.peptide_file.replace('summary_9', 'summary_TMT_9')
if new_file == peptide_summary.peptide_file:
    print('error creating new protein file name')
    sys.exit()
tmt_headers = []
for sample in protein_summary.sample_list:
    tmt_headers += [x+'_'+sample for x in TMT]
with open(os.path.join(protein_summary.results_path, new_file), 'w') as fout:
    peptide_summary.prefix.insert(0, 'New PAW TMT results file!')
    for line in peptide_summary.prefix:
        print(line, file=fout)
    print('\t'.join(peptide_summary.headers + tmt_headers), file=fout)
    for line in peptide_summary.table:
        print(line, file=fout)
    for line in peptide_summary.suffix:
        print(line, file=fout)

log_obj.close()
# end


    
                             

    
