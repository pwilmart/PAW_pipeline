"""PAW_lib.py: Written by Phil Wilmarth and Billy Rathje, OHSU.
Library of support functions and classes for PAW pipeline programs.

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
# Updated for Python 3, Aug. 2017 -PW
# added a sample name parsing class, -PW 10/21/2017
# recoded older routines like "amino_acid_count" for better Comet support. -PW 10/27/2017

import os
import sys
import glob
import fnmatch
import gzip
import time
import re
import copy
from pprint import pprint
from collections import OrderedDict

import tkinter
from tkinter import filedialog

import pandas as pd
import numpy as np

# this is only used in a debugging block
##import matplotlib.pyplot as pyplot    # this seems incompatible with  standard IDLE (OK with Anaconda)

###################### standard dialog boxes ###########################
# updated from fasta_lib.py -PW 9/16/2017

def get_folder(default_location, title_string=None):
    """Dialog box to browse to a folder.  Returns folder path.

    Usage: full_folder_name = get_folder(default_location, [title]),
        where "default_location" is a starting folder location,
        "title" is an optional message to list in the dialog box,
        and "full_folder_name" is the complete selected folder name.
    Written by Phil Wilmarth, 2008, 2016
    """
    # set up GUI elements
    root = tkinter.Tk()
    root.withdraw()
    try:
        root.tk.call('console', 'hide')
    except:
        pass
    
    # set default title string and location if not passed
    if title_string is None:   
        title_string = 'Select a folder with desired files/dirs'
    if not default_location:
        default_location = os.getcwd()
    
    # create dialog box for folder selection
    root.update()   # helps make sure dialog box goes away after selection
    full_folder_name = filedialog.askdirectory(parent=root, initialdir=default_location, 
                                               title=title_string, mustexist=True)    
    # return full folder name
    return full_folder_name    

def get_file(default_location, ext_list, title_string=None):
    """Dialog box to browse to a file.  Returns full file name.

    Usage: full_file_name = get_file(default_location, ext_list, [title]),
        where "default_location" is a starting folder location,
        ext_list is a list of (label, pattern) tuples,
        e.g. ext_list = [('Text files', '*.txt')],
        "title" is an optional message to list in the dialog box, and
        "full_file_name" is the complete name of the selected file.
    Written by Phil Wilmarth, OHSU, 2008, 2016.
    """
    # set up GUI elements
    root = tkinter.Tk()
    root.withdraw()
    try:
        root.tk.call('console', 'hide')
    except:
        pass
    
    # set default title string and ext list if not passed
    if title_string is None:   
        title_string = 'Select a single FILE'
    if not ext_list:
        ext_list =  [('All files', '*.*')]
    if not default_location:
        default_location = os.getcwd()
    
    # create dialog box for file selection
    root.update()   # helps make sure dialog box goes away after selection
    filename = filedialog.askopenfilename(parent=root, initialdir=default_location, 
                                          filetypes=ext_list, title=title_string)    
    # return full filename
    return filename      

def save_file(default_location, ext_list, default_file='', title_string=None):
    """Dialog box to save a file.  Returns full name of desired file.

    Usage: full_file_name = save_file(def_loc, ext_list, [def_file], [title]),
        where "def_loc" is a starting folder location,
        ext_list is a list of (label, pattern) tuples,
        e.g. ext_list = [('Text files', '*.txt')],
        "def_file" is an optional default filename,
        "title" is an optional message to list in dialog box, and
        "full_file_name" is the complete name of the desired file.
    Written by Phil Wilmarth, OHSU, 2009, 2016.
    """
    # set up GUI elements
    root = tkinter.Tk()
    root.withdraw()
    try:
        root.tk.call('console', 'hide')
    except:
        pass
    
    # set default title string if not passed
    if title_string is None:   
        title_string = 'Select a single FILE'
    if not ext_list:
        ext_list =  [('All files', '*.*')]
    if not default_location:
        default_location = os.getcwd()
        
    # create dialog box for file selection
    root.update()   # helps make sure dialog box goes away after selection
    filename = filedialog.asksaveasfilename(parent=root, initialdir=default_location, 
                                            initialfile=default_file, filetypes=ext_list, 
                                            title=title_string)
    # return full filename
    return filename    
    
    
def get_files(default_location, ext_list, title_string=None):
    """Dialog box to browse for files.  Returns a tuple of file names.

    Usage: file_name_list = get_files(default_location, ext_list, [title]),
        where "default_location" is a starting folder location,
        ext_list is a list of (label, pattern) tuples,
        e.g. ext_list = [('Text files', '*.txt')],
        "title" is an optional message to list in the dialog box, and
        "file_name_list" is a tuple of file name(s).
    Written by Phil Wilmarth, OHSU, 2010, 2016.
    """
    # set up GUI elements
    root = tkinter.Tk()
    root.withdraw()
    
    # set default title string if not passed
    if title_string is None:   
        title_string = 'Select one or more FILE(s)'
    if not ext_list:
        ext_list =  [('All files', '*.*')]
    if not default_location:
        default_location = os.getcwd()
        
    # create dialog box for file selection
    root.update()   # helps make sure dialog box goes away after selection
    filenames = filedialog.askopenfilenames(parent=root, initialdir=default_location, 
                                            filetypes=ext_list, multiple=True, 
                                            title=title_string)
    return filenames   

def get_string(title, prompt='Enter a string', initial=''):
    """Function to wrapper tkSimpleDialog.askstring function
    Written by Phil Wilmarth, OHSU, 2010.
    """
    from tkinter.simpledialog import askstring
    return askstring(title, prompt, initialvalue=initial)

    # end

################## some support functions/classes for PAW pipeline use #####################
# updated for 2017 Comet compatibility (new mod formats) -PW 10/27/2017

class Peptide:
    """An object for Comet peptide strings.
    """
    def __init__(self, sequence, delim='.', enzyme='Trypsin'):
        self.full_seq = sequence    # original string
        self.enzyme = enzyme
        self.prefix = None          # preceeding residue string
        self.seq = None             # actual peptide sequence
        self.suffix = None          # following residue string
        self.base_seq = None        # actual peptide sequence without any mods
        self.net = None             # number of enzymatic termini (given enzyme)
        self.length = None          # number of amino acids in base sequence

        # compile a couple of regex
        self.new_mods = re.compile(r'\[[-+]?([0-9]+(\.[0-9]*)?|\.[0-9]+)\]')
        self.old_mods = re.compile(r'[*#@^~$%!+nc\[\]\{\}\(\)]')

        # load attributes
        self.split_peptide(delim)
        self.compute_net(enzyme)

    def split_peptide(self, delim):
        """This splits SEQUEST/Comet peptide strings into prefix, sequence, and suffix.
        Computes some things and sets some attributes; supports the new bracketed
        floating point modification format (Comet 2017 and newer).
        """
        # this removes new Comet modification notation (bracketed floating points)
        base_seq = self.new_mods.sub('', self.full_seq)

        # probably have bounding residues delimited by periods
        items = base_seq.split(delim)
        if len(items) == 3:
            self.prefix, middle, self.suffix = items
            self.seq = self.full_seq[len(self.prefix) + 1: -(len(self.suffix) + 1)]
        elif len(items) == 1:
            self.prefix, self.suffix = 'X', 'X'
            middle = items[0]
            self.seq = self.full_seq
        else:
            print('WARNING: malformed peptide string:', self.full_seq)

        # remove older style modification symbols: *#@^~$%!+[](){} and 'n', 'c'
        self.base_seq = self.old_mods.sub('', middle)
        self.length = len(self.base_seq)
        return

    def _N_side_cleavage(self, prefix, prefix_pattern, nterm, nterm_pattern, suffix, suffix_pattern):
        """Computes number of termini constent with protease cleavage for N-terminal side cutters."""
        self.net = 0
        if (prefix in prefix_pattern) or (nterm in nterm_pattern):
            self.net += 1
        if suffix in suffix_pattern:
            self.net += 1

    def _C_side_cleavage(self, prefix, prefix_pattern, cterm, cterm_pattern, suffix, suffix_pattern, noP=True):
        """Computes number of termini constent with protease cleavage for C-terminal side cutters."""
        self.net = 0
        ct_okay = False
        if prefix in prefix_pattern:
            self.net += 1
        if (cterm in cterm_pattern) or (suffix in suffix_pattern):
            self.net += 1
            ct_okay = True
        if noP and (suffix == 'P') and (self.net > 0) and ct_okay:   # trypsin strict
            self.net -= 1
    
    def compute_net(self, enzyme):
        """Figures out the number of peptide termini consistent with the enzyme cleavage.
        Written by Phil Wilmarth, OHSU, 2008, rewritten 2017.
        """
        # valid amino acid characters
        amino_acids = set(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
                           'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'Y', 'Z', 'X'])
        
        # get the prefix amino acid residue
        i = len(self.prefix) - 1
        while self.prefix[i] not in amino_acids:
            i = i - 1
            if i < 0:
                break
        if i >= 0:
            prefix = self.prefix[i]
        else:
            prefix = 'X'
        # get suffix amino acid residue
        i = 0
        while self.suffix[i] not in amino_acids:
            i = i + 1
            if i >= len(self.suffix):
                break
        if i < len(self.suffix):
            suffix = self.suffix[i]
        else:
            suffix = 'X'
        
        cterm = self.base_seq[-1]  # last amino acid in sequence
        nterm = self.base_seq[0]   # first amino acid in sequence
        print(prefix, nterm, cterm, suffix)
        
        # determine number of enzymatic termini, net
        """need to support different enzymes and deal with proline.
        Seems Comet deals with premature stop codons as sequence breaks (* might be in prefix or suffix)."""
        if enzyme == 'Trypsin':         # cleaves at C-side of K, R (except if P)
            self._C_side_cleavage(prefix, 'KR-*', cterm, 'KR', suffix, '-*', noP=True)
        elif enzyme == 'Trypsin/P':     # cleaves at C-side of K, R
            self._C_side_cleavage(prefix, 'KR-*', cterm, 'KR', suffix, '-*', noP=False)
        elif enzyme == 'Lys_C':         # cleaves at C-side of K (except if P)
            self._C_side_cleavage(prefix, 'K-*', cterm, 'K', suffix, '-*', noP=True)
        elif enzyme == 'Lys_N':         # cleaves at N-side of K
            self._N_side_cleavage(prefix, '-*', nterm, 'K', suffix, 'K-*')
        elif enzyme == 'Arg_C':         # cleaves at C-side of R (except if P)
            self._C_side_cleavage(prefix, 'R-*', cterm, 'R', suffix, '-*', noP=True)
        elif enzyme == 'Asp_N':         # cleaves at N-side of D
            self._N_side_cleavage(prefix, '-*', nterm, 'D', suffix, 'D-*')
        elif enzyme == 'CNBr':          # cleaves at C-side of M
            self._C_side_cleavage(prefix, 'M-*', cterm, 'M', suffix, '-*', noP=False)
        elif enzyme == 'Glu_C':         # cleaves at C-side of D, E (except if P)
            self._C_side_cleavage(prefix, 'DE-*', cterm, 'DE', suffix, '-*', noP=True)
        elif enzyme == 'PepsinA':       # cleaves at C-side of F, L (except if P)
            self._C_side_cleavage(prefix, 'FL-*', cterm, 'FL', suffix, '-*', noP=True)
        elif enzyme == 'Chymotrypsin':  # cleaves at C-side of FWYL (except if P)
            self._C_side_cleavage(prefix, 'FWYL-*', cterm, 'FWYL', suffix, '-*', noP=True)
        elif enzyme == 'No_enzyme':
            self.net = 2
        else:
            print('WARNING: unknown enzyme specified', enzyme)
            self.net = 0

    def mask_base(self):
        """Masks I and L to j in base_seq."""
        return re.sub(r'[IL]', 'j', self.base_seq)
        
def find_peptide(peptide, proteins, mask=True, verbose=True):
    """Finds peptides in protein sequences.  Returns list of match tuples.
    This version requires additional attributes for protein objects.

    Usage: List = find_peptide(peptide, proteins, [limit=999]),
        where "peptide" is an uppercase peptide sequence string,
        "proteins" is a list of FASTAProtein objects,
        optional "limit" is the maximum number of protein matches (def=something large)
        and "List" is the returned list of match tuples.
            tuples: (accession, beginning res. #, ending res. #, full sequence)

    Written by Phil Wilmarth, OHSU, 9/4/08, 10/28/2015.
    """
    import types

    # make sure "proteins" is iterable
    protein_list = []
    if isinstance(proteins, list):
        protein_list = proteins
    elif isinstance(proteins, Protein):
        protein_list.append(proteins)
    else:
        print('FIND_PEPTIDE WARNING: "proteins" was not a list or Protein object.')
        
    matches = []
    for p in protein_list:
        matches += p.findPeptide(peptide, mask, pad_count=1)
    
    # if no matches, print warning (wrong database?)
    if verbose and len(matches) == 0:        
        print('FIND_PEPTIDE WARNING: "%s" not found in protein sequences.' % (peptide,))
        if len(protein_list) <= 20:
            for p in protein_list:
                print('...', p.accession)

    # return the match list (empty list if no matches)
    return matches 
    # end

class CometParams(object):
    """Container for Comet parameters values."""
    def __init__(self):
        """Creates the attribute placeholder/defaults."""
        
        # the most important attributes (initiate with defaults)
        """Do values need to be converted from strings to floats or ints?"""
        self.folder = None
        self.params_dict = {}                               # dictionary (parsed file)
        self.data_base = ''
        self.database_path = None
        self.peptide_mass_tolerance = 1.25
        self.peptide_mass_units = 'Da'
        self.parent_ion_mass_type = 'Mono'
        self.fragment_bin_tolerance = 1.0005
        self.fragment_bin_offset = 0.4
        self.fragment_ion_mass_type = 'Mono'
        self.ion_series = self.default_ion_series()         # dictionary
        self.search_enzyme_number = 1
        self.search_enzyme = 'Trypsin'
        self.num_enzyme_termini = 2
        self.missed_cleavages = 2
        self.variable_mods = self.default_variable_mods()   # dictionary
        self.max_variable_mods_in_peptide = 5
        self.static_mods = self.default_static_mods()       # dictionary
        self.enzyme_table = self.default_enzyme_table()     # dictionary
        
        self.output_lines = 20
        self.expect_score = True
        self.sample_enzyme_number = 1
        self.sample_enzyme = 'Tryp'
        
    def load_from_folder(self, folder):
        """Find any params files and parse comet first, sequest if no comet."""
        self.folder = folder    # location of the params file
        params_list = self.find_params_files(folder)
        self.which_one = None
        if 'comet.params' in [x.lower() for x in params_list]:
            self.which_one = 'comet'
            params_file = os.path.join(folder, 'comet.params')
            print('comet:', params_file)
            fin = open(params_file, 'rt')
            contents = fin.readlines()
            fin.close()
        elif 'sequest.params' in [x.lower() for x in params_list]:
            self.which_one = 'sequest'
            params_file = os.path.join(folder, 'sequest.params')
            print('sequest:', params_file)
            fin = open(params_file, 'rt')
            contents = fin.readlines()
            fin.close()
        else:
            print('   WARNING: no parameter files were found!')
            return
        
        # parse the params file and update the attributes
        if self.which_one == 'comet':
            self.parse_comet_params(contents) 
        else:
            self.parse_sequest_params(contents)
        self.load_attributes()
        return
    
    def load_from_string(self, string):
        """Find any params files and parse comet first, sequest if no comet."""
        # parse the params file and update the attributes
        self.which_one = 'comet'
        self.parse_comet_params(string.splitlines())
        self.load_attributes()
        return

    def load_attributes(self):
        """Populates the generic atributes from the respective parsed params dictionaries."""
        
        if self.which_one == 'comet':
            self.data_base = self.params_dict['database_name']
            self.database_path = os.path.dirname(self.data_base)
            self.peptide_mass_tolerance = float(self.params_dict['peptide_mass_tolerance'])
            self.peptide_mass_units = ['Da', 'milliDa', 'PPM'][int(self.params_dict['peptide_mass_units'])]
            self.parent_ion_mass_type = ['Ave', 'Mono'][int(self.params_dict['mass_type_parent'])]
            self.fragment_bin_tolerance = float(self.params_dict['fragment_bin_tol'])
            self.fragment_bin_offset = float(self.params_dict['fragment_bin_offset'])
            self.fragment_ion_mass_type = ['Ave', 'Mono'][int(self.params_dict['mass_type_fragment'])]
            self.search_enzyme_number = int(self.params_dict['search_enzyme_number'])
            self.search_enzyme = self.enzyme_table[self.search_enzyme_number][0]
            self.num_enzyme_termini = int(self.params_dict['num_enzyme_termini']) # need to set this for filtering protein matches
            if self.num_enzyme_termini in [1, 2, 8, 9]:
                if self.num_enzyme_termini in [8, 9]:
                    self.num_enzyme_termini = 1
            else:
                self.num_enzyme_termini = 0
            self.missed_cleavages = int(self.params_dict['allowed_missed_cleavage'])
            self.max_variable_mods_in_peptide = int(self.params_dict['max_variable_mods_in_peptide'])
            self.output_lines = int(self.params_dict['num_output_lines'])
            self.expect_score = [False, True][int(self.params_dict['print_expect_score'])]
            self.sample_enzyme_number = int(self.params_dict['sample_enzyme_number'])
            self.sample_enzyme = self.enzyme_table[self.sample_enzyme_number][0]
                                                          
        elif self.which_one == 'sequest':                                                          
            self.data_base = self.params_dict['first_database_name']
            self.database_path = os.path.dirname(self.data_base)
            self.peptide_mass_tolerance = float(self.params_dict['peptide_mass_tolerance'])
            self.peptide_mass_units = ['Da', 'milliDa', 'PPM'][int(self.params_dict['peptide_mass_units'])]
            self.parent_ion_mass_type = ['Ave', 'Mono'][int(self.params_dict['mass_type_parent'])]
            self.fragment_bin_tolerance = float(self.params_dict['fragment_ion_tolerance'])
            self.fragment_bin_offset = 0.0
            self.fragment_ion_mass_type = ['Ave', 'Mono'][int(self.params_dict['mass_type_fragment'])]
            if self.params_dict['enzyme_info'] == 'No_Enzyme 0 0 - -':
                self.search_enzyme_number = 0
            elif self.params_dict['enzyme_info'] == 'Trypsin(KR) 1 1 KR':
                self.search_enzyme_number = 1
            else:
                print('   WARNING: enzyme info was not no enzyme or trypsin - setting to trypsin')
                self.search_enzyme_number = 1
            self.search_enzyme = self.enzyme_table[self.search_enzyme_number][0]
            self.missed_cleavages = int(self.params_dict['max_num_internal_cleavage_sites'])
            self.max_variable_mods_in_peptide = int(self.params_dict['max_num_differential_per_peptide'])
            self.output_lines = int(self.params_dict['num_output_lines'])

        return

    def parse_comet_params(self, contents):
        """Parses comet.params file and loads params: values into a dictionary."""
        for i, line in enumerate(contents):
            if '[COMET_ENZYME_INFO]' in line:   # save the enzyme table parsing for later
                enzyme_start = i + 1
                break
            line = line.split('#')[0].strip()   # lines that start with "#" are comment lines
            if line:
                key = line.split('=')[0].strip()    # parameters and values separated by "="
                value = line.split('=')[1].strip()
                self.params_dict[key] = value

        # parse the enzyme table
        self.enzyme_table = {}
        for line in contents[enzyme_start:]:
            items = [x.strip() for x in line.split()]
            if len(items) == 5:
                key = int(items[0].replace('.', ''))
                value = tuple(items[1:])
                self.enzyme_table[key] = value

        # make ion series dictionary
        self.ion_series = {}
        keys = [k for k in self.params_dict.keys() if k.startswith('use_')]
        for k in keys:
            self.ion_series[k] = (True if self.params_dict[k] == '1' else False)
        
        # make variable mods dictionary
        self.variable_mods = {}
        keys = [k for k in self.params_dict.keys() if k.startswith('variable_')]
        for k in keys:
            self.variable_mods[k] = tuple([x.strip() for x in self.params_dict[k].split()])

        # make static mods dictionary
        self.static_mods = {}
        keys = [k for k in self.params_dict.keys() if k.startswith('add_')]
        for k in keys:
            self.static_mods[k] = self.params_dict[k]
                                       
    def parse_sequest_params(self, contents):
        """Parses SEQUEST.PARAMS files and loads params: values into a dictionary"""
        for i, line in enumerate(contents):
            line = line.split('#')[0].strip()   # lines that start with "#" are comment lines
            if line:
                line = line.split(';')[0].strip()
                if line and '=' in line:
                    key = line.split('=')[0].strip()
                    value = line.split('=')[1].strip()
                    self.params_dict[key] = value
    
        # make Comet-style ion series dictionary
        """
        1--> 0 or 1 whether neutral losses of series A should be honoured. (1 = yes)
        2--> 0 or 1 whether neutral losses of series B should be honoured. (1 = yes)
        3--> 0 or 1 whether neutral losses of series Y should be honoured. (1 = yes)
        4 -> factor for series A
        5--> factor for series B
        6--> factor for series C
        7--> factor for series D
        8--> factor for series V
        9--> factor for series W
        10-> factor for series X
        11-> factor for series Y
        12-> factor for series Z
        """
        values = [float(x) for x in self.params_dict['ion_series'].split()]

        # see if any neutal loss ions were turned on (value=1)
        if len([x for x in values[:3] if x > 0]) > 0:
            self.ion_series['use_NL_ions'] = True
        else:
            self.ion_series['use_NL_ions'] = False

        # set the main ion series flags
        self.ion_series['use_A_ions'] = (True if values[3] > 0.0 else False)
        self.ion_series['use_B_ions'] = (True if values[4] > 0.0 else False)
        self.ion_series['use_C_ions'] = (True if values[5] > 0.0 else False)
        self.ion_series['use_X_ions'] = (True if values[9] > 0.0 else False)
        self.ion_series['use_Y_ions'] = (True if values[10] > 0.0 else False)
        self.ion_series['use_Z_ions'] = (True if values[11] > 0.0 else False)
        
        # make variable mods dictionary
        mods = self.params_dict['diff_search_options'].split()
        self.mod_tuples = list(zip(mods[::2], mods[1::2]))
        
        # term_diff_search_options parameter is C-term deltamass, then N-term deltamass
        self.mod_tuples += list(zip(self.params_dict['term_diff_search_options'].split(), ['c', 'n']))
        self.mod_tuples = [tuple(list(x) + ['0', '3', '-1', '0', '0']) for x in self.mod_tuples]

        """Need to remove any deltamass of zero entries"""

        # load into dictionary (Comet style mods)
        for i, tuple_ in enumerate(self.mod_tuples):
            self.variable_mods['variable_mod0%d' % (i+1,)] = tuple_

        # make static mods dictionary
        self.static_mods = {}
        keys = [k for k in self.params_dict.keys() if k.startswith('add_')]
        for k in keys:
            self.static_mods[k] = self.params_dict[k]
        return

    def find_params_files(self, folder):
        return [p for p in os.listdir(folder) if p.endswith('.params')]

    def default_ion_series(self):
        """Sets default ion series (B, Y and neutral loss)."""
        return OrderedDict([('use_A_ions', False),
                            ('use_B_ions', False),
                            ('use_C_ions', False),
                            ('use_X_ions', False),
                            ('use_Y_ions', False),
                            ('use_Z_ions', False),
                            ('use_NL_ions', False)])
    
    def default_variable_mods(self):
        """Sets default variable mods."""
        return OrderedDict([('variable_mod01', ('15.9949', 'M', '0', '3', '-1', '0', '0')),
                            ('variable_mod02', ('0.0', 'X', '0', '3', '-1', '0', '0')),
                            ('variable_mod03', ('0.0', 'X', '0', '3', '-1', '0', '0')),
                            ('variable_mod04', ('0.0', 'X', '0', '3', '-1', '0', '0')),
                            ('variable_mod05', ('0.0', 'X', '0', '3', '-1', '0', '0')),
                            ('variable_mod06', ('0.0', 'X', '0', '3', '-1', '0', '0')),
                            ('variable_mod07', ('0.0', 'X', '0', '3', '-1', '0', '0')),
                            ('variable_mod08', ('0.0', 'X', '0', '3', '-1', '0', '0')),
                            ('variable_mod09', ('0.0', 'X', '0', '3', '-1', '0', '0'))])

    def default_static_mods(self):
        """Sets default static modifications."""
        return OrderedDict([('add_Cterm_peptide', '0.0000'),
                            ('add_Nterm_peptide', '0.0000'),
                            ('add_Cterm_protein', '0.0000'),
                            ('add_Nterm_protein', '0.0000'),
                            ('add_G_glycine', '0.0000'),
                            ('add_A_alanine', '0.0000'),
                            ('add_S_serine', '0.0000'),
                            ('add_P_proline', '0.0000'),
                            ('add_V_valine', '0.0000'),
                            ('add_T_threonine', '0.0000'),
                            ('add_C_cysteine', '57.021464'),
                            ('add_L_leucine', '0.0000'),
                            ('add_I_isoleucine', '0.0000'),
                            ('add_N_asparagine', '0.0000'),
                            ('add_D_aspartic_acid', '0.0000'),
                            ('add_Q_glutamine', '0.0000'),
                            ('add_K_lysine', '0.0000'),
                            ('add_E_glutamic_acid', '0.0000'),
                            ('add_M_methionine', '0.0000'),
                            ('add_O_ornithine', '0.0000'),
                            ('add_H_histidine', '0.0000'),
                            ('add_F_phenylalanine', '0.0000'),
                            ('add_U_selenocysteine', '0.0000'),
                            ('add_R_arginine', '0.0000'),
                            ('add_Y_tyrosine', '0.0000'),
                            ('add_W_tryptophan', '0.0000'),
                            ('add_B_user_amino_acid', '0.0000'),
                            ('add_J_user_amino_acid', '0.0000'),
                            ('add_X_user_amino_acid', '0.0000'),
                            ('add_Z_user_amino_acid', '0.0000')])

    def default_enzyme_table(self):
        return OrderedDict([(0, ('No_enzyme', '0', '-', '-')),
                            (1, ('Trypsin', '1', 'KR', 'P')),
                            (2, ('Trypsin/P', '1', 'KR', '-')),
                            (3, ('Lys_C', '1', 'K', 'P')),
                            (4, ('Lys_N', '0', 'K', '-')),
                            (5, ('Arg_C', '1', 'R', 'P')),
                            (6, ('Asp_N', '0', 'D', '-')),
                            (7, ('CNBr', '1', 'M', '-')),
                            (8, ('Glu_C', '1', 'DE', 'P')),
                            (9, ('PepsinA', '1', 'FL', 'P')),
                            (10, ('Chymotrypsin', '1', 'FWYL', 'P'))])

    def _snoop(self):
        """Diagnostic console dump of attributes."""
        print('data_base:', self.data_base)
        print('database_path:', self.database_path)
        print('peptide_mass_tolerance:', self.peptide_mass_tolerance)
        print('peptide_mass_units:', self.peptide_mass_units)
        print('parent_ion_mass_type:', self.parent_ion_mass_type)
        print('fragment_bin_tolerance:', self.fragment_bin_tolerance)
        print('fragment_bin_offset:', self.fragment_bin_offset)
        print('fragment_ion_mass_type:', self.fragment_ion_mass_type)
        print('search_enzyme_number:', self.search_enzyme_number)
        print('search_enzyme:', self.search_enzyme)
        print('missed_cleavages', self.missed_cleavages)
        print('max_variable_mods_in_peptide:', self.max_variable_mods_in_peptide)

        print('output_lines:', self.output_lines)
        print('expect_score:', self.expect_score)
        print('sample_enzyme_number:', self.sample_enzyme_number)
        print('sample_enzyme:', self.sample_enzyme)
        
        print('ion_series:')    # dictionary
        pprint(self.ion_series)
        print('variable_mods:')   # dictionary
        pprint(self.variable_mods)
        print('static_mods:')       # dictionary
        pprint(self.static_mods)       # dictionary        
        return
    
    # end class

class PeptideInfo:
    """Data structure for some basic peptide information."""
    def __init__(self, sequence='', begin=0, end=0, mass=0, missed=0):
        self.seq = sequence
        self.beg = begin
        self.end = end
        self.mass = mass
        self.missed = missed
        return            

class Protein:
    """Object to hold protein accession numbers, descriptions, and sequences.

    Methods:
        __init_:standard constructor, no parameters.
        readProtein: returns next protein from "fasta_reader"
        printProtein: prints sequence in FASTA format
        reverseProtein: reverses sequences and modifies accession/descriptions
        molwtProtein: computes average MW of sequence
        frequencyProtein: returns aa composition dictionary
        seqlenProtein: returns aa sequence length
        findPeptide: finds location of peptide in protein sequence
        coverage: calculates coverage and aa counts from peptide list
        enzymaticDigest: theroetical enzymatic digest of protein sequence
        
    Written by Phil Wilmarth, OHSU, 2009, 2016.
    
    Updated for new Comet mod formats -PW 10/27/2017
    Removed any parsing of accessions and descriptions methods and attributes -PW 20180711
    """
    def __init__(self):
        """Basic constructor, no parameters.
        """
        # bare bones __init__ function
        self.accession = 'blank'
        self.description = 'blank'
        self.sequence = ''
        self.sequence_padded = None
        self.sequence_masked = None
        self.pad_count = None
        self.length = 0
        self.peptides = []
        return
    
    def readProtein(self, fasta_reader):
        """Gets the next FASTA protein entry from FastaReader object.

        Usage: Boolean = object.readProtein(fasta_reader),
            where "object" is an instance of a Protein object and
            "fasta_reader" is an instance of a FastaReader object.
            Return value is "False" when EOF encountered.
        """
        status = fasta_reader.readNextProtein(self)
        return status
    
    def printProtein(self, file_obj=None, length=80):
        """Prints FASTA protein entry to file (stdout is default).

        Usage: object.printProtein([file_obj=None, length=80]),
            where "object" is an instance of a Protein object, and
            "file_obj" is a file object (a value of None will print
            to standard out stream.  Optional "length" is number of
            characters per line for the protein sequence.
        """
        if file_obj == None:
            file_obj = sys.stdout
        
        # print new accession and new descriptor on first line
        if self.description == '':
            print('>'+self.accession, file=file_obj)
        else:
            print('>'+self.accession, self.description, file=file_obj)
        
        # initialize some things
        char_count = 0
        char_line = ''
        
        # build up sequence line with "length" characters per line
        for char in self.sequence:
            if char_count < length:  # do not have "width" chars yet
                char_line += char
                char_count += 1
            else:                   # line is "width" long so print and reset
                print(char_line, file=file_obj)
                char_line = char
                char_count = 1
        
        # print last sequence line (often less than "width" long) and return
        if len(char_line):
            print(char_line, file=file_obj)
        return

    def reverseProtein(self, decoy_string):
        """Reverses protein sequence and returns new Protein object.

        Usage: rev_prot = object.reverseProtein(decoy_string),
            where "object" is a Protein object, "decoy_string" is the
            unique identifier text to add to the beginning of the 
            protein accesion number, and "rev_prot" is new Protein object.
        """
        # make sure decoy_string ends with an undescore
        if not decoy_string.endswith('_'):
            decoy_string = decoy_string + '_'
        
        # create a new Protein instance
        rev_prot = Protein() 
        
        # prefix the decoy_string to desired parts of accession
        if self.accession.startswith('CONT_'):
            new_acc = decoy_string + self.accession.split('|')[0]
        else:
            new_acc = decoy_string + self.accession.replace('|', '&') # best to remove "|"        
        rev_prot.accession = new_acc
        
        # change the desciptions, too.
        rev_prot.description = 'REVERSED'
        
        # reversed the protein sequence and return new protein object
        rev_prot.sequence = self.sequence[::-1]
        return rev_prot

    def molwtProtein(self, show_errs=True):
        """Returns protein molecular weight as the sum of average aa masses.
        If "show_errs" flag set, invalid amino acid characters are reported.
        Does not add any modification delta masses (fixed or variable).
        """        
        # start with water then add aa masses
        self.setMasses()
        bad_char = {}
        molwt = self.ave_masses['water']
        for aa in self.sequence:
            try:
                molwt += self.ave_masses[aa]
            except:     # keep track of bad characters
                bad_char[aa] = True
        
        bad_char = sorted(bad_char.keys())
        if len(bad_char) > 0 and show_errs:     # report bad chars if desired
            print('   WARNING: unknown symbol(s) (%s) in %s:\n%s' %
                  (''.join(bad_char), self.accession, self.sequence))
        return molwt

    def frequencyProtein(self, show_errs=True):
        """Returns aa frequency distrubution as a dictionary.
        If "show_errs" flag set, invalid amino acid characters are reported.
        """
        freq = {'X':0, 'G':0, 'A':0, 'S':0, 'P':0, 'V':0, 'T':0,
                'C':0, 'L':0, 'I':0, 'J':0, 'N':0, 'O':0, 'B':0,
                'D':0, 'Q':0, 'K':0, 'Z':0, 'E':0, 'M':0, 'H':0,
                'F':0, 'R':0, 'Y':0, 'W':0, 'U':0, '*':0, '-':0 }
        
        # count the amino acids for all residues in sequence
        bad_char = {}
        for aa in self.sequence:
            try:
                freq[aa] += 1
            except:     # keep track of bad characters
                bad_char[aa] = True
        
        bad_char = sorted(bad_char.keys())
        if len(bad_char) > 0 and show_errs: # report any bad chars, if desired
            print('   WARNING: unknown symbol(s) (%s) in %s:\n%s' %
                  (''.join(bad_char), self.accession, self.sequence))
        return freq
    
    def seqlenProtein(self):
        """Calculates protein sequence length.
        """        
        self.length = len(self.sequence)
        return self.length

    def split_peptide(self, sequence):
        """Splits peptide assuming that there might be single preceeding and following residues with periods."""
        if re.match(r'[-A-Z]\..+\.[-A-Z]', sequence):
            return sequence[0], sequence[2:-2], sequence[-1]
        else:
            if min(sequence.count('['), sequence.count(']')) != sequence.count('.'):
                print('   WARNING: possible malformed peptide string:', sequence)
            return '', sequence, ''

    def peptide_decorations(self, sequence):
        """Separate modifications from amino acid residues so that mods can be put back later."""
        residues = []
        decorations = []
        char_count = 0
        decoration = ''
        for char in sequence:
            if re.match('[A-Z]', char):
                residues.append(char)
                decorations.append(decoration)
                decoration = ''
            else:
                decoration += char
        # might have C-terminal mod
        residues.append('')
        decorations.append(decoration)
        
        return residues, decorations

    def redecorate_peptide(self, peptide, decorations):
        """Redecorates a peptide sequence with mods."""
        residues = list(peptide + '')
        return ''.join(['' + x + y for (x, y) in zip(decorations, residues)])

    def base_peptide_sequence(self, sequence, mask=True):
        """Returns the peptide amino acid residues from SEQUEST peptide strings
        """
        # remove bounding residues (SEQUEST/Comet format: A.BCD.E)
        prefix, peptide, suffix = self.split_peptide(sequence)

        # remove the 2017 Comet style mod strings
        peptide = re.sub(r'\[[-+]?[0-9]*(.)?[0-9]*\]', '', peptide)        
        # remove modification symbols: '*', '#', '@', '^', '~', '$', '%', '!', '+', 'n', 'c', '[', ']', "(', ')', '{', '}'
        peptide = re.sub(r'[*#@^~$%!+nc\[\]\{\}\(\)]', '', peptide)
        
        # mask I/L if needed:
        if mask:
            return re.sub(r'[IL]', 'j', peptide)
        else:
            return peptide

    def findPeptide(self, peptide, mask=True, pad_count=1):
        """Calculates location of all 'peptide' matches in 'self.sequence.'
        Returns a match tuple list.
        Match tuples: (accession, beginning res. #, ending res. #, peptide sequence with context)
        """
        matches = []

        # remove bounding residues (SEQUEST/Comet format: A.BCD.E)
        prefix, middle, suffix = self.split_peptide(peptide)

        # get a clean peptide sequence to lookup (retain mods)
        residues, decorations = self.peptide_decorations(middle)
        base_pep_masked = ''.join(residues)
        if mask:
            base_pep_masked = re.sub(r'[IL]', 'j', base_pep_masked)
        
        # fix the protein sequence for peptide lookups (pad and mask I/L). Save the results to save time
        if (not self.sequence_masked) or (pad_count != self.pad_count):
            self.sequence_padded = ('-' * pad_count) + self.sequence + ('-' * pad_count) # add bounding symbols
            if mask:
                self.sequence_masked = re.sub(r'[IL]', 'j', self.sequence_padded)
            else:
                self.sequence_masked = self.sequence_padded
            self.pad_count = pad_count
        
        # find all matches of base_pep_masked to protein sequence (padded and masked)
        search_string = '(?=(%s))' % base_pep_masked    # need to do this to get overlapping matches (new python regex not yet ready: overlapped=True flag)
        for match in re.finditer(search_string, self.sequence_masked):
            
            start = match.span()[0] # NOTE: look ahead matching does not give ending position of string (beg=end)
            end = start + len(base_pep_masked)
            start_prot, end_prot = start - self.pad_count + 1, end - self.pad_count

            # add bounding AAs, periods, and put back modification special chars
            pre = self.sequence_padded[start-self.pad_count:start]
            post = self.sequence_padded[end:end+self.pad_count]
            middle = self.redecorate_peptide(self.sequence_padded[start:end], decorations)
            full_seq = pre + '.' + middle + '.' + post
            """might want to create a match object instead of tuple."""
            matches.append((self.accession, start_prot, end_prot, full_seq))
        
        # return the match list (empty list if no matches)
        return matches
        
    def calcCoverage(self, peptide_list):
        """Calculates % coverage and aa frequency map of matched peptides.
        "peptide_list" is list of sequences with optional counts (as tuples).
        """
        freq_dict = {}
        try:    # see if peptide_list is a list of tuples or not
            for peptide, count in peptide_list:
                for (acc, beg, end, seq) in self.findPeptide(peptide):
                    for key in [str(i) for i in range(beg, end+1)]:
                        if freq_dict.get(key, False):
                            freq_dict[key] = freq_dict[key] + count
                        else:
                            freq_dict[key] = count
        except ValueError:
            for peptide in peptide_list:
                for (acc, beg, end, seq) in self.findPeptide(peptide):
                    for key in [str(i) for i in range(beg, end+1)]:
                        if freq_dict.get(key, False):
                            freq_dict[key] = freq_dict[key] + 1
                        else:
                            freq_dict[key] = 1
                            
        coverage = 100.0*float(len(freq_dict))/float(len(self.sequence))
        coverage_map = []
        for i, aa in enumerate(self.sequence):
            coverage_map.append((str(i+1), aa, freq_dict.get(str(i+1), 0)))
        return (coverage, coverage_map)
    
    def enzymaticDigest(self, enzyme_regex=None, low=400.0, high=6000.0, length=6, missed=3, mass='mono'):
        """Performs a tryptic digest of a protein sequence. This does not
        do any modifications to residues except for reduction/alkylation of
        cys residues (C+57). Mass filters should be relaxed.
        
        Returns a list of digested peptides.
        enzyme_regex is a compiled re object for the enzyme cleavage
            (if enzyme_regex not defined, do tryptic digest by default)
        low, high - mass limits for peptides.
        length - minimum amino acid length
        missed - maximum number of missed cleavages.
        mass - 'ave' average or 'mono' monoisotopic masses.
        """

        """Regular expression digestion table:
        trypsin from: http://stackoverflow.com/questions/18364380/python-3-cut-peptide-regular-expression

        regex = re.compile(r".")                        # no enzyme
        regex = re.compile(r".(?:(?<![KR](?!P)).)*")    # trypsin strict
        regex = re.compile(r".(?:(?<![KR]).)*")         # trypsin with cleavage at P
        regex = re.compile(r".(?:(?<![K](?!P)).)*")     # Lys-C strict
        regex = re.compile(r".(?:(?<![K]).)*")          # Lys-C with cleavage at P
        regex = re.compile(r".(?:(?![K]).)*")           # Lys-N
        regex = re.compile(r".(?:(?<![R](?!P)).)*")     # Arg-C strict
        regex = re.compile(r".(?:(?![D]).)*")           # Asp-N
        regex = re.compile(r".(?:(?<![M]).)*")          # CnBr
        regex = re.compile(r".(?:(?<![DE](?!P)).)*")    # Glu-C
        regex = re.compile(r".(?:(?<![FL](?!P)).)*")    # PepsinA
        regex = re.compile(r".(?:(?<![FWYL](?!P)).)*")  # chymotrypsin
        """   
        # skip if there is no sequence to digest
        if len(self.sequence) == 0:
            return []

        # tryptic digestion is the default
        if not enzyme_regex:
            enzyme_regex = re.compile(r".(?:(?<![KR](?!P)).)*")
        
        # set up masses, default is alkylated cysteine. No mechanism for other modifications yet.
        self.setMasses()
        if mass == 'ave':
            masses = copy.deepcopy(self.ave_masses)
            masses['C'] = 160.197
        elif mass == 'mono':
            masses = copy.deepcopy(self.mono_masses)
            masses['C'] = 160.03065
        else:
            print('...WARNING: masses must be "ave" or "mono"')

        # digest the sequence
        digest_matches = [x for x in enzyme_regex.finditer(self.sequence)] # list of re match objects

        # get info from match objects into PeptideInfo object attributes
        digest = [PeptideInfo(mass=masses['water']) for x in digest_matches]
        for i, match in enumerate(digest_matches):
            digest[i].seq = match.group()
            digest[i].beg, digest[i].end = match.span()
            digest[i].beg += 1
            for aa in match.group():
                digest[i].mass += masses[aa]
            
        # test peptides and missed cleavage peptides for mass ranges and min length
        valid_digest = []
        for i in range(len(digest)):
            
            # check if peptide is within the mass range and meets min length
            if (low <= digest[i].mass <= high) and (len(digest[i].seq) >= length):
                valid_digest.append(digest[i])
                
            # create and check missed cleavages
            for j in range(1, missed+1):
                if (i+j) > len(digest)-1:
                    continue
                temp = PeptideInfo(begin=100000)    # a peptide object for missed cleavages
    
                # calculate running sums for each number of missed cleavages
                for k in range(j+1):
                    if (i+k) > len(digest)-1:
                        continue
                    temp.seq += digest[i+k].seq
                    temp.beg = min(temp.beg, digest[i+k].beg)
                    temp.end = max(temp.end, digest[i+k].end)
                    temp.mass += (digest[i+k].mass - masses['water'])
                temp.mass += masses['water']
                temp.missed = k
                
                # check missed cleavage peptide for valid mass range and length
                if (low <= temp.mass <= high) and (len(temp.seq) >= length):
                    valid_digest.append(temp)
        
        # return the list of digested peptides
##        self.peptides = valid_digest  # this saves the digest results for each protein (uses memeory)
        return valid_digest

    def setMasses(self):
        """Set average and monoisotopic mass dictionaries."""
        self.ave_masses = {'X':  0.0000, 'G': 57.0513, 'A': 71.0779, 'S': 87.0773, 'P': 97.1152,
                           'V': 99.1311, 'T':101.1039, 'C':103.1429, 'L':113.1576, 'I':113.1576,
                           'J':113.1576, 'N':114.1026, 'O':114.1472, 'B':114.5950, 'D':115.0874,
                           'Q':128.1292, 'K':128.1723, 'Z':128.6216, 'E':129.1140, 'M':131.1961,
                           'H':137.1393, 'F':147.1739, 'R':156.1857, 'Y':163.1733, 'W':186.2099,
                           'U':150.0379, '*': 0.00000, '-': 0.00000, 'water':18.02}
        self.mono_masses = {'X':  0.000000, 'G': 57.021464, 'A': 71.037114, 'S': 87.032028, 'P':97.052764,
                            'V': 99.068414, 'T':101.047679, 'C':103.009185, 'L':113.084064, 'I':113.084064,
                            'J':113.084064, 'N':114.042927, 'O':114.147200, 'B':114.595000, 'D':115.026943,
                            'Q':128.058578, 'K':128.094963, 'Z':128.621600, 'E':129.042593, 'M':131.040485,
                            'H':137.058912, 'F':147.068414, 'R':156.101111, 'Y':163.063320, 'W':186.079313,
                            'U':150.953630, '*':  0.000000, '-':  0.000000, 'water':18.01057}
        return                                             
    # end class

class FastaReader:
    """Reads FASTA entries from a file-like object.

    methods:
    __init__: basic constructor, no parameters.
    
    readProtein: reads one FASTA entry from a file object (text or zipped)
        arguments are "next_protein" and "file_obj"
        returns True (next protein) or False (EOF or not FASTA).

    written by Phil Wilmarth, OHSU, 2009.
    """
    def __init__(self, fasta_file):
        """Basic constructor function.  No parameters

        self._last_line retains the previous '>' line and
        self._valid is a dictionary of valid protein FASTA chars.
        """        
        # attribute to save last line from previous read
        self._last_line = 'start value'
        self._file_obj = None
        self._fasta_file = fasta_file
        
        # list of valid amino acid characters
        self._valid = {'X':True, 'G':True, 'A':True, 'S':True, 'P':True,\
                      'V':True, 'T':True, 'C':True, 'L':True, 'I':True,\
                      'J':True, 'N':True, 'O':True, 'B':True, 'D':True,\
                      'Q':True, 'K':True, 'Z':True, 'E':True, 'M':True,\
                      'H':True, 'F':True, 'R':True, 'Y':True, 'W':True,\
                      'U':True, '*':True, '-':True }
        
        # get file object and save as attribute
        if not os.path.exists(fasta_file):
            ext_list = [('FASTA files', '*.fasta'), ('Zipped FASTA files', '*.gz'), ('All files', '*.*')]
            fasta_file = get_file(default_location, extension_list, title_string="Select FASTA file")
        try:
            if fasta_file.endswith('.gz'):
                self._file_obj = gzip.open(fasta_file, 'rt')
            else :
                self._file_obj = open(fasta_file, 'rt')
        except IOError:
            print('   WARNING: Fasta database could not be opened!')
            raise
        return

    def readNextProtein(self, next_protein, check_for_errs=False):
        """Loads one FASTA protein text entry into a Protein object.

        Returns True (protein entry found) or False (end of file).
        If "check_for_errs" flag is set, amino acid chars are checked.

        Written by Phil Wilmarth, OHSU, 2009.
        """
        # at first call, start reading lines
        if self._last_line == 'start value':
            self._last_line = self._file_obj.readline()
            if not self._last_line:
                self._file_obj.close()
                return(False)
            self._last_line = self._last_line.strip()
        
        # get next protein's info from _last_line
        if self._last_line.startswith('>'):
            next_protein.accession = self._last_line.split()[0][1:]
            start = len(next_protein.accession)+2
            next_protein.description = self._last_line[start:]
        
        # return if empty line (EOF) or non-description line
        else:
            self._file_obj.close()
            return(False)                    
        
        # reset variables and read in next entry
        next_protein.sequence = ""
        line = self._last_line
        self._last_line = ""
        bad_char = {}
        while line:
            line = self._file_obj.readline()
            if not line:
                break
            else:
                testline = line.strip()
            if testline == '':
                continue
            
            # stop reading at next descriptor line (and save line)
            if line.startswith('>'):
                self._last_line = line.strip()
                
                # report bad characters if conditions were met
                bad_char = sorted(bad_char.keys())
                if len(bad_char) > 0 and check_for_errs:
                    print('   WARNING: unknown symbol(s) (%s) in %s' %
                          (''.join(bad_char), next_protein.accession))
                break
            
            # add next sequence line to protein's sequence
            else:
                line = line.rstrip()
                line = line.upper()
                if check_for_errs: # checking chars slows down the program
                    for char in line:
                        if self._valid.get(char, False):
                            next_protein.sequence += char
                        else:
                            bad_char[char] = True                
                else: # blindly adding the line is faster...
                    next_protein.sequence += line
        
        # return (protein info retained in next_protein)
        return True

    # end class

class PAWShell(object):
    """Command line loop to define sample names. Oct. 2017 -PW"""
    
    def __init__(self, folder):
        """Set attributes."""
        self.folder = folder
        self.prompt = '\n(Auto Each Help List Pattern Reset Show Done) ? '
        self.samples = {}
        self.which_sample = {}
        self._get_filtered_files()
        self.remaining = list(self.files) # so we get a copy
        self.quit = False   # 

    def cmdloop(self, message):
        """Main command line loop."""
        print(message)
        while True:
            # get user response and parse
            cmd_line = input(self.prompt)
            cmd = cmd_line.split()[0].upper()
            arg = cmd_line[len(cmd):].lstrip()

            # repsond to command
            if not cmd:
                continue
            elif cmd[0] == 'A':
                self.do_auto(arg)
            elif cmd[0] == 'E':
                self.do_each(arg)
            elif cmd[0] == 'H':
                self.do_help(arg)
            elif cmd[0] == 'L':
                self.do_list(arg)
            elif cmd[0] == 'P':
                self.do_pattern(arg)
            elif cmd[0] == 'R':
                self.do_reset(arg)
            elif cmd[0] == 'S':
                self.do_show(arg)
            elif cmd[0] == 'D':
                self.do_done(arg)

            # see if done
            if self.quit:
                break

        # return the sample mappings                    
        return self.samples, self.which_sample

    def do_auto(self, arg=None):
        """Automatically parse the sample names.

        Uses multi-replacement from: https://stackoverflow.com/questions/6116978/python-replace-multiple-strings.
        """
        for txt_name in self.files:            
            # see if project code (e.g. ABC-999) is at beginning
            project = None
            m = re.match(r'[A-Z]{3,4}[-]?[0-9]{2,4}', txt_name)
            if m:
                project = m.group()
                to_mask = txt_name[len(project):]
                
            # remove some common labels so fraction numbers get processed correctly (use regex here - see above)   
            remove = ['mm_', 'KCL_', 'inject_', 'FT_', '_LT', '_VE2', '_VE', '_OT', '_QE', '_SCXFrac', '_SCX',
                      'pctACN_', '%ACN_', 'cm_', 'min_', 'ug_', '_fr', '_frac', '_fx', '_TMT', '_Fxn']
            rep = {k: '_' for k in remove}

            # this does the replacements (see above citation)
            rep = dict((re.escape(k), v) for k, v in rep.items())
            pattern = re.compile("|".join(rep.keys()), flags=re.IGNORECASE)
            sample_name = pattern.sub(lambda m: rep[re.escape(m.group(0))], txt_name)

            # might have gel bands
            sample_name = re.sub(r'1Dband|band|_filtered.txt.gz|_filtered.txt', '', sample_name, flags=re.IGNORECASE)
            
            # remove small numberic strings (less than 4 digits) and replace common separators with underscore
            done = False
            while not done:
                old_name = sample_name
                sample_name = re.sub(r'[0-9]pt[0-9]|[._-][0-9]{1,3}[._-]?|[_]+', '_', sample_name, flags=re.IGNORECASE)
                if old_name == sample_name:
                    done = True
            if sample_name.endswith('_'):
                sample_name = sample_name[:-1]            
            if sample_name.startswith('_'):
                sample_name = sample_name[1:]

            # add back the project code (if any)
            if project:
                sample_name = project + '_' + sample_name
            
            # update samples dictionary (sample: list of filenames)            
            self.which_sample[txt_name] = sample_name
            if sample_name in self.samples:
                self.samples[sample_name].append(txt_name)
            else:
                self.samples[sample_name] = [txt_name]

        # show results
        self.do_show()
        self.remaining = []
        return

    def do_each(self, arg=None):
        """Each file will be a separate sample."""
        sample_names = [x.split('_filtered.txt')[0] for x in self.files]
        self.samples = {k: [v] for (k, v) in zip(sample_names, self.files)}
        self.which_sample = {k: v for (k, v) in zip(self.files, sample_names)}
        self.do_show()
        self.remaining = []
        return

    def do_help(self, arg=None):
        """Prints help."""
        print("""\nCommands: Auto, Each, Help, List, Pattern, Reset, Show, Quit.
(commands are not case sensitive and first letter is sufficient.)

...Auto: automatically parse the sample names from the filenames,
...Each: treat each filename as a separate sample,
...Help: prints this message,
...List: lists all of the (remaining) filename(s),
...Pattern pattern: glob-style search pattern to get a subset of filenames
      ("*" or "*.*" is all files,
       "*.txt" is all files that have a "txt" extension,
       "a*" is all files that start with the letter "a",
       "*string*" is all files that contain "string", etc.),
...Reset: start over,
...Show: print the current samples and associated files [verbose],
...Done: exit from the command loop and return to processing.""")
        return
       
    def do_list(self, arg=None):
        """List the files."""
        self._list()
        return

    def do_pattern(self, arg=None):
        """Find files matching glob pattern and assign to sample name."""
##        selected = fnmatch.filter(self.remaining, arg)    # case sensitive
        selected = fnmatch.filter([x.upper() for x in self.remaining], arg.upper())
        selected = [x for x in self.remaining if x.upper() in selected]
        if not selected:
            print('  WARNING: pattern did not match any files')
            return
        self._list(selected, message='  Selected files:')

        # get the name for the sample set
        sample = ''
        while not sample:
            sample = input('PAW> Sample name? ')
        if sample not in self.samples:
            self.samples[sample] = selected
        else:
            self.samples[sample] = sorted(self.samples[sample] + selected)
        for f in selected:
            self.which_sample[f] = sample

        # remove the selected files from the remaining list
        [self.remaining.remove(x) for x in selected]
        return

    def do_reset(self, arg=None):
        """Resets dictionaries and file list."""
        self.remaining = list(self.files)   # make a copy
        self.samples = {}
        self.which_sample = {}
        return

    def do_show(self, arg=None):
        """Show the current status of sample definitions."""
        print('  %s samples are defined:' % len(self.samples))
        for s in self.samples:
            print('  ..%s has %d files:' % (s, len(self.samples[s])))
            for i, f in enumerate(self.samples[s]):
                print('  ....(%2d) %s' % (i+1, f))

        if arg and arg.lower() == 'verbose':
            print('\n  LC_run/TXT_file to sample mappings:')
            for i, (f, s) in enumerate(self.which_sample.items()):
                print('  ..(%2d) %s: %s' % (i+1, f, s))
        return
    
    def do_done(self, arg=None):
        """Quit the loop and return the dictionaries."""
        self.quit = True
        return
    
    def _get_filtered_files(self):
        """Gets list of filtered TXT files."""
        save = os.getcwd()
        os.chdir(self.folder)
        self.files = glob.glob('*_filtered.txt*')
        os.chdir(save)
        return

    def _list(self, files=None, message='  Remaining files:'):
        """List file names."""
        if not files:
            files = self.remaining
        if files == self.files:
            message = '  All files:'
        print(message)
        for i, f in enumerate(files):
            print('    (%2d) %s' % (i+1, f))
        return
    # end class

def amino_acid_count(sequence_string, enzyme='Tryp', return_base_pep=False):
    """Counts amino acids in peptides.  Returns (length, ntt) tuple.

    Usage: (length, ntt) = amino_acid_count(sequence_string),
        where "sequence_string" is a peptide sequence with bounding residues,
        "enzyme" is a string for the specific protease used,
        "length" is the returned number of amino acids, and
        "ntt" is the number of tryptic termini.

    Written by Phil Wilmarth, OHSU, 2008.

    THIS NEEDS TO BE REWRITTEN!!!
    """
    import re
    # supported enzymes are: 'Tryp', 'GluC', 'AspN', and 'LysC'
    # 
    # This routine removes bounding amino acids, removes special characters
    # (mods), is now case sensitive, and computes number of enzymatic termini.
    # Assumes periods are used to separate bounding AAs from peptide.
    # Bounding AAs can be more than one character ("-" for N-term or C-term).
    # Modifications are Comet/SEQUEST format: special characters, "n", and "c";
    #   and are within the bounding periods (if present).
    #
    # Fixed bug in ntt caclulation, 4/30/07 -PW
    # Added support for different enzymes, 7/6/2010 -PW
    # Supports Comet PTM format ("n" and "c" for termini), 6/9/2015 -PW
    # Simplified ntt calculations
    
    # find the string between the bounding periods '.'
    parts = len(sequence_string.split('.'))
    if parts == 3: # we have bounding residues
        start = sequence_string.index('.') + 1   # start is after first period
        temp = sequence_string[::-1] # reverse string
        end = temp.index('.')+1     # find first period in reversed string
        end = len(sequence_string) - end     # end is past the last period
    elif parts == 1:
        start = 0
        end = len(sequence_string)
    else:
        print('...amino_acid_count WARNING: number of "periods" was not 2 or 0', sequence_string)
        if return_base_pep:
            return (0, 0, "")
        else:
            return(0, 0)
    sequence = sequence_string[start:end]

    # remove any modification symbols:
    # '*', '#', '@', '^', '~', '$', '%', '!', '+', 'n', 'c', '[', ']' (Current Comet along with old style nt, ct)
    splitter = re.compile(r'[*#@^~$%!+nc\[\]]')
    base_seq = ''.join(splitter.split(sequence))

##    # remove any special characters from the sequence string
##    base_seq = ''
##    for c in sequence:
##        if c.isalpha() and c.isupper():
##            base_seq += c
##    
    # get the prefix and suffix amino acids
    prefix = sequence_string[start-2:start-1]
    if (prefix == "") or (start == 0):
        prefix = "X"    # no bounding residue info so unknown AA
    cterm = base_seq[-1]  # last amino acid in sequence
    nterm = base_seq[0]   # first amino acid in sequence
    suffix = sequence_string[end+1:end+2]
    if suffix == "":
        suffix = "X"    # no bounding residue info so unknown AA

    # determine number of enzymatic termini, ntt
    """need to support all enzymes and deal with proline
    Seems Comet deals with premature stop codons as sequence breaks (* in prefix or suffix)
    """
    ntt = 0
    if enzyme.upper() == 'TRYP':  # cleaves at c-side of K, R
        if (prefix in 'KR-*'):
            ntt += 1
        if (cterm in 'KR') or (suffix in '-*'):
            ntt += 1
## # this interferes with semi-tryptic option in Comet
##        if suffix == 'P' and ntt > 0:   # trypsin strict???
##            ntt -= 1
    elif enzyme.upper() == 'GLUC':  # cleaves at c-side of D, E
        if prefix in 'DE-*':
            ntt += 1
        if (cterm in 'DE') or (suffix in '-*'):
            ntt += 1
    elif enzyme.upper() == 'ASPN': # cleaves at n-side of D
        if (prefix in '-*') or (nterm == 'D'):
            ntt += 1
        if suffix in 'D-*':
            ntt += 1
    elif enzyme.upper() == 'LYSC':
        if (prefix in 'K-*'):
            ntt += 1
        if (cterm in 'K') or (suffix in '-*'):
            ntt += 1        
    else:
        print('   amino_acid_count WARNING: unknown enzyme specified', enzyme)
    
    # return length, number of tryptic termini, and (optional) base peptide sequence
    if return_base_pep:
        return (len(base_seq), ntt, base_seq)
    else:
        return (len(base_seq), ntt)

def get_base_peptide_sequence(sequence, mask=True):
    """Returns the amino acid sequence from SEQUEST peptide sequences
    """
    # get rid of bounding residues, if any
    try:
        peptide = sequence.split('.')[1]
    except IndexError:
        peptide = sequence
    
    # remove any modification symbols and mask I/L:
    # '*', '#', '@', '^', '~', '$', '%', '!', '+', 'n', 'c', '[', ']' (Current Comet along with old style nt, ct)
    splitter = re.compile(r'[*#@^~$%!+nc\[\]]')
    base_pep = ''.join(splitter.split(peptide))
    if mask:
        base_pep = re.sub(r'[IL]', 'j', base_pep)
    return base_pep    

################## end support functions for PAW pipeline use ##########################    


################################# support for histogram GUI stuff ######################
    
class FileLoader:
    """Load all TXT file in folder into a pandas dataframe. 
    Returns dataframe and list of TXT info objects keyed to TxtIdx column.
    Written by Phil Wilmarth, OHSU, 2014. Extended by Billy Rathje, 2014.
    Added support for TXT file without PSR Core instrument designations 7/21/14 -PW
    """
    def __init__(self, file_list):
        self.peptideMassTol = 1.25  # search parent ion tolerance default (plus/minus in Da)
        self.modStrings = []        # list of variable modification symbols specified in search
        self.enzyme = True          # False if a no enzyme search was used (default=True)

        # get the list of TXT files and figure out instrument type
        folder = os.path.dirname(file_list[0])

        # parse the params file and get relevant settings
        self.params = CometParams()
        self.params.load_from_folder(folder)
        self.peptideMassTol = self.params.peptide_mass_tolerance
        
        # Construct list of differential mods
        self.modStrings = self.generateModStrings()
        
        self.enzyme = self.parseEnzymeInfo()
        
        self.fileList = file_list
    
##        # pick columns to use (drop ISBDisc column) and their data types
##        use_cols = ['start', 'end', 'Z', 'expM', 'theoM', 'SpRank', 'Xcorr', 'Sequence',
##                    'Loci', 'deltaCN', 'NewDeltaCN', 'NewDisc', 'ntt', 'ForR']
##        col_types = {'start':np.int32, 'end':np.int32, 'Z':np.int32, 'expM':np.float64,
##                    'theoM':np.float64, 'SpRank':np.int32, 'Xcorr':np.float64, 'Sequence':'object',
##                    'deltaCN':np.float64, 'Loci':'object', 'NewDeltaCN':np.float64, 'NewDisc':np.float64,
##                    'ntt':np.int32, 'ForR':'object'}
    
        # loop over files and read into pandas DataFrames (save in list)
        self.txt_info_list = []
        self.frame = pd.DataFrame()
        print('\nProcessing all TXT files in:', os.path.basename(folder), time.asctime())
        for i, file_name in enumerate(file_list):
            file_obj = open(file_name, 'rU')
            name = os.path.basename(file_name)[:-4]
##            frame = pd.read_csv(file_obj, sep='\t', usecols=use_cols, dtype=col_types)
            frame = pd.read_csv(file_obj, sep='\t')
            info = TxtInfo(name)
            self.txt_info_list.append(info)
            frame['TxtIdx'] = i
            
            # save the frame's contents
            self.frame = pd.concat([self.frame, frame])
            print('...%s had %s lines' % (info.basename, len(frame)))

        print('...%s total lines read in' % len(self.frame))

    def getFrame(self):
        return self.frame
        
    def getParams(self):
        return self.params

    def getPeptideMassTol(self):
        return self.peptideMassTol
                
    def generateModStrings(self):
        mod_list = []
        last_deltamass_seen = 0
        for k, v in self.params.variable_mods.items():
            if v[1] == 'X':    # skip invalid residues
                continue
            elif float(v[0]) == 0.0:     # skip if deltamass equals zero
                continue
            last_deltamass_seen = '%+0.4f' % float(v[0])
            mod_list.append(v[1] + last_deltamass_seen)
        return mod_list
            
    def parseEnzymeInfo(self):
        if self.params.search_enzyme_number == 0:
            return False
        else:
            return True
                                                           
class Plot:
    ''' Plot is a class for individual histograms. It keeps histogram info (counts and bins) as well as
        running remainders, fdrs, etc. There are several plots for each FigureGenerator object. '''
    def __init__(self, z, ntt, data):
        self.mod = None
        self.dm = 0
        self.threshold = 0.0
        self.z = z
        self.ntt = ntt
        self.data = data  
        self.histo = pd.DataFrame()    
        self.forward, self.reverse = self.make_histograms(data[data.ForR == 'F'].drop_duplicates(['start', 'end', 'Z', 'TxtIdx']), 
                                                          data[data.ForR == 'R'].drop_duplicates(['start', 'end', 'Z', 'TxtIdx']))
        self.generateColumns(self.forward, self.reverse)
        self.Smforward, self.Smreverse = self.smooth(self.forward), self.smooth(self.reverse)
    
    def make_histograms(self, forward, reverse):
        # Plots
        mybins = np.linspace(-8, 12, 301)
        counts, bins = np.histogram(forward['NewDisc'], bins=mybins)
        rcounts, bins = np.histogram(reverse['NewDisc'], bins=mybins)
        return counts, rcounts
        
    def generateColumns(self, counts, rcounts):
        self.histo['DiscScore'] = np.linspace(-8, 12, 301)[0:300]  # This matches Phil's spreadsheets.
        self.histo['Forward'] = counts
        self.histo['Reverse'] = rcounts
        self.histo['RRForward'] = self.runningRemainder(counts)
        self.histo['RRReverse'] = self.runningRemainder(rcounts)
        self.histo['SmRRForward'] = self.runningRemainder(self.smooth(counts))
        self.histo['SmRRReverse'] = self.runningRemainder(self.smooth(rcounts))
        self.histo['FDR'] = (self.histo['RRReverse'] / self.histo['RRForward']) * 100.0
        self.histo['SmFDR'] = ((self.histo['SmRRReverse'] / self.histo['SmRRForward']) * 100.0)
        
        # Convert to float format
        self.histo.Forward = self.histo.Forward.map('{:.2f}'.format)
        self.histo.Reverse = self.histo.Reverse.map('{:.2f}'.format)        
        self.histo.SmRRForward = self.histo.SmRRForward.map('{:.2f}'.format)
        self.histo.SmRRReverse = self.histo.SmRRReverse.map('{:.2f}'.format)
        self.histo.FDR = self.histo.FDR.map('{:.2f}'.format)
        self.histo.SmFDR = self.histo.SmFDR.map('{:.2f}'.format)

    def runningRemainder(self, lst):
        retlst = []
        total = np.sum(lst)
        run = 0.0
        for x in np.nditer(lst):
            run += x
            retlst.append(total - run)    
        return np.array(retlst)
        
    def fdr(self, RRForward, RRReverse):
        return (RRReverse/RRForward) * 100
        
    def smooth(self, x, window_len=11, window='hanning'):
        # taken from: http://wiki.scipy.org/Cookbook/SignalSmooth
        """smooth the data using a window with requested size.
        
        This method is based on the convolution of a scaled window with the signal.
        The signal is prepared by introducing reflected copies of the signal 
        (with the window size) in both ends so that transient parts are minimized
        in the begining and end part of the output signal.
        
        input:
            x: the input signal 
            window_len: the dimension of the smoothing window; should be an odd integer
            window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
                flat window will produce a moving average smoothing.    
        output:
            the smoothed signal
            
        example:    
        t=linspace(-2,2,0.1)
        x=sin(t)+randn(len(t))*0.1
        y=smooth(x)
        
        see also:         
        numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
        scipy.signal.lfilter
    
        TODO: the window parameter could be the window itself if an array instead of a string
        NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
        """    
        if x.ndim != 1:
            raise ValueError("smooth only accepts 1 dimension arrays.")  
        if x.size < window_len:
            raise ValueError("Input vector needs to be bigger than window size.")      
        if window_len < 3:
            return x      
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
            raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
    
        s = np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]     # s is longer than x by 2*(window_len-1)
        if window == 'flat': #moving average
            w = np.ones(window_len,'d')
        else:
            w = eval('np.'+window+'(window_len)')
    
        y = np.convolve(w/w.sum(),s,mode='valid')    # y is trimmed by half window width (beginning and ending are skipped)
        trim = (window_len -1)//2    # added by PW. Original returned a longer vector
        return y[trim:-trim]   
        
class DeltaMassPlot(Plot):
    """Do we need some mechanism to use PPM instead of Da?
    """
    def __init__(self, smoothed, z, dm, dmData, ACCURATE_MASS, dmRange):
        self.smoothed = smoothed
        self.z = z
        self.dm = dm
        self.dmRange = dmRange
        self.mybins = None          # Number of bins
        self.histo = pd.DataFrame() # full range deltamass tables
        self.ACCURATE_MASS = ACCURATE_MASS
        self.LOW_MASS_BINS = 1000
                                   
        self.forwardDeltaMass, self.reverseDeltaMass =  self.make_histograms(dmData[(dmData.ForR == 'F') & (dmData.Z == self.z)].drop_duplicates(['start', 'end', 'Z', 'TxtIdx']), 
                                                                             dmData[(dmData.ForR == 'R') & (dmData.Z == self.z)].drop_duplicates(['start', 'end', 'Z', 'TxtIdx']),
                                                                             -self.dmRange, self.dmRange)
        self.smForwardDeltaMass, self.smReverseDeltaMass = self.smooth(self.forwardDeltaMass), self.smooth(self.reverseDeltaMass)
        
        if self.dm == 0:
            self.forwardDeltaMassZero, self.reverseDeltaMassZero =  self.make_histograms(dmData[(dmData.ForR == 'F') &
                                                                                                (dmData.Z == self.z) &
                                                                                                (dmData.dmassDa >= -0.05) & 
                                                                                                (dmData.dmassDa <= 0.05)].drop_duplicates(['start', 'end', 'Z', 'TxtIdx']), 
                                                                                        dmData[(dmData.ForR == 'R') &
                                                                                                (dmData.Z == self.z) &
                                                                                                (dmData.dmassDa >= -0.05) & 
                                                                                                (dmData.dmassDa <= 0.05)].drop_duplicates(['start', 'end', 'Z', 'TxtIdx']),
                                                                                                -0.05, 0.05)
            self.smForwardDeltaMassZero, self.smReverseDeltaMassZero = self.smooth(self.forwardDeltaMassZero), self.smooth(self.reverseDeltaMassZero)
        
        if self.dm == 1:    
            self.forwardDeltaMassOne, self.reverseDeltaMassOne =  self.make_histograms(dmData[(dmData.ForR == 'F') &
                                                                                            (dmData.Z == self.z) &
                                                                                            (dmData.dmassDa >= 0.90) & 
                                                                                            (dmData.dmassDa <= 1.10)].drop_duplicates(['start', 'end', 'Z', 'TxtIdx']), 
                                                                                    dmData[(dmData.ForR == 'R') &
                                                                                            (dmData.Z == self.z) &
                                                                                            (dmData.dmassDa >= 0.90) & 
                                                                                            (dmData.dmassDa <= 1.10)].drop_duplicates(['start', 'end', 'Z', 'TxtIdx']),
                                                                                            0.90, 1.10)                                                                                           
            self.smForwardDeltaMassOne, self.smReverseDeltaMassOne = self.smooth(self.forwardDeltaMassOne), self.smooth(self.reverseDeltaMassOne)
            
        self.generateColumns()
        self.thresholdLow = None
        self.thresholdHigh = None
        if(dm == 2):
            self.thresholds = self.findPeakWindows()         
        elif(ACCURATE_MASS):
            self.thresholdLow, self.thresholdHigh = self.findPeakWindows()
        else:
            self.dm = 'ALL'           
                  
        # Truncate data frame for these dm values - hard coded for now
#        """Maybe we should make copies of the filtered frames and plot from those?
#        """          
        #if dm == 0 and ACCURATE_MASS:
        #    self.histo = self.histo[(self.histo.deltaMass >= -0.05) & (self.histo.deltaMass <= 0.05)]
        #    #self.zero_histo = self.histo[(self.histo.deltaMass >= -0.05) & (self.histo.deltaMass <= 0.05)].copy()
        #    self.histo = self.histo.reset_index(drop=True)
        #    #self.histo = self.histo.ReverseDeltaMass[(self.histo.ReverseDeltaMass >= -0.04) & (self.histo.ReverseDeltaMass <= 0.04)]
        #if dm == 1:
        #    self.histo = self.histo[(self.histo.deltaMass >= 0.9) & (self.histo.deltaMass <= 1.1)]
        #    #self.one_histo = self.histo[(self.histo.deltaMass >= 0.9) & (self.histo.deltaMass <= 1.1)].copy()
        #    self.histo = self.histo.reset_index(drop=True)
        #    #self.histo = self.histo.ReverseDeltaMass[(self.histo.ReverseDeltaMass >= 0.09) & (self.histo.ReverseDeltaMass <= 1.1)]


    def make_histograms(self, forward, reverse, low, high):
        # Plots
        if not self.ACCURATE_MASS:
            self.mybins = np.linspace(-self.dmRange, self.dmRange, self.LOW_MASS_BINS+1)
            mybins = self.mybins
        else:
            self.mybins = np.linspace(-self.dmRange, self.dmRange, ((2*self.dmRange)/0.0005)+1)
            mybins = np.linspace(low, high, ((high-low)/0.0005)+1)
        counts, bins = np.histogram(forward['dmassDa'], bins=mybins)
        rcounts, bins = np.histogram(reverse['dmassDa'], bins=mybins)
        return counts, rcounts
        
    def make_histograms_ppm(self, forward, reverse):
        """NOTE: ppm units are different scale than Da. Left, right, and bin width
        would need to be different. Something like -50 to +50 in 0.01 chunks for accurate mass
        and -500 to +500 in 1.0 chunks for low-res data.
        """
        # Plots
        if not self.ACCURATE_MASS:
#            self.mybins = np.linspace(-self.dmRange, self.dmRange, self.LOW_MASS_BINS+2)
            self.dmRange = 500.0 # may want another variable here
            self.mybins = np.linspace(-self.dmRange, self.dmRange, (2*self.dmrange)+1)
        else:
            self.dmRange = 50.0 # may want another variable here
            self.mybins = np.linspace(-self.dmRange, self.dmRange, ((2*self.dmRange)/0.01)+1)
        counts, bins = np.histogram(forward['dmassPPM'], bins=self.mybins)
        rcounts, bins = np.histogram(reverse['dmassPPM'], bins=self.mybins)
        return counts, rcounts
                   
    def generateColumns(self):
        """May need to add PPM support here
        """
        if not self.ACCURATE_MASS:
            self.histo['deltaMass'] = np.linspace(-self.dmRange, self.dmRange, self.LOW_MASS_BINS+1)[:-1]
        else:
            self.histo['deltaMass'] = np.linspace(-self.dmRange, self.dmRange, ((2*self.dmRange)/0.0005)+1)[:-1]
        self.histo['ForwardDeltaMass'] = self.forwardDeltaMass
        self.histo['ReverseDeltaMass'] = self.reverseDeltaMass
        self.histo['SmForwardDeltaMass'] = self.smForwardDeltaMass
        self.histo['SmReverseDeltaMass'] = self.smReverseDeltaMass
        
    def findPeakWindows(self, thresh=300, width=5):
        offset = -self.dmRange
        gain = (2*self.dmRange) / (((2*self.dmRange)/0.0005)+1)   # Two times the range / number of bins
        
        if(self.smoothed):
            target_histo = self.smForwardDeltaMass
            decoy_histo = self.smReverseDeltaMass
        else:
            target_histo = self.forwardDeltaMass
            decoy_histo = self.reverseDeltaMass
        mid = target_histo.size//2
        cent = target_histo[mid-200:mid+200]
        cent = int(np.where(cent == target_histo[mid-200:mid+200].max())[0][0] + (mid-200))
        start = cent - 50
        stop = cent + 50
        low = start
        high = stop
        for idx in range(stop-start+1):
            i = start + idx
            t = target_histo[i-width:i+width].sum()
            d = decoy_histo[i-width:i+width].sum()
            try:
                rel = 100.0*t/float(d)
            except ZeroDivisionError:
                rel = 100.0
            if rel <= thresh and i <= cent:
                low = i
            if i > cent and rel >= thresh:
                high = i
            #dm = (i + 0.5)*histo.gain + histo.offset
    ##        print('%i %0.4f %0.2f %0.2f %0.2f' % (i, dm, t, d, rel))
    
        # prevent long tails by invoking a maximum asymmetry
        short = min([abs(cent-low), abs(cent-high)])
        if (cent-low) > 2*short:
            low = cent - 2*short
        if (high-cent) > 2*short:
            high = cent + 2*short
        #print('For charge %s+:' % (self.z+1,))
        net = target_histo[low:high].sum() - decoy_histo[low:high].sum()
        #print('   cent: %s, low: %s, high: %s, net: %0.0f' % (cent, low, high, net))    
    
        # print delta mass windows in Da
        cent = (cent+0.5)*gain + offset     # calculate at center of bin
        low = (low)*gain + offset           # calculate at left of bin
        high = (high+1.0)*gain + offset     # calculate at right of bin
        left_mass_width = cent - low
        right_mass_width = high - cent
        #print('   cent: %0.4f, low: %0.4f, high: %0.4f' % (cent, low, high))
        if self.dm == 0:
            return (low, high)
        elif self.dm == 1:
            # +1 Da regions
            one_low = cent + 0.9840 - left_mass_width
            one_high = cent + 1.0033 + right_mass_width
            return (one_low, one_high)
        elif self.dm == 2:
            one_low = cent + 0.9840 - left_mass_width
            one_high = cent + 1.0033 + right_mass_width
            return ((-self.dmRange, low), (high, one_low), (one_high, self.dmRange))

class TxtInfo:
    """Container for information and stats about main data.
    """
    def __init__(self, full_txt_name):
        self.path = os.path.dirname(full_txt_name)  # path name of folder containing the TXT files
        self.basename = os.path.basename(full_txt_name) # basename of the TXT file (no extension)
        self.target_top_hits = 0        # total number of target top hits
        self.decoy_top_hits = 0         # total number of decoy top hits
        self.target_scans = 0           # total number of target scans
        self.decoy_scans = 0            # total number of decoy scans
        self.target_matches_top = None     # multi-dimensional collection of counters for target matches
        self.decoy_matches_top = None      # multi-dimensional collection of counters for decoys
        self.target_matches_scan = None     # multi-dimensional collection of counters for target matches
        self.decoy_matches_scan = None      # multi-dimensional collection of counters for decoys
        
class BRTxtInfo:
    """Container for information and stats about each TXT file.
    """
    def __init__(self, full_txt_name):
        """full_txt_name should have file extension removed
        """
        self.path = os.path.dirname(full_txt_name)  # path name of folder containing the TXT files
        self.basename = os.path.basename(full_txt_name) # basename of the TXT file (no extension)
        self.target_top_hits = 0        # total number of target top hits
        self.decoy_top_hits = 0         # total number of decoy top hits
        self.target_scans = 0           # total number of target scans
        self.decoy_scans = 0            # total number of decoy scans
        
        self.target_matches_top = None     # multi-dimensional collection of counters for target matches
        self.decoy_matches_top = None      # multi-dimensional collection of counters for decoys
        self.target_matches_scan = None     # multi-dimensional collection of counters for target matches
        self.decoy_matches_scan = None      # multi-dimensional collection of counters for decoys

        self.target_subclass = None     # multi-dimensional collection of counters for target matches
        self.decoy_subclass = None      # multi-dimensional collection of counters for decoys
        self.target_filtered = 0
        self.decoy_filtered = 0
        self.min_length = 7             # minimum peptide length (should be passed in)
        self.maxMods = 3               # maximum number of mods per peptide (should be passed in)
 
    def getStats(self, frame, dm_list, z_list, ntt_list, mod_list, masses, scores):
        """Computes several stats on numbers of target and decoy matches
        masses are the deltamass windows
        scores are the score thresholds
        """
        # restrict by minimunm peptide length first
        len_frame = frame[(frame.Length >= self.min_length) & (frame.NumMods <= self.maxMods)]
        
        # get the global stats first
        self.target_top_hits = len(len_frame[len_frame.ForR == 'F'])
        self.decoy_top_hits = len(len_frame[len_frame.ForR == 'R'])
        #   print(self.basename, self.target_top_hits, self.decoy_top_hits, self.target_top_hits - self.decoy_top_hits)
        self.target_scans = len(len_frame[len_frame.ForR == 'F'].drop_duplicates(['start', 'end', 'Z']))
        self.decoy_scans = len(len_frame[len_frame.ForR == 'R'].drop_duplicates(['start', 'end', 'Z']))
        #   print(self.basename, self.target_scans, self.decoy_scans, self.target_scans - self.decoy_scans)
    
        # get the peptide subclass stats
        #print(self.basename)
        for dm in range(len(dm_list)):
            mass_frame = len_frame[(len_frame.dmassDa >= masses[dm]['low']) &
                                (len_frame.dmassDa <= masses[dm]['high'])]
            #mass_frame = len_frame
            for z in range(len(z_list)):    # z is one less than the charge state
                for ntt in range(len(ntt_list)):
                    for mod in range(len(mod_list)):
                        subclass_frame = mass_frame[(mass_frame.Z == z+1) &
                                                    (mass_frame.ntt == ntt) &
                                                    (mass_frame.ModsStr == mod_list[mod])]
                                 
                        s = scores[dm][z][ntt][mod]
                        threshold = s.histo.DiscScore[s.threshold]
                            
                        self.target_subclass[dm][z][ntt][mod] = len(subclass_frame[(subclass_frame.ForR == 'F') &
                                                                                (subclass_frame.NewDisc >= threshold)])
                        self.target_filtered += len(subclass_frame[(subclass_frame.ForR == 'F') &
                                                                (subclass_frame.NewDisc >= threshold)])
                        self.decoy_subclass[dm][z][ntt][mod] = len(subclass_frame[(subclass_frame.ForR == 'R') &
                                                                                (subclass_frame.NewDisc >= threshold)])
                        self.decoy_filtered += len(subclass_frame[(subclass_frame.ForR == 'R') &
                                                                (subclass_frame.NewDisc >= threshold)])
            #        print(self.target_subclass[dm][z][ntt][mod], z, ntt, mod)
            #print(self.target_filtered, self.decoy_filtered, self.target_filtered - self.decoy_filtered)
            
class FigureGenerator:
    ''' FigureGenerator operates on a set of TXT file. It generates extra calculated columns and makes several histograms
        based on instrument type and user parameters. Each FigureGenerator object should keep a list (figures) of several
        Plot objects corresponding to individual histograms. This list is passed to the GUI object for display.
        '''
    def __init__(self, files, accurateMass=True, smoothed=True, dmList=['0 Da', '1 Da', 'out'], zList=[1, 2, 3, 4],
                 nttList=[0, 1, 2], modString=' *', minLength=7, maxMods=2):
        
        
        # Main containers
        self.minLength = minLength
        self.maxMods = maxMods
        self.globalStats = None
        self.txtStats = []
        self.txtObjects = []

        # Loading in file attributes
        f = FileLoader(files)                   # load files
        self.f = f                               # pointer to initial file object
        x = self.f.getFrame()
        self.accurateMass = accurateMass
        self.smoothed = smoothed
        self.peptideMassTol = self.f.getPeptideMassTol()
        
        self.order = [' ', '*', '#', '@', '^', '~', '$', '[', ']']
        
        # Main container dimensions
        dmList, zList, nttList = self.getLists(self.accurateMass, self.f.enzyme)
        
        # setup mod list
        modList = [c for c in self.f.modStrings]
        modList.insert(0, 'Unmodified')
        
        # Make a list of special characters present in the dataset in addition
        # to the list of mod strings (eg. M+15.99...)
        specialCharsList = []
        for i, mod in enumerate(modList):
            if mod.startswith('nt'):
                specialCharsList.append(']')
                continue
            elif mod.startswith('ct'):
                specialCharsList.append('[')
                continue
            else:
                specialCharsList.append(self.order[i])
            
        #self.modDict = {}
        #for i, mod in enumerate(ModList):
        #    if mod.startswith('nterm'):
        #        modDict['['] = mod
        #    elif mod.startswith('cterm'):
        #        modDict[']'] = mod
        #    else:
        #        modDict[self.order[i]] = mod
        
        # check lists
#        print('enzyme check:', self.f.enzyme)
        print('\ndm list:', dmList)
        print('Z list:', zList)
        print('NTT list:', nttList)
        print('Mod List:', modList)
        print('Special chars:', specialCharsList, '\n')
        
        self.container = [[[[None for mod in modList] for ntt in nttList] for z in zList] for dm in dmList]
        self.dmContainer = [[None for z in zList] for dm in dmList]
        
        # Add calculated columns
        self.generateCalculatedColumns(x)
        
        # Make histograms
        self.generatePlots(x, dmList, zList, nttList, modList)
        
        # Store these specifically for calculating stats later
        self.modList = modList
        self.specialCharsList = specialCharsList
        self.nttList = nttList
        self.zList = zList
        self.dmList = dmList
        self.frame = x
                
        # Test
        print("Calculating stats...")
        #self.get_stats_helper()
        self.getTXTFileObjectsHelper()
                
    def generatePlots(self, full_frame, dmList, zList, nttList, modList):
        """Generates the histogram plot data.
        full_frame => pandas dataframe of TXT file contents
        dmList => list of strings describing the delta mass windows
        zList => list of integers spanning the charge state range (continuous range)
        nttList => list of integers spanning the ntt range (continous range)
        """
        z_offset = min([int(z) for z in zList])
        ntt_offset = min([int(ntt) for ntt in nttList])
        data = full_frame[(full_frame.Length >= self.minLength) & (full_frame.NumMods <= self.maxMods)]   # apply global constraints first
        for dm in range(len(dmList)):
            print('Loading: ', dmList[dm], '...')
            for z in zList:
                print('\tLoading: %d+ ...' % z)                
                dmPlot =  DeltaMassPlot(self.smoothed, z, dm, data, self.accurateMass, self.peptideMassTol)
                if not self.accurateMass:
                    dmPlot.thresholdLow = -self.peptideMassTol
                    dmPlot.thresholdHigh = self.peptideMassTol
                self.dmContainer[dm][z-z_offset] = dmPlot
                #for ntt in nttList:
                #    print('\t\t Loading Ntt =', ntt, '...')
                #    for mod in range(len(modList)):
                #        if dm == 2:
                #            dataToHist = data[((data.dmassDa >= dmPlot.thresholds[0][0]) & (data.dmassDa <= dmPlot.thresholds[0][1]) |
                #                              (data.dmassDa >= dmPlot.thresholds[1][0]) & (data.dmassDa <= dmPlot.thresholds[1][1]) |
                #                              (data.dmassDa >= dmPlot.thresholds[2][0]) & (data.dmassDa <= dmPlot.thresholds[2][1])) &
                #                              (data.Z == z) & (data.ntt == ntt) & (data.ModsStr == self.order[mod])]
                #        else:
                #            dataToHist = data[(data.dmassDa >= dmPlot.thresholdLow) & (data.dmassDa <= dmPlot.thresholdHigh) &
                #                              (data.Z == z) & (data.ntt == ntt) & (data.ModsStr == self.order[mod])]
                #                                
                #        self.container[dm][z-z_offset][ntt-ntt_offset][mod] = Plot(z, ntt, dataToHist)
                #                
                #        self.container[dm][z-z_offset][ntt-ntt_offset][mod].dm = dmList[dm]
                #        self.container[dm][z-z_offset][ntt-ntt_offset][mod].mod = modList[mod]
                        
    def regenerateScorePlots(self):
        z_offset = min([int(z) for z in self.zList])
        ntt_offset = min([int(ntt) for ntt in self.nttList])
        
        data = self.frame[(self.frame.Length >= self.minLength) & (self.frame.NumMods <= self.maxMods)] 
        
        for dm in range(len(self.dmList)):
            for z in self.zList:
                for ntt in self.nttList:
                    for imod, mod in enumerate(self.specialCharsList):
                        if dm == 2:
                            dataToHist = data[((data.dmassDa >= self.dmContainer[2][z-z_offset].thresholds[0][0]) & (data.dmassDa <= self.dmContainer[2][z-z_offset].thresholds[0][1]) |
                                              (data.dmassDa >= self.dmContainer[2][z-z_offset].thresholds[1][0]) & (data.dmassDa <= self.dmContainer[2][z-z_offset].thresholds[1][1]) |
                                              (data.dmassDa >= self.dmContainer[2][z-z_offset].thresholds[2][0]) & (data.dmassDa <= self.dmContainer[2][z-z_offset].thresholds[2][1])) &
                                              (data.Z == z) & (data.ntt == ntt) & (data.ModsStr == mod)]
                        else:
                            dataToHist = data[(data.dmassDa >= self.dmContainer[dm][z-z_offset].thresholdLow) & (data.dmassDa <= self.dmContainer[dm][z-z_offset].thresholdHigh) &
                                              (data.Z == z) & (data.ntt == ntt) & (data.ModsStr == mod)]
                                                
                        self.container[dm][z-z_offset][ntt-ntt_offset][imod] = Plot(z, ntt, dataToHist)
                                
                        self.container[dm][z-z_offset][ntt-ntt_offset][imod].dm = self.dmList[dm]
                        self.container[dm][z-z_offset][ntt-ntt_offset][imod].mod = self.modList[imod]
            
           
    def get_stats_helper(self): 
        self.get_stats(self.frame, self.f.fileList, self.dmList, self.zList, self.nttList, self.modList, [{'low':-2, 'high':2}, {'low':-2, 'high':2}, {'low':-2, 'high':2}], self.container)
    
    def getLists(self, accurateMass, enzyme):
        # order:   dm, z. ntt
        nttList = []
        if enzyme:
            nttList = [1, 2]
        else:
            nttList = [0, 1, 2]
        if accurateMass:
            return ['0 Da', '1 Da', 'out'], [2, 3, 4], nttList
        else:
            return ['All'], [1, 2, 3], nttList
    
    def getTXTFileObjectsHelper(self): 
        self.getTXTFileObjects(self.frame, self.f.fileList, self.dmList, self.zList, self.nttList, self.modList, [{'low':-2, 'high':2}, {'low':-2, 'high':2}, {'low':-2, 'high':2}], self.container)
            
    def getTXTFileObjects(self, frame, fileList, dmList, zList, nttList, modList, masses, container):
        for i, filename in enumerate(fileList):
            fileFrame = frame[(frame['TxtIdx'] == i)]
            t = BRTxtInfo(filename.split('.txt')[0])              # object for global stats and stats by text file
            self.txtObjects.append(t)
      
    def get_stats(self, frame, fileList, dmList, zList, nttList, modList, masses, container):
        # global stats
        self.globalStats = BRTxtInfo(self.f.folder)              # object for global stats and stats by text file
        # multidimensional containers
        self.globalStats.target_subclass = [[[[None for mod in modList] for ntt in nttList] for z in zList] for dm in dmList]
        self.globalStats.decoy_subclass = [[[[None for mod in modList] for ntt in nttList] for z in zList] for dm in dmList]
        self.globalStats.getStats(frame, dmList, zList, nttList, modList, masses, container)

        #for i, filename in enumerate(fileList):
        #    fileFrame = frame[(frame['TxtIdx'] == i)]
        #    t = BRTxtInfo(filename)              # object for global stats and stats by text file
        #    # multidimensional containers
        #    t.target_subclass = [[[[None for mod in modList] for ntt in nttList] for z in zList] for dm in dmList]
        #    t.decoy_subclass = [[[[None for mod in modList] for ntt in nttList] for z in zList] for dm in dmList]
        #    self.globalStats.getStats(fileFrame, dmList, zList, nttList, modList, masses, container)
        #    self.txtStats.append(t)
                        
    def debugPlot(self):
        for fig in self.container[0]:
            for fig in fig[1:]:
                fig = fig[0]
                mybins = np.linspace(-8, 12, 301)
                # found the center calculation online, may want to double check. corrects x axis
                center = (mybins[:-1] + mybins[1:]) / 2
                plot.fill(center, fig.forward, color='b')
                plot.fill(center, fig.reverse, color='r')
                plot.title(str(fig.z) + "+ charge state and " + str(fig.ntt) + " tryptic termini.")
                plot.figure()
##        plot.show() # this may be out-of-date
        plot.draw()
    
    def generateCalculatedColumns(self, x):
        # generate calculated columns
        
        # deltaMasses (Da and ppm))
        x['dmassDa'] = x['expM'] - x['theoM']
        x['dmassPPM'] = (10000000*(x['expM'] - x['theoM'])  / (np.sqrt(x['expM'] * (x['theoM'] ))))
        
        # peptide length
        temp = x['Sequence'].str.split('.')

        # count number of amino acids in peptide
        x['Length'] = temp.map(lambda y: ''.join([c for c in str(y[1]) if c.isalpha()]))
        x['Length'] = x['Length'].str.len()

        # count number of variable mods in peptide
        x['NumMods'] = temp.map(lambda y: ''.join([c for c in str(y[1]) if not c.isalpha()]))
        x['NumMods'] = x['NumMods'].str.len()

        # get list of modification symbols in order (leading space for unmodified amino acids)
        x['ModsStr'] = temp.map(lambda y: self.mods_str(y[1])) # with luck, mods are in order with leading space

    def mods_str(self, seq):
        """Returns an odered list of modification types present in peptide.
        Includes a space character for unmodified amino acids.
        """
        s = list(set([char for char in seq if not char.isalpha()]))
        s = ''.join([x for x in self.order if x in s])
        if not s:
            s = ' '
        return s     # want leading space character for unmodified amino acids  
        
class DataInfoAndFilter:
    """Container for global stats on all loaded data.
    Uses TxtInfo objects to keep stats for each TXT file.
    Has aggregate counters and totaling method. Counters
    track both before filtering and post filtering counts.
    """    
    def __init__(self, folder, frame, txt_info_list, dm_list, z_list, ntt_list, mod_list, 
                  min_length=7, max_mods=2, parent_tol=2.5):
        import copy
        
        # main data passed in
        self.folder = folder    # full path to TXT files
        self.frame = frame      # pandas dataframe of TXT file contents and some extras        
        self.pre_filter = txt_info_list                              # list of TxtInfo objects for pre-filter stats
        self.post_filter = [copy.deepcopy(x) for x in txt_info_list] # list of TxtInfo objects for post-filter stats
        self.dm_list = dm_list      # list of delta mass window names
        self.z_list = z_list        # list of allowed charge states (contiguous range)
        self.ntt_list = ntt_list    # list of number of tryptic termini (o, 1, 2)
        self.mod_list = mod_list    # full, ordered list of variable modification symbols starting with a space for unmodified residues

        # some restrictions and limits        
        self.min_length = min_length        # minimum peptide length (should be passed in)
        self.max_mods = max_mods            # maximum number of mods per peptide (should be passed in)
        self.z_offset = min([int(z) for z in z_list]) # to map peptide charge to z-axis index
        self.z_max = max([int(z) for z in z_list])    # maximum peptide charge
        self.ntt_offset = min([int(ntt) for ntt in ntt_list]) # to map ntt range to index range
        self.ntt_max = max([int(ntt) for ntt in ntt_list])    # maximum ntt value
        self.parent_tol = parent_tol        # parent ion tolerance (should be passed in)
        
        # data structures for counting statistics
        self.short_target_top_hits = 0      # total number of target top hits below min length
        self.short_decoy_top_hits = 0       # total number of decoy top hits below min length
        self.short_target_scans = 0         # total number of target scans below min length
        self.short_decoy_scans = 0          # total number of decoy scans below min length
        self.pre_target_top_hits = 0        # total number of target top hits
        self.pre_decoy_top_hits = 0         # total number of decoy top hits
        self.pre_target_scans = 0           # total number of target scans
        self.pre_decoy_scans = 0            # total number of decoy scans
        self.pre_target_matches_top = None  # multi-dimensional collection of counters for target matches
        self.pre_decoy_matches_top = None   # multi-dimensional collection of counters for decoys
        self.pre_target_matches_scan = None # multi-dimensional collection of counters for target matches
        self.pre_decoy_matches_scan = None  # multi-dimensional collection of counters for decoys
        self.post_target_top_hits = 0       # total number of target top hits
        self.post_decoy_top_hits = 0        # total number of decoy top hits
        self.post_target_scans = 0          # total number of target scans
        self.post_decoy_scans = 0           # total number of decoy scans
        self.post_target_matches_top = None # multi-dimensional collection of counters for target matches
        self.post_decoy_matches_top = None  # multi-dimensional collection of counters for decoys
        self.post_target_matches_scan = None# multi-dimensional collection of counters for target matches
        self.post_decoy_matches_scan = None # multi-dimensional collection of counters for decoys

        # set up for log file
##        self.sqt_container = os.path.dirname(self.folder)
##        self.filtered_folder = os.path.join(self.sqt_container, 'filtered_files')
        self.filtered_folder = os.path.join(os.path.dirname(self.folder), 'filtered_files')
        self.log_file = open(os.path.join(self.folder, os.path.basename(self.folder) + '_PAW.log'), 'a')
        self.write = [None, self.log_file]

        return

    def aggregate_pre_global_stats(self):
        """Sums up the per TXT stats for pre-filtered data
        """
        dm_list = ['All']   # we don't care about mass windows for aggregate stats
        
        # initialize counters and counter containers
        self.pre_target_top_hits = 0
        self.pre_decoy_top_hits = 0
        self.pre_target_scans = 0
        self.pre_decoy_scans = 0
        self.pre_target_matches_top = np.array([[[[0 for mod in self.mod_list] for ntt in self.ntt_list] for z in self.z_list] for dm in dm_list])
        self.pre_decoy_matches_top = np.array([[[[0 for mod in self.mod_list] for ntt in self.ntt_list] for z in self.z_list] for dm in dm_list])
        self.pre_target_matches_scan = np.array([[[[0 for mod in self.mod_list] for ntt in self.ntt_list] for z in self.z_list] for dm in dm_list])
        self.pre_decoy_matches_scan = np.array([[[[0 for mod in self.mod_list] for ntt in self.ntt_list] for z in self.z_list] for dm in dm_list])
        
        # sum over the individual TXT file data
        for txt_info in self.pre_filter:
            self.pre_target_top_hits += txt_info.target_top_hits
            self.pre_decoy_top_hits += txt_info.decoy_top_hits
            self.pre_target_scans += txt_info.target_scans
            self.pre_decoy_scans += txt_info.decoy_scans
            self.pre_target_matches_top += txt_info.target_matches_top
            self.pre_decoy_matches_top += txt_info.decoy_matches_top
            self.pre_target_matches_scan += txt_info.target_matches_scan
            self.pre_decoy_matches_scan += txt_info.decoy_matches_scan
        return

    def aggregate_post_global_stats(self):
        """Sums up the per TXT stats after filtering
        """
        dm_list = ['All']   # we don't care about mass windows for aggregate stats

        # initialize counters and counter containers
        self.post_target_top_hits = 0
        self.post_decoy_top_hits = 0
        self.post_target_scans = 0
        self.post_decoy_scans = 0
        self.post_target_matches_top = np.array([[[[0 for mod in self.mod_list] for ntt in self.ntt_list] for z in self.z_list] for dm in self.dm_list])
        self.post_decoy_matches_top = np.array([[[[0 for mod in self.mod_list] for ntt in self.ntt_list] for z in self.z_list] for dm in self.dm_list])
        self.post_target_matches_scan = np.array([[[[0 for mod in self.mod_list] for ntt in self.ntt_list] for z in self.z_list] for dm in self.dm_list])
        self.post_decoy_matches_scan = np.array([[[[0 for mod in self.mod_list] for ntt in self.ntt_list] for z in self.z_list] for dm in self.dm_list])

        # sum across the filtered TXT file data
        for txt_info in self.post_filter:
            self.post_target_top_hits += txt_info.target_top_hits
            self.post_decoy_top_hits += txt_info.decoy_top_hits
            self.post_target_scans += txt_info.target_scans
            self.post_decoy_scans += txt_info.decoy_scans
            self.post_target_matches_top += txt_info.target_matches_top
            self.post_decoy_matches_top += txt_info.decoy_matches_top
            self.post_target_matches_scan += txt_info.target_matches_scan
            self.post_decoy_matches_scan += txt_info.decoy_matches_scan
        return

    def get_pre_stats(self):
        """Computes numbers of target and decoy matches in various sliced and diced ways
        """
        import pprint

        for obj in self.write:
            print('\nCompiling pre-filter stats', time.asctime(), file=obj)
        # get stats on short peptides first
        self.short_target_top_hits = len(self.frame[(self.frame.Length < self.min_length) & (self.frame.ForR == 'F')])
        self.short_decoy_top_hits = len(self.frame[(self.frame.Length < self.min_length) & (self.frame.ForR == 'R')])
        self.short_target_scans = len(self.frame[(self.frame.Length < self.min_length) & (self.frame.ForR == 'F')].drop_duplicates(['start', 'end', 'Z', 'TxtIdx']))
        self.short_decoy_scans = len(self.frame[(self.frame.Length < self.min_length) & (self.frame.ForR == 'R')].drop_duplicates(['start', 'end', 'Z', 'TxtIdx']))
        for obj in self.write:
            print('...short peptide top hits: %s (%s)' % (self.short_target_top_hits, self.short_decoy_top_hits), file=obj)
            print('...short peptide scans: %s (%s)\n' % (self.short_target_scans, self.short_decoy_scans), file=obj)
        
        # restrict by minimum peptide length, maximum number of mods, charge state range and ntt range first
        frame = self.frame[(self.frame.Length >= self.min_length) & (self.frame.NumMods <= self.max_mods) &
                           ((self.frame.Z >= self.z_offset) & (self.frame.Z <= self.z_max)) &
                           ((self.frame.ntt >= self.ntt_offset) & (self.frame.ntt <= self.ntt_max))]

        print('*** pre-filter ***:', len(self.frame), len(frame))
        
        dm_list = ['all']
        for i, txt_info_obj in enumerate(self.pre_filter):
            # create multidimensional counters for each TXT file
            self.pre_filter[i].target_matches_top = np.array([[[[0 for mod in self.mod_list] for ntt in self.ntt_list] for z in self.z_list] for dm in dm_list])
            self.pre_filter[i].decoy_matches_top = np.array([[[[0 for mod in self.mod_list] for ntt in self.ntt_list] for z in self.z_list] for dm in dm_list])  
            self.pre_filter[i].target_matches_scan = np.array([[[[0 for mod in self.mod_list] for ntt in self.ntt_list] for z in self.z_list] for dm in dm_list])
            self.pre_filter[i].decoy_matches_scan = np.array([[[[0 for mod in self.mod_list] for ntt in self.ntt_list] for z in self.z_list] for dm in dm_list])  

            # get the global stats first (this prevents scans from being overcounted; top hit ties can have different ntt or mods)
            self.pre_filter[i].target_top_hits = len(frame[(frame.TxtIdx == i) & (frame.ForR == 'F')])
            self.pre_filter[i].decoy_top_hits = len(frame[(frame.TxtIdx == i) & (frame.ForR == 'R')])
            self.pre_filter[i].target_scans = len(frame[(frame.TxtIdx == i) & (frame.ForR == 'F')].drop_duplicates(['start', 'end', 'Z']))
            self.pre_filter[i].decoy_scans = len(frame[(frame.TxtIdx == i) & (frame.ForR == 'R')].drop_duplicates(['start', 'end', 'Z']))
            y = self.pre_filter[i]
            for obj in self.write:
                print('...%s: all top hits %s (%s), all scans %s (%s), fdr %0.2f' %
                      (y.basename, y.target_top_hits, y.decoy_top_hits, y.target_scans, y.decoy_scans,
                       100.0 * y.decoy_scans / y.target_scans), file=obj)
    
            # get the peptide subclass stats
            for dm in range(len(dm_list)): # different thresholds for each mass window
                for charge in self.z_list:    # different thresholds for each charge, z is one less than the charge state
                    z = int(charge) - self.z_offset
                    mass_frame = frame[(frame.dmassDa >= (-1)*self.parent_tol) & (frame.dmassDa <= self.parent_tol)]
                    for ntt in self.ntt_list: # different thresholds for each number of tryptic termini
                        ntt_idx = int(ntt) - self.ntt_offset
                        for mod in range(len(self.mod_list)):    # and different thresholds for each homogeneous modification state
                            subclass_frame = mass_frame[(mass_frame.TxtIdx == i) & (mass_frame.Z == int(charge)) & (mass_frame.ntt == int(ntt)) &
                                                        (mass_frame.ModsStr.map(lambda s: s[0]) == self.mod_list[mod])] # first get basic data subclass (need first char of ModsStr to get all)
                            self.pre_filter[i].target_matches_top[dm][z][ntt_idx][mod] = len(subclass_frame[(subclass_frame.ForR == 'F')]) 
                            self.pre_filter[i].decoy_matches_top[dm][z][ntt_idx][mod] = len(subclass_frame[(subclass_frame.ForR == 'R')])  
                            self.pre_filter[i].target_matches_scan[dm][z][ntt_idx][mod] += len(subclass_frame[(subclass_frame.ForR == 'F')].drop_duplicates(['start', 'end', 'Z']))
                            self.pre_filter[i].decoy_matches_scan[dm][z][ntt_idx][mod] += len(subclass_frame[(subclass_frame.ForR == 'R')].drop_duplicates(['start', 'end', 'Z']))
 
        # get aggregate stats
        self.aggregate_pre_global_stats()
        for obj in self.write:
            try:
                print('\n...Aggregate top hits: %s (%s) %0.2f' %
                      (self.pre_target_top_hits, self.pre_decoy_top_hits, 100.0*self.pre_decoy_top_hits/self.pre_target_top_hits), file=obj)
            except ZeroDivisionError:
                print('\n...Aggregate top hits: %s (%s) %0.2f' % (self.pre_target_top_hits, self.pre_decoy_top_hits, 0.0), file=obj)
                # raise ZeroDivisionError # will this cause the program to terminate?
            try:
                print('...Aggregate scans: %s (%s) %0.2f\n' %
                      (self.pre_target_scans, self.pre_decoy_scans, 100.0*self.pre_decoy_scans/self.pre_target_scans), file=obj)
            except ZeroDivisonError:
                print('...Aggregate scans: %s (%s) %0.2f\n' % (self.pre_target_scans, self.pre_decoy_scans, 0.0), file=obj)
                # raise ZeroDivisionError
            print('All target subclass matches (all top hits):', file=obj)
            pprint.pprint(self.pre_target_matches_top, stream=obj)
            print('All decoy subclass matches (all top hits):', file=obj)
            pprint.pprint(self.pre_decoy_matches_top, stream=obj)
            print('All net subclass matches (all top hits):', file=obj)
            pprint.pprint(self.pre_target_matches_top - self.pre_decoy_matches_top, stream=obj)
            print('All target subclass matches (scans):', file=obj)
            pprint.pprint(self.pre_target_matches_scan, stream=obj)
            print('All decoy subclass matches (scans):', file=obj)
            pprint.pprint(self.pre_decoy_matches_scan, stream=obj)
            print('All net subclass matches (scans):', file=obj)
            pprint.pprint(self.pre_target_matches_scan - self.pre_decoy_matches_scan, stream=obj)

    def return_threshold(self, index_tuple, z_offset, ntt_offset, mod_list, dm_scores):
        """Returns the largest threshold associated with any modifications.
        (index_tuple) => (Z, ntt, mods_str);
            Z -> charge state; nnt -> number of tryptic termini; mods_str -> list of modification symbols in peptide
        z_offset => maps charge state to z-index,
        mod_list => full, ordered list of variable modifications specified in search
        dm_scores => subset of thresholds for respective deltamass window.
        """
        thresholds = [dm_scores[index_tuple[0]-z_offset][index_tuple[1]-ntt_offset][mod_list.index(mod)] for mod in index_tuple[2]]
        return max(thresholds)

    def copy_params_files(self):
        """Copies any params files to filtered folder."""
        import shutil
        params_list = [p for p in os.listdir(self.folder) if p.endswith('.params')]
        print('params_list:', params_list)
        for param in params_list:
            shutil.copy2(os.path.join(self.folder, param), os.path.join(self.filtered_folder, param))
 
    def filter_with_stats(self, mass_thresholds, score_thresholds):
        """Filters and computes numbers of target and decoy matches passing thresholds over the peptide subclasses
        """
        import pprint

        for obj in self.write:
            print('\nFiltering data and compiling stats', time.asctime(), file=obj)
        scores = np.array(score_thresholds)
        masses = mass_thresholds
        all_filter_frame = pd.DataFrame()

        # print out global limits, various lists, and threshold values
        for obj in self.write:
            print('...Minimum peptide length:', self.min_length, file=obj)
            print('...Maximum number of mods per peptide:', self.max_mods, file=obj)
            print('...DeltaMass list:', self.dm_list, file=obj)
            print('...Charge state list:', self.z_list, file=obj)
            print('...NTT list:', self.ntt_list, file=obj)
            print('...Modifications list:', self.mod_list, file=obj)
            print('...Delta mass windows:', file=obj)
            for i, z in enumerate(self.z_list):
                print('......Z = %s' % z, file=obj)
                for j, dm in enumerate(self.dm_list[:-1]):
                    print('.........DeltaMass %s: %0.4f to %0.4f' % (dm, masses[i][j].low, masses[i][j].high), file=obj)
            print('...Conditional score thresholds:', file=obj)
            pprint.pprint(scores, stream=obj)
            print(file=obj)

        # lets see what columns get loaded from TXT files
        print('Frame columns at start of filtering:')
        for col in self.frame.columns:
            print(col, self.frame[col].dtype)
        
        # restrict by minimum peptide length first
        frame = self.frame[(self.frame.Length >= self.min_length) & (self.frame.NumMods <= self.max_mods) &
                           ((self.frame.Z >= self.z_offset) & (self.frame.Z <= self.z_max)) &
                           ((self.frame.ntt >= self.ntt_offset) & (self.frame.ntt <= self.ntt_max))].copy()
        
        for i, txt_info_obj in enumerate(self.pre_filter):  # loop over all TXT file data
            filter_frame = pd.DataFrame()
            
            # create multidimensional counters for each TXT file
            self.post_filter[i].target_matches_top = np.array([[[[0 for mod in self.mod_list] for ntt in self.ntt_list] for z in self.z_list] for dm in self.dm_list])
            self.post_filter[i].decoy_matches_top = np.array([[[[0 for mod in self.mod_list] for ntt in self.ntt_list] for z in self.z_list] for dm in self.dm_list])  
            self.post_filter[i].target_matches_scan = np.array([[[[0 for mod in self.mod_list] for ntt in self.ntt_list] for z in self.z_list] for dm in self.dm_list])
            self.post_filter[i].decoy_matches_scan = np.array([[[[0 for mod in self.mod_list] for ntt in self.ntt_list] for z in self.z_list] for dm in self.dm_list])  
    
            # get the peptide subclass stats
            for dm in range(len(self.dm_list)): # different thresholds for each mass window
                for charge in self.z_list:    # different thresholds for each charge, z is z_offset less than the charge state
                    z = int(charge) - self.z_offset
                    if dm < 2:  # these are simple inclusive windows
                        mass_frame = frame[(frame.dmassDa >= masses[z][dm].low) &
                                           (frame.dmassDa <= masses[z][dm].high) &
                                           (frame.TxtIdx == i) & (frame.Z == int(charge))].copy()
                    else:   # this is everything outside of the previous two windows
                        mass_frame = frame[(((frame.dmassDa >= masses[z][dm].low) & (frame.dmassDa < masses[z][0].low)) |
                                            ((frame.dmassDa > masses[z][0].high) & (frame.dmassDa < masses[z][1].low)) |
                                            ((frame.dmassDa > masses[z][1].high) & (frame.dmassDa <= masses[z][dm].high))) &
                                           (frame.TxtIdx == i) & (frame.Z == int(charge))].copy()

                    # this finds the most stringent threshold to test heterogeneously modified peptides against
                    mass_frame['IndexTuple'] = list(zip(mass_frame.Z, mass_frame.ntt, mass_frame.ModsStr))
                    mass_frame['TestValue'] = mass_frame.IndexTuple.map(lambda index_tuple: self.return_threshold(index_tuple, self.z_offset, self.ntt_offset, self.mod_list, scores[dm]))

                    for ntt in self.ntt_list: # different thresholds for each number of tryptic termini
                        ntt_idx = int(ntt) - self.ntt_offset
                        for mod in range(len(self.mod_list)):    # and different thresholds for each modification state
                            subclass_frame = mass_frame[(mass_frame.ntt == int(ntt)) & (mass_frame.ModsStr == self.mod_list[mod]) &
                                                        (mass_frame.NewDisc >= mass_frame.TestValue)] # get any data subclass hits above threshold
                            filter_frame = pd.concat([filter_frame, subclass_frame])    # save the filtered peptide subclass frame

                            # save some stats on the filtered subclass hits
                            self.post_filter[i].target_matches_top[dm][z][ntt_idx][mod] = len(subclass_frame[subclass_frame.ForR == 'F']) 
                            self.post_filter[i].decoy_matches_top[dm][z][ntt_idx][mod] = len(subclass_frame[subclass_frame.ForR == 'R'])
                            self.post_filter[i].target_matches_scan[dm][z][ntt_idx][mod] += len(subclass_frame[subclass_frame.ForR == 'F'
                                                                                                          ].drop_duplicates(['start', 'end', 'Z']))
                            self.post_filter[i].decoy_matches_scan[dm][z][ntt_idx][mod] += len(subclass_frame[subclass_frame.ForR == 'R'
                                                                                                         ].drop_duplicates(['start', 'end', 'Z']))
                            
            # get the global stats last from the filtered frame
            self.post_filter[i].target_top_hits = self.post_filter[i].target_matches_top.sum()
            self.post_filter[i].decoy_top_hits = self.post_filter[i].decoy_matches_top.sum()
            self.post_filter[i].target_scans = self.post_filter[i].target_matches_scan.sum()
            self.post_filter[i].decoy_scans = self.post_filter[i].decoy_matches_scan.sum()
            y = self.post_filter[i]
            for obj in self.write:
                print('...%s: filtered top hits %s (%s), filtered scans %s (%s), fdr %0.2f' %
                      (y.basename, y.target_top_hits, y.decoy_top_hits, y.target_scans,
                       y.decoy_scans, 100.0 * y.decoy_scans / y.target_scans), file=obj)
            # write the filtered scans to TXT file
            self.write_data(filter_frame, txt_info_obj.basename)
            all_filter_frame = pd.concat([all_filter_frame, filter_frame]) # we may not need to keep this dataframe in memory

        # get aggregate stats
        self.aggregate_post_global_stats()
        for obj in self.write:
            try:
                print('\n...Aggregate filtered top hits: %s (%s) %0.2f' %
                      (self.post_target_top_hits, self.post_decoy_top_hits, 100.0*self.post_decoy_top_hits/self.post_target_top_hits), file=obj)
            except ZeroDivisionError:
                print('\n...Aggregate filtered top hits: %s (%s) %0.2f' % (self.post_target_top_hits, self.post_decoy_top_hits, 0.0), file=obj)
            try:
                print('...Aggregate filtered scans: %s (%s) %0.2f\n' %
                      (self.post_target_scans, self.post_decoy_scans, 100.0*self.post_decoy_scans/self.post_target_scans), file=obj)
            except ZeroDivisionError:
                print('...Aggregate filtered scans: %s (%s) %0.2f\n' % (self.post_target_scans, self.post_decoy_scans, 0.0), file=obj)
            print('Filtered target subclass matches (all top hits):', file=obj)
            pprint.pprint(self.post_target_matches_top, stream=obj)
            print('Filtered decoy subclass matches (all top hits):', file=obj)
            pprint.pprint(self.post_decoy_matches_top, stream=obj)
            print('Filtered net subclass matches (all top hits):', file=obj)
            pprint.pprint(self.post_target_matches_top - self.post_decoy_matches_top, stream=obj)
            print('Filtered target subclass matches (scans):', file=obj)
            pprint.pprint(self.post_target_matches_scan, stream=obj)
            print('Filtered decoy subclass matches (scans):', file=obj)
            pprint.pprint(self.post_decoy_matches_scan, stream=obj)
            print('Filtered net subclass matches (scans):', file=obj)
            pprint.pprint(self.post_target_matches_scan - self.post_decoy_matches_scan, stream=obj)

        return all_filter_frame

    def write_data(self, frame, basename):
        """Writes a pandas dataframe to a TXT file.
        """
        if not os.path.exists(self.filtered_folder):
            os.mkdir(os.path.join(self.filtered_folder))
        file_name = os.path.join(self.filtered_folder, basename + '_filtered.txt')
##        cols = ['start', 'end', 'Z', 'expM', 'SpRank', 'theoM', 'deltaCN', 'Xcorr', 'Sequence', 'Loci',
##                'NewDeltaCN', 'NewDisc', 'ntt', 'ForR', 'dmassDa', 'dmassPPM', 'Length',
##                'NumMods', 'ModsStr', 'IndexTuple', 'TestValue']
##        frame.to_csv(file_name, sep='\t', index=False, columns=cols)
        frame.to_csv(file_name, sep='\t', index=False)
        for i, obj in enumerate(self.write):
            print('......filtered TXT file written', file=obj)

        # get the list of passing scans and extract the SQT and MS2 information
        scan_list = list(zip(frame.start, frame.end, frame.Z))
        self.write_sqt(basename, scan_list)
        self.write_ms2(basename, scan_list)

    def write_sqt(self, basename, scan_list):
        """Writes a filtered SQT file to filtered_files folder
        """
        # see if there is anything to process
        sqt_out = os.path.join(self.filtered_folder, basename + '_filtered.sqt')
        if os.path.exists(os.path.join(self.folder, basename + '.sqt')):
            gz_flag = False
        elif os.path.exists(os.path.join(self.folder, basename + '.sqt.gz')):
            gz_flag = True
        else:
            print('in "write_sqt" fall through:')
            print('sqt_folder:', self.folder)
            print('basename:', basename)
            print(os.path.exists(os.path.join(self.folder, basename + '.sqt')))
            return
        if not scan_list:
            return

        # open the original SQT file and the output file
        if gz_flag:
            sqt_in = gzip.open(os.path.join(self.folder, basename + '.sqt.gz'))
        else:
            sqt_in = open(os.path.join(self.folder, basename + '.sqt'), 'rU')
        sqt_out = open(os.path.join(self.filtered_folder, basename + '_filtered.sqt'), 'w')

        # pass header lines and scan blocks that match scans in scan_list
        scan_set = set(scan_list)
        sqt_out.write('H\tComment\tFiltered SQT file created ' + time.asctime() + '\n')
        buff = []
        passes = False
        for line in sqt_in:
            line = line.strip()
            items = line.split('\t')
            if items[0] == 'H':
                sqt_out.write(line + '\n')
            if items[0] == 'S' and buff:
                for new in buff:
                    sqt_out.write(new + '\n')
                buff = []
                passes = False
            if items[0] == 'S' and (int(items[1]), int(items[2]), int(items[3])) in scan_set:
                passes = True
            if passes:
                buff.append(line)

        # need to process possible passing last scan
        if buff:
            for new in buff:
                sqt_out.write(new + '\n')
        for obj in self.write:
            print('......filtered SQT file written', file=obj)

        # close files and return
        sqt_in.close()
        sqt_out.close()
        return
        

    def write_ms2(self, basename, scan_list):
        """Writes a filtered MS2 file to 'path' folder
        """
        # see if there is anything to process
        ms2_out = os.path.join(self.filtered_folder, basename + '_filtered.ms2')
        if os.path.exists(os.path.join(self.folder, basename + '.ms2')):
            gz_flag = False
        elif os.path.exists(os.path.join(self.folder, basename + '.ms2.gz')):
            gz_flag = True
        else:
            return
        if not scan_list:
            return

        # open the original MS2 file and the output file
        if gz_flag:
            ms2_in = gzip.open(os.path.join(self.folder, basename + '.ms2.gz'))
        else:
            ms2_in = open(os.path.join(self.folder, basename + '.ms2'), 'r')
        ms2_out = open(os.path.join(self.filtered_folder, basename + '_filtered.ms2'), 'w')

        # pass header lines and passing scan blocks
        scan_set = set([(x[0], x[1]) for x in scan_list])    # don't want Z for MS2 files
        ms2_out.write('H\tComment\tFiltered MS2 file created ' + time.asctime() + '\n')
        buff = []
        passes = False
        for line in ms2_in:
            line = line.strip()
            items = line.split('\t')
            if items[0] == 'H':
                ms2_out.write(line + '\n')
            if items[0] == 'S' and buff:
                for new in buff:
                    ms2_out.write(new + '\n')
                buff = []
                passes = False
            if items[0] == 'S' and (int(items[1]), int(items[2])) in scan_set:
                passes = True
            if passes:
                buff.append(line)

        # have to worry about possible passing last scan
        if buff:
            for new in buff:
                ms2_out.write(new + '\n')
        for obj in self.write:
            print('......filtered MS2 file written', file=obj)

        # close files and return
        ms2_in.close()
        ms2_out.close()
        return

    # end class
               
class Threshold:
    """Hold low and high deltamass thresholds (Da).
    """
    def __init__(self):
        self.low = -5.0
        self.high = +5.0
