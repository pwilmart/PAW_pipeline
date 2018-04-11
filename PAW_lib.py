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
        
        # determine number of enzymatic termini, nct
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
        parseIPI: cleans up IPI entries
        parseNCBI: cleans up nr entries
        parseUniProt: cleans up Sprot/Trembl entries
        parseCONT: cleans up Contaminant entries
        reverseProtein: reverses sequences and modifies accession/descriptions
        molwtProtein: computes average MW of sequence
        frequencyProtein: returns aa composition dictionary
        seqlenProtein: returns aa sequence length
        findPeptide: finds location of peptide in protein sequence
        coverage: calculates coverage and aa counts from peptide list
        enzymaticDigest: theroetical enzymatic digest of protein sequence
        
    Written by Phil Wilmarth, OHSU, 2009, 2016.
    
    Updated for new Comet mod formats -PW 10/27/2017
    """
    def __init__(self):
        """Basic constructor, no parameters.
        """
        # bare bones __init__ function
        self.accession = 'blank'
        self.new_acc = 'blank'
        self.description = 'blank'
        self.new_desc = 'blank'
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
        self.new_acc = self.accession
        self.new_desc = self.description
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
        if self.new_desc == '':
            print('>'+self.new_acc, file=file_obj)
        else:
            print('>'+self.new_acc, self.new_desc, file=file_obj)
        
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
        rev_prot.new_acc = rev_prot.accession
        
        # change the desciptions, too.
        rev_prot.description = 'REVERSED'
        rev_prot.new_desc = 'REVERSED'
        
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
            next_protein.new_acc = next_protein.accession
            start = len(next_protein.accession)+2
            next_protein.description = self._last_line[start:]
            next_protein.new_desc = next_protein.description
        
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
        self.prompt = '\n(Auto Each Help List Pattern Reset Show Quit) ? '
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
            elif cmd[0] == 'Q':
                self.do_quit(arg)

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
...Quit: quit the command loop.""")
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
    
    def do_quit(self, arg=None):
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

################## end support functions for PAW pipeline use ##########################    
