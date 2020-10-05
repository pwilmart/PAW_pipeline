"""program "PAW_results.py"
Produces a protein results table and a peptide results table.

This is the protein inference step in the PAW analysis pipeline for Comet searches
of extremely large data sets.

Written by Phil Wilmarth, 2007-2017, OHSU.

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
# added tracking and reporting of peptide and protein FDRs (1/10/11 -PW)
# removed lots of legacy code (SEQUEST and DTASelect support) -PW 5/11/2016
# made tables easier to read into R -PW 2/13/2019

"""To do:
Keep dta_name variable??? Yes, dta_name is a handy key...

ntt determination should be from routine in PAW_lib
ACTUALLY SHOULD BE FROM NTT ATTRIBUTES. That needs some work to do -
may need to use tuples for dictionary values (currently counts, need to add ntt)

Make same set grouping and subset removal into functions

Make report table generations into three functions
"""
VERSION = 'v1.0.9'

# set protein ID criteria here
minimum_peptide_per_protein = 2     # min. number distinct peptides/protein/sample2
# set peptide requirements here
minimum_ntt_per_peptide = 2        # how many ntt for distinct peptides?


# seldom changed parameters:
# turns on or off calculating sum of MS2 fragment ion intensities
calc_ms2_int = False
max_num_peaks = 50
# turn on reporting of all peptide copies (True) or just best scoring spectrum (False)
full_peptide_list = True    # needs to be True for TMT data processing
# True for parsimony filtering, Note: subset proteins will not have any unique peptides
remove_subsets = True
# others:
minimum_unique_per_protein = 0      # min. number UNIQUE peptides/protein/sample (this should usually be zero!)
multiple_charge_states_ok = False   # are different charge states distinct?
modifications_ok = True             # are modified peptides distinct?
allow_prot_nterm_acytl = False      # +42 on M1 or Res2 considered as NTT=2
# NOTE: all parameters above are global in scope

###########################################################################################
######### coding starts here and below ####################################################
###########################################################################################
# global imports
import os
import time
import gzip
import sys
import re
import copy
import glob
import fnmatch
import PAW_lib
import PAW_protein_grouper

# some configuration parameters (also global in scope)
decoy_string = 'REV'
default_location = 'F:\PSR_Core_Analysis'

if not os.path.exists(default_location):
    default_location = os.getcwd()


def is_peptide_valid(peptide):
    """Tests if peptides meet sufficient criteria.

    Usage boolean = is_peptide_valid(peptide)
        where "peptide" is a peptide sequence (with bounding amino acids).
    Returns True if peptide criteria met and False otherwise.
    """    
    # see if N-term acetylation candidate
    pep = peptide
    nterm = False
    if ']' in pep and allow_prot_nterm_acytl == True and (pep.startswith('M') or pep.startswith('-')):
        nterm = True
        pep = pep.replace(']', '') # remove nt mod symbol
    
    # get ntt value
    """Should see if ntt needs to be replaced"""
    NTT = ntt(pep.split('_')[0])        
    if nterm and pep.startswith('M'):
        NTT = NTT + 1
    
    # see if modifications are OK 
    if (len(peptide_mods(pep.split('_')[0])) > 0) and not modifications_ok:
        print('failed mod test')
        return False
    
    # see if ntt criteria is satisfied (THIS IS ONLY CORRECT FOR TRYPSIN)
    if NTT >= minimum_ntt_per_peptide:
        return True
    else:
        return False

def is_protein_valid(prot, distinct_dict, unique_dict):
    """Test proteins for sufficient per sample evidence.

    Usage boolean = is_protein_valid(prot, distinct_dict, unique_dict)
        where "prot" is a PAW_Protein instance, "distinct_dict" is a
        list of full peptide sequences (with bounding amino acids) plus
        appended charge state, and "unique_dict" is ...
        
    Returns True if protein criteria met and False otherwise.
    """   
    # remember that "distinct" keys have charge appended to sequence   
    valid_peptides = {}
    unique_peptides = {}
    for peptide in distinct_dict.keys():
        if not is_peptide_valid(peptide):
            continue
        if not multiple_charge_states_ok: # are two different charge states OK
            peptide = peptide[:-2]
        valid_peptides[peptide] = 1
    for hit in prot.tophit_list:
        unique_peptides[hit.Filename] = unique_dict.get(hit.Filename, 0)
    
    # count number of distinct, valid peptides
##    # this filters out one-hit wonders only
##    if ((sum(valid_peptides.values()) >= minimum_peptide_per_protein) and
##        (sum(unique_peptides.values()) >= minimum_unique_per_protein) and
##        sum(valid_peptides.values()) < 2):
    if ((sum(valid_peptides.values()) >= minimum_peptide_per_protein) and
        (sum(unique_peptides.values()) >= minimum_unique_per_protein)):
        prot.valid.append(True)
    else:
        prot.valid.append(False)
    return

class PAW_Protein:
    """Holds information about distinct protein groups, no fields per se.
    """
    def __init__(self):
        self.locus_list = []    # list of PAW_Loci classes for protein group (may be one or more)
        self.tophit_list = []   # list of PAW_TopHit classes for each MS/MS spectra
        self.supersets = []     # list of protein group supersets (if any)
        self.subsets = []       # list of proteins with peptide subsets contained in this set
        self.spc_list = []      # list of total spectral count and total MS2 intensity per sample
        self.valid = []         # list of valid protein flags per sample

    def peptides(self):
        """Returns a set of peptides mapped to protein.
        """
        peptides = {}
        
        # loop over DTA files mapped to protein
        for hit in self.tophit_list:
            
            # see if different charge states of same peptide are OK
            if multiple_charge_states_ok:
                peptide = hit.distinct
            else:
                peptide = hit.distinct[:-2]
            
            # peptide count is the dictionary value
            if peptides.get(peptide, False):
                peptides[peptide] = peptides[peptide] + 1
            else:
                peptides[peptide] = 1
        
        # dictionary keys are the peptide set elements
        return set(peptides.keys())

    def peptides_strict(self):
        """Returns a set of strictly defined peptides mapped to protein.
        """
        peptides = {}
        
        # loop over DTA files mapped to protein
        for hit in self.tophit_list:
            if not is_peptide_valid(hit.Sequence):
                continue
            
            # see if different charge states of same peptide are OK
            if multiple_charge_states_ok:
                peptide = hit.distinct.split('.')[1] + hit.distinct[-2:]
            else:
                peptide = hit.distinct.split('.')[1] 
            
            # take care of I/L ambiguity
            peptide = re.sub(r'[IL]', r'j', peptide)
            
            # peptide count is the dictionary value
            if peptide in peptides:
                peptides[peptide] = peptides[peptide] + 1
            else:
                peptides[peptide] = 1
        
        # dictionary keys are the peptide set elements
        return set(peptides.keys())

class PAW_Locus:
    """Holds protein (locus) information.
    """
    def __init__(self):
        self.ID = 'Acession_Number'
        self.DTACount = '2'
        self.SpectrumCount = '2'
        self.Coverage = '100.0'
        self.SeqLength = '250'
        self.MW = '25000'
        self.pI = '5.7'
        self.Description = 'A long text string'

    def txt_parse_locus(self, line_items, col_map):
        """Loads locus data structure from PAW txt file line.
        """
        self.ID = line_items[col_map['Loci']]
        self.DTACount = 1
        self.SpectrumCount = 1
        self.Coverage = 0.0
        self.SeqLength = 0
        self.MW = 0.0
        self.pI = 0.0
        self.Description = ''
        return

class PAW_TopHit:
    """Holds top hit peptide information.
    """
    def __init__(self):
        self.Filename = 'my_dta.100.100.2'
        self.XCorr = '2.500'
        self.DeltCN = '0.200'
        self.NewDisc = '2.500'
        self.PrecursorMass = '1500.0'
        self.CalculatedMass = '1500.0'
        self.TotalIntensity = '5000.0'
        self.SpRank = '1'
        self.FragmentIonPercentage = '50.0'
        self.CopyCount = '1'
        self.Sequence = 'K.VTDDFGHR.A'
        self.Unique = 'FALSE'
        # calculated fields:
        self.DBUnique = 'FALSE'
        self.distinct = 'K.VTDDFGHR.A_2'
        self.ntt = 2
        self.sample = 'sample'
        self.beg = 0
        self.end = 0

    def txt_parse_tophit(self, line_items, dta_name, col_map):
        """Loads a tophit data structure from PAW txt line.
        """
        self.Filename = dta_name
        self.XCorr = line_items[col_map['Xcorr']]
        self.DeltCN = line_items[col_map['deltaCN']]
        self.NewDisc = line_items[col_map['NewDisc']]
        self.PrecursorMass = line_items[col_map['expM']]
        self.CalculatedMass = line_items[col_map['theoM']]
        self.TotalIntensity = '0.0'
        self.SpRank = line_items[col_map['SpRank']]
        self.FragmentIonPercentage = '0.0'
        self.CopyCount = 1
        self.Sequence = line_items[col_map['Sequence']]
        self.Unique = 'FALSE'
        self.distinct = self.Sequence + '_' + self.Filename[-1]
        self.ntt = line_items[col_map['ntt']]
        return

class FDR_Counter:
    """Data structure of counters for various target and decoy matches.
    """
    def __init__(self):
        self.FR = ['F', 'R']   # target or decoy
        self.NTT = ['0', '1', '2']     # number of tryptic termini
        self.Z = ['1', '2', '3', '4']    # peptide charge states
        """need to handle new terminla mod definitions"""
        self.mods = ['*', '#', '@', '^', '~', '$', '[', ']']
        self.unmod = [ [ [0 for z in self.Z] for ntt in self.NTT] for fr in self.FR]  # unmodified peptide counters
        self.mod = [ [ [0 for z in self.Z] for ntt in self.NTT] for fr in self.FR]   # modified peptide counters
        self.distinct = {}  # dictionary for distinct decoy sequences (conditioned on valid peptide criteria)
        self.tot_target = 0
        self.tot_valid_target = 0
        self.tot_decoy = 0
        self.tot_valid_decoy = 0

    def increment(self, seq, z, ntt, fr):
        """Increments counters for z, ntt, fr class.
        """        
        # test for valid z, ntt, and fr values
        if z not in self.Z:
            print('...FDR_counter WARNING: %s out of charge state range' % (z,))
            return
        if ntt not in self.NTT:
            print('...FDR_counter WARNING: %s out of NTT range' % (ntt,))
            return
        if fr not in self.FR:
            print('...FDR_counter WARNING: %s is not "F" or "R"' % (fr,))
            return
        
        # check seq for mods and if it meets valid peptide criteria
        unmod = False
        mod = False
        num_mods = 0
        valid = is_peptide_valid(seq)
        temp = seq.split('.')   # remove bounding residues if present
        if len(temp) > 0:
            seq = temp[1]        
        num_mods = len(peptide_mods(seq))
        if num_mods == 0: # unmodified peptides
            unmod = True
        else:
            mod = True
        
        # increment the appropriate counters
        Z = int(z) - 1
        NTT = int(ntt)
        if fr == 'F':
            FR = 0
            self.tot_target += 1
            if valid:
                self.tot_valid_target += 1
        else:
            FR = 1
            self.tot_decoy += 1
            if valid:
                self.tot_valid_decoy += 1
                self.distinct[seq] = True
        if unmod:
            self.unmod[FR][NTT][Z] += 1
        else:
            self.mod[FR][NTT][Z] += 1
        return

    def report(self, write):
        """Prints out a report with the FDR information.
        """
        for obj in write:
            print('\n########### FDR REPORT ############', file=obj)
            print('\nunmodified peptides:', file=obj)
            for ntt in [2, 1, 0]:
                string = 'ntt=%s: ' % (ntt,)
                for z in [0, 1, 2, 3]:
                    try:
                        rate = 100 * float(self.unmod[1][ntt][z]) / float(self.unmod[0][ntt][z])
                    except ZeroDivisionError:
                        rate = 0.0
                    string += '%s(%s)%0.2f, ' % (self.unmod[0][ntt][z], self.unmod[1][ntt][z], rate)
                print(string, file=obj)
            print('\nmodified peptides:', file=obj)
            for ntt in [2, 1, 0]:
                string = 'ntt=%s: ' % (ntt,)
                for z in [0, 1, 2, 3]:
                    try:
                        rate = 100 * float(self.mod[1][ntt][z]) / float(self.mod[0][ntt][z])
                    except ZeroDivisionError:
                        rate = 0.0
                    string += '%s(%s)%0.2f, ' % (self.mod[0][ntt][z], self.mod[1][ntt][z], rate)
                print(string, file=obj)
            print('\ndistinct decoy sequences: %s' % (len(self.distinct),), file=obj)
            try:
                rate = 100 * float(self.tot_decoy) / float(self.tot_target)
            except ZeroDivisionError:
                rate = 0.0
            print('total matches: %s (%s) %0.2f' % (self.tot_target, self.tot_decoy, rate), file=obj)
            try:
                rate = 100 * float(self.tot_valid_decoy) / float(self.tot_valid_target)
            except ZeroDivisionError:
                rate = 0.0
            print('total valid matches: %s (%s) %0.2f' % (self.tot_valid_target, self.tot_valid_decoy, rate), file=obj)

    # end class FDR_Counter

"""Only need database name! Do not need a data container."""        
def get_database(folder):
    """Get FASTA database name from a Comet parameters file or an SQT file header.
    """
    database = ''
    if os.path.exists(os.path.join(folder, 'comet.params')): # get DB from params file
        for line in open(os.path.join(folder, 'comet.params')):
            line = line.strip()
            if line.startswith('database_name'):
                database = line.split('= ')[1]
    elif os.path.exists(os.path.join(folder, 'sequest.params')): # get DB from params file
        for line in open(os.path.join(folder, 'sequest.params')):
            line = line.strip()
            if line.startswith('first_database_name'):
                database = line.split('= ')[1]

    # try getting database from SQT files            
    else: 
        hbuff = []
        try:
            # get header lines
            if glob.glob('*.sqt')[0]:
                for line in open(glob.glob('*.sqt')[0]):
                    if line.startswith('H'):
                        hbuff.append(line.strip())
                    if line.startswith('S'):
                        break
            elif glob.glob('*.sqt.gz')[0]:
                for line in gzip.open(glob.glob('*.sqt.gz')[0]):
                    if line.startswith('H'):
                        hbuff.append(line.strip())
                    if line.startswith('S'):
                        break

            # extract DB from header
            for line in hbuff:
                if line.startswith('H\tDatabase'):
                    database = line.split('\t')[2]
                    
        except IndexError:
            print('...WARNING: no params file or SQT files')

        if not os.path.exists(database):    # browse to database if path not found
            ext_list = [('FASTA files', '*.fasta'), ('Zipped files', '*.gz'), ('All files', '*.*')]
            title = 'Please locate the FASTA file'
            database = PAW_lib.get_file(default_location, ext_list, title)

    return database

def load_results_from_txt_files(txt_file_list, txt_to_sample, write):
    """Loads results from PAW "txt" files.
    """
    # build a (partial) PAW_Protein list from the filtered txt files
    matches = {}
    proteins = []
    fdr = FDR_Counter()
    for txt_file in txt_file_list:

        # get the header line and make the column map
        col_map = {}
        txt_base_name = os.path.basename(txt_file)
        base_name = txt_base_name.replace('.txt.gz', '') # do the longer possible extension first
        base_name = base_name.replace('.txt', '') # in case the extension was not '.txt.gz'
        try:
            if txt_file.endswith('.gz'):
                contents = [x.strip() for x in gzip.open(txt_file).readlines()]
            else:
                contents = [x.strip() for x in open(txt_file,'r').readlines()]
            for obj in write:
                print('...processing', os.path.split(txt_file)[1], file=obj)
        except:
            for obj in write:
                print('...WARNING: TXT and SQT file mis-match', file=obj) # TXT list built from SQT list
            continue
        try:
            if contents[0].startswith('start\tend'): # test for PAW header line
                for i, item in enumerate(contents[0].split('\t')):
                    col_map[item] = i
            else:
                print('...WARNING: non-PAW TXT file:', base_name)
                continue
        except IndexError:
            print('...WARNING: empty TXT file:', base_name)
            continue
        
        for line in contents[1:]:
            temp = line.split('\t')
            fdr.increment(temp[col_map['Sequence']], temp[col_map['Z']], temp[col_map['ntt']], temp[col_map['ForR']]) # count target, decoy matches
            name_list = [base_name] + [str(int(x)) for x in temp[:3]]
            dta_name = '.'.join(name_list)
            prot = temp[col_map['Loci']]
            
            # if protein already in "matches" dictionary, add the DTA information
            if matches.get(prot, False):
                old_prot = matches[prot]
                new_tophit = PAW_TopHit()
                new_tophit.txt_parse_tophit(temp, dta_name, col_map)
                new_tophit.sample = txt_to_sample[txt_base_name]
                old_prot.tophit_list.append(copy.deepcopy(new_tophit))
                old_prot.locus_list[0].DTACount += 1
                old_prot.locus_list[0].SpectrumCount += 1
                matches[prot] = old_prot
            
            # if not, add a new protein match to "matches" dictionary
            else:
                new_prot = PAW_Protein()
                new_locus = PAW_Locus()
                new_locus.txt_parse_locus(temp, col_map)
                new_prot.locus_list.append(copy.deepcopy(new_locus))
                new_tophit = PAW_TopHit()
                new_tophit.txt_parse_tophit(temp, dta_name, col_map)
                new_tophit.sample = txt_to_sample[txt_base_name]
                new_prot.tophit_list.append(copy.deepcopy(new_tophit))
                matches[prot] = copy.deepcopy(new_prot)
    
    # print out the FDR information
    fdr.report(write)
    
    # need strict list of peptides for each protein: minimum ntt, allow modifications or not, I/L, etc.
    prot_pep_list = {}
    for prot in matches.keys():
        prot_pep_list[prot] = matches[prot].peptides_strict()

#############################################################
    """this should be a function. we do this more than once"""
    # make list of proteins sorted by (strict) number of peptides
    prot_list = []
    for prot in matches.keys():
        prot_list.append((len(matches[prot].peptides_strict()), prot))
    prot_list.sort()
##############################################################
    
    # skip proteins with too few peptides to speed things up
    try:
        skip = [x[0] for x in prot_list].index(minimum_peptide_per_protein)
    except ValueError:
        skip = 0
    for obj in write:
        print('\n################ PARSIMONY REPORT ###############', file=obj)
        print('\n   there were', len(matches), 'protein matches', file=obj)
        print('   there were', skip, 'proteins with too few potential peptide(s)', file=obj)
    
    # clean up "matches" by removing skipped proteins
    for x, acc in prot_list[:skip]: 
        del matches[acc]
    for obj in write:
        print('   new matches length is:', len(matches), file=obj)

##############################################################
    """this should be a function
    should have a protein sorting function with option for ascending or descending"""
    # find redundant proteins to group together
    redundant_group = 1
    redundants = {}
    for i in range(skip, len(prot_list)):
        for j in range(i+1, len(prot_list)):
            s1 = prot_pep_list[prot_list[i][1]]
            s2 = prot_pep_list[prot_list[j][1]]
            if s1.issubset(s2) and len(s1) == len(s2):
                if redundants.get(prot_list[i][1], False):
                    redundants[prot_list[j][1]] = redundants[prot_list[i][1]]
                else:
                    redundants[prot_list[i][1]] = redundant_group
                    redundants[prot_list[j][1]] = redundant_group
                    redundant_group += 1
    
    # collaspe the redundant groups
    to_group = {}
    for (prot, group) in redundants.items():
        if to_group.get(group, False):
            temp = to_group[group]
            temp.append(prot)
            temp.sort()
            temp.insert(0, temp.pop())
            to_group[group] = temp
        else:
            to_group[group] = [prot]
    summary = list(to_group.items())
    summary.sort()
    for group_number, group_list in summary:
        group_list.sort()   # this determines protein group order
        keeper = group_list[0]
        for redundant in group_list[1:]:
            temp = matches[redundant].locus_list[0]
            matches[keeper].locus_list.append(copy.deepcopy(temp))
            
            # print out any redundant sets that are not exactly identical
            keep_set = matches[keeper].peptides()
            redun_set = matches[redundant].peptides()
            if keep_set != redun_set:
                keep_keys = sorted(list(keep_set))
                redun_keys = sorted(list(redun_set))
                for obj in write:
                    print('\n   redundant "mis-match": %s and %s' %
                          (matches[keeper].locus_list[0].ID, matches[redundant].locus_list[0].ID), file=obj)
                for i, first in enumerate(keep_keys):
                    try:
                        for obj in write:
                            if first != redun_keys[i]:
                                print('      %s %s %s' % (i, first, redun_keys[i]), file=obj)
                    except IndexError:
                        for obj in write:
                            print('      %s %s' % (i, first), file=obj)
            #
            del matches[redundant]
        for obj in write:
            print('   (%s) redundant group: %s' % (group_number, group_list), file=obj)
    for obj in write:
        print('\n   now matches is this long:', len(matches), file=obj)
    
    # build the list of protein objects sorted by accession
    proteins_list = []
    for i in range(skip, len(prot_list)):
        try:
            proteins_list.append((prot_list[i][1], matches[prot_list[i][1]]))
        except KeyError:    # redundant proteins have been deleted from "matches" 
            pass
    proteins_list.sort()
    for prot in proteins_list:
        proteins.append(prot[1])
        
################################################
################################################
    """This should be a function"""    
    # identify any protein/groups with peptide sets that are subsets
    # need to rebuild sets keyed to index in "proteins"
    prot_pep_list = {}
    for i in range(len(proteins)):
        prot_pep_list[i] = proteins[i].peptides_strict()
    
    # make list of proteins sorted by number of peptides
    prot_list = []
    keys = list(prot_pep_list.keys())
    keys.sort()
    for prot in keys:
        prot_list.append((len(prot_pep_list[prot]), prot))
    prot_list.sort()

    # test for subsets
    for i in range(len(prot_list)):
        for j in range(i+1, len(prot_list)):
            s1 = prot_pep_list[prot_list[i][1]]
            s2 = prot_pep_list[prot_list[j][1]]
            if s1.issubset(s2):
                proteins[prot_list[i][1]].supersets.append(prot_list[j][1])
                proteins[prot_list[j][1]].subsets.append(prot_list[i][1])
    
    # replace list of indexes with list of accessions for each protein
    for prot in proteins:
        superset_temp = []
        for index in prot.supersets:
            temp = []
            for locus in proteins[index].locus_list:
                temp.append(locus.ID)
            if len(temp) != 0:
                superset_temp.append(temp)
        if len(superset_temp) != 0:
            prot.supersets = superset_temp
        #
        subset_temp = []
        for index in prot.subsets:
            temp = []
            for locus in proteins[index].locus_list:
                temp.append(locus.ID)
            if len(temp) != 0:
                subset_temp.append(temp)
        if len(subset_temp) != 0:
            prot.subsets = subset_temp
    
    # print information about subsets before removing them from "proteins"
    subset_list = []
    for obj in write:
        print(file=obj)
    for i, prot in enumerate(proteins):
        if len(prot.supersets) != 0:
            subset_list.append(i)
            subset = [x.ID for x in prot.locus_list]
            for obj in write:
                print('   %s- %s subset of %s' % (i, subset, prot.supersets), file=obj)
    for obj in write:
        print(file=obj)
    subset_list.reverse()   # have to delete from top to bottom
    if remove_subsets:
        for i in subset_list:
            del proteins[i]
    
##    # print information from the supersets data
##    for i, prot in enumerate(proteins):
##        if len(prot.subsets) != 0:
##            superset = [x.ID for x in prot.locus_list]
##            print('   %s- %s superset of %s' % (i, superset, prot.subsets))
##    print()
    #
    for obj in write:
        print('   length of proteins is: %s\n' % (len(proteins),), file=obj)
    return proteins
######################################################################

def reversed_hit(locus_list, decoy_string):
    """Checks if any proteins are reversed (decoy) entries.
    """
    rev = False
    for loci in locus_list:
        if decoy_string in loci.ID:
            rev = True
    return rev

def ntt(sequence):
    """Counts number of tryptic termini for peptides with bounding residues.
    """
    """This has already been computed in the TXT files. Need to use that instead! (otherwise need to know enzyme)"""
    ntt = 0
    parts = sequence.split('.')
    if len(parts) < 3:
        print('   ### WARNING: need bounding residues to determine ntt! ###')
    else:
        length, ntt = PAW_lib.amino_acid_count(sequence) # this is going to be trypsin!
    return ntt
    
def get_database_proteins(fasta_full_name, proteins, write):
    """Get FASTA database entries for the identified proteins.

    Usage: (DBProteins, DB_map, fasta_file_name, DB_total) = get_database_proteins(database, proteins, write),
        where "DBProteins" is the list of FASTA objects, "DB_map" is a dictionary of accessions to indices,
        "fasta_file_name" is the basename for the database, "DB_total" is the toal number
        of proteins in the FASTA database, "database" is the path to the FASTA database,
        "proteins" is the list of identified proteins, and "write" is for console and log file use.
    """
    # if database path can't be found, then browse to database
    if not os.path.exists(fasta_full_name):
        fasta_name = os.path.basename(fasta_full_name)
        # bugger in an alternative database location:
        if os.path.exists(os.path.join('E:\Carr_plasma\databases', fasta_name)):
            fasta_full_name = os.path.join('E:\Carr_plasma\databases', fasta_name)
        else:
            print('...select the FASTA file')
            extensions = [('FASTA files', '*.fasta')]
            title = 'Select the %s database' % (os.path.basename(database),)
            fasta_full_name = PAW_lib.get_file(default, extensions, title)
            if not fasta_full_name: sys.exit()     # cancel button repsonse
    fasta_name = os.path.basename(fasta_full_name)
    
    # make a dictionary of identified loci so we don't have to load all proteins
    prot_id_map = {}
    for prot in proteins:
        for loci in prot.locus_list:
            prot_id_map[loci.ID] = True
    
    # open FASTA reader and create a Protein instance.
    for obj in write:
        print('...reading:', fasta_name, file=obj)
    f = PAW_lib.FastaReader(fasta_full_name)
    p = PAW_lib.Protein()
    
    # start looping over FASTA entries and keep the ones we need
    DBProteins = []
    DB_map = {}
    DB_len = 0
    DB_total = 0
    while f.readNextProtein(p):  # lets get everything copied into the data structures
        DB_total += 1
        if p.accession in prot_id_map:
            DBProteins.append(copy.deepcopy(p))
            DB_map[p.accession] = DB_len
            DB_len += 1
    for obj in write:
        print('...closing database (%s entries)...' % (DB_total,), file=obj)
    
    # compute coverage, length, and MW (needed when reading PAW TXT files)
    for prot in proteins:
        for loci in prot.locus_list:
            seq_dict = {}
            for hit in prot.tophit_list:
                seq_dict[hit.Sequence] = True
            try:
                loci.Coverage = '%0.01f' % (DBProteins[DB_map[loci.ID]].calcCoverage(seq_dict.keys())[0],)
                loci.SeqLength = DBProteins[DB_map[loci.ID]].seqlenProtein()
                loci.MW = '%0.0f' % (DBProteins[DB_map[loci.ID]].molwtProtein(),)
            except KeyError: # phrog DB
                loci.Coverage = '%0.01f' % (DBProteins[DB_map[loci.ID.split()[0]]].calcCoverage(seq_dict.keys())[0],)
                loci.SeqLength = DBProteins[DB_map[loci.ID.split()[0]]].seqlenProtein()
                loci.MW = '%0.0f' % (DBProteins[DB_map[loci.ID.split()[0]]].molwtProtein(),)
    
    # loop over proteins and get description strings
    for prot in proteins:
        for loci in prot.locus_list:
            try:
                loci.Description = DBProteins[DB_map[loci.ID]].description
            except KeyError: # phrog DB
                loci.Description = DBProteins[DB_map[loci.ID.split()[0]]].description
    #
    for obj in write:
        print('...done updating accessions and descriptions...', file=obj)
    return (DBProteins, DB_map, fasta_name, DB_total)

def compute_ms2_intensity(proteins, all_ms2, folder, write):
    """Computes total MS2 intensity from DTA files.

    Usage: compute_ms2_intensity(proteins, all_ms2, folder)
        where "proteins" is list of PAW_Protein objects, "all_ms2"
        is a dictionary of all protein matches of each MS2 file (by protein
        list index and tophit_list index), and "folder" is the
        folder containing the filtered MS2 format files.
    """
    # Assumes that MS2 files are in the same location as the TXT files,
    # open each MS2 file, build the filename for each scan, and check
    # it against the 'all_ms2' list.  If so, calculate DTA fragment ion sum
    # and replace the field(s) in the corresponding (protein, tophit) data structure.
    # max_num_peaks is a global set near top of module.
    os.chdir(folder)
    ms2_list = glob.glob('*.ms2')
    if not ms2_list:
        ms2_list = glob.glob('*.ms2.gz')
    if not ms2_list:
        print('...WARNING: MS2 file(s) were not found!')
        return
    for ms2 in ms2_list:
        if ms2.endswith('.gz'):
            ms2_file = gzip.open(ms2)
        else:
            ms2_file = open(ms2, 'r')
        buff = []   # temporary buffer to hold data for one MS2 scan
        for line in ms2_file:
            if line.startswith('H'):    # skip header lines
                continue
            line = line.rstrip()
            buff.append(line)   # save lines in scan buffer
            if line.startswith('S\t') and len(buff) > 1: # process previous scan block
                file_name = []
                start, end = buff[0].split('\t')[1], buff[0].split('\t')[2]  # start, stop scan numbers
                start = str(int(start))
                end = str(int(end))
                for top in buff[:15]:       # parse enough top lines to get charge states
                    if top.startswith('Z'):   # get the charge state from the Z line(s)
                        z = top.split('\t')[1]   
                        dta_name = ms2[:-3] + start + '.' + end + '.' + z    # build original DTA filename
                        if dta_name in all_ms2:
                            file_name.append(dta_name)
                for f in file_name:     # file_name list might be empty most of the time
                    frag_int = []
                    for i in range(len(buff)): # get fragment ion intensities
                        if buff[i].startswith('S\t') or buff[i].startswith('Z\t') or buff[i].startswith('I\t'):
                            continue
                        frag_int.append(float(buff[i].split()[1]))
                    frag_int.sort()
                    frag_int.reverse()
                    int_sum = sum(frag_int[:max_num_peaks])     # sum up desired number of peaks
                    for i in range(len(all_ms2[f])):     # DTA may match to more than one protein
                        prot, hit = all_ms2[f][i]   # get protein, tophit indices for new TotalIntensity
                        proteins[prot].tophit_list[hit].TotalIntensity = int_sum
                # reset scan buffer
                buff = []
                buff.append(line)
        
        # process the last buffer
        if len(buff) > 1:
            file_name = []
            start, end = buff[0].split('\t')[1], buff[0].split('\t')[2]  # start, stop scan numbers
            for j in buff[:15]:     # parse enough top lines to get charge states
                if j.startswith('Z'):
                    z = j.split('\t')[1]    # get the charge state from the Z line(s)
                    dta_name = ms2[:-3] + start + '.' + end + '.' + z # build original DTA filename
                    if dta_name in all_ms2:
                        file_name.append(dta_name)                         
            for f in file_name:
                frag_int = []
                for i in range(len(buff)): # get fragment ion intensities
                    if buff[i].startswith('S\t') or buff[i].startswith('Z\t') or buff[i].startswith('I\t'):
                        continue
                    frag_int.append(float(buff[i].split()[1]))
                frag_int.sort()
                frag_int.reverse()
                int_sum = sum(frag_int[:max_num_peaks])     # sum up desired number of peaks
                for i in range(len(all_ms2[f])):    # DTA may match to more than one protein
                    prot, hit = all_ms2[f][i]   # get protein, tophit indices for new TotalIntensity
                    proteins[prot].tophit_list[hit].TotalIntensity = int_sum
    #
    for obj in write:
        print('...done calculating fragment ion intensity sums...', file=obj)
    #
    return

def peptide_mods(peptide):
    """Looks for modification symbols in peptides.  Returns list of mod symbols.

    THIS NEEDS TO BE CHANGED TO HANDLE NEW COMET MODS

    """
    
    # see if there are bounding residues
    temp = peptide.split('.')
    if len(temp) > 1:
        peptide = temp[1]
    
    # check for Comet modification symbols (also has older style n-term and c-term symbols)
    mod_list = []
    for char in peptide:
        if char in ['*', '#', '@', '^', '~', '$', '%', '!', '+', 'n', 'c', '[', ']']:
            mod_list.append(char)
    #
    return mod_list            

def get_filter_tag(accession, decoy_string, count):
    """Checks if accession number corresponds to a decoy or contaminant.
    """
    tag = ''
    if 'CONT|' in accession or 'CONT_' in accession:
        tag = 'contaminant'
    if decoy_string in accession:
        tag = 'reversed'
    if count > 0:
        tag = 'redundant'
    return tag

def calc_split_count(prot, proteins, all_ms2_sample, s):
    """Calculates total SpC after spliting shared peptide counts.
    Also normalizes counts (adds 0.15 [or 100.0] to all zero values).
    """
    other_loci = {}
    split_count_total = 0.0

    # get correct values if using intensities
    if calc_ms2_int:
        self_unique = float(prot.spc_list[s][3])
    else:
        self_unique = float(prot.spc_list[s][1])
        
    for hit in prot.tophit_list:   # loop over all MS2 for protein "prot"
        if calc_ms2_int:
            value = float(hit.TotalIntensity)
        else:
            value = 1.0 # what we are splitting (a count or MS2 intensity)
        
        other_unique = 0.0
        if hit.Filename in all_ms2_sample:
            if len(all_ms2_sample[hit.Filename]) > 1:
                
                # split counts based on relative unique counts per sample
                for (p, x) in all_ms2_sample[hit.Filename]:
                    if calc_ms2_int:
                        other_unique += float(proteins[p].spc_list[s][3])
                    else:
                        other_unique += proteins[p].spc_list[s][1]
                try:
                    split_count_total += value * (self_unique / other_unique)
                except ZeroDivisionError: # do equal splitting if zero unique total
                    split_count_total += value / float(len(all_ms2_sample[hit.Filename]))
            else:
                split_count_total += value
                
    if calc_ms2_int:
        try:
            if split_count_total < .0001:
                split_count_total = 0.0
            value = split_count_total
            if value < 0.0:
                value = 0.0
        except:
            print('   calc_split_count WARNING: calc_ms2_int try except failed')
            value = 0.0
    else:
        value = split_count_total
    
    return ('%0.03f' % (value,))
        
#######################################################################################
# Main program to compile protein results from groups of filtered SQT/TXT files.
# Groups should be a collection of related samples comprising a single proteome
# The result will be a single list of valid proteins
# 5/2/08 -PW
# 6/10/08 - Added full accession number lookup since SEQUEST truncates Acc. No.s -PW
# 7/5/08 - trap duplicate DTAfile entries in the DTASelect-filter.xml file
#           - added MS2 Intensity table in addition to the spectral count table. -PW
# 4/23/09 - Improved reports and added partitioning SpC by unique peptide counts. -PW
# 1/1/2010 - Direct processing of PAW TXT files added (DTASelect no longer needed). -PW
# 3/30/2010 - Checked Parsimony logic, added more reporting of subsets, rejected proteins. -PW
# 8/2014 - Added dynamic parsing of TXT files so new filtered files are handled.
#               - Removed some out-of-date code sections. -PW
# 5/11/2016 - Removing legacy support sections and cleaning up comment style -PW
# 8/29/2017 - Updated for Python 3 -PW
# 10/04/2020 - Added sample mapping from file support -PW
#######################################################################################

print('=====================================================')
print(' program PAW_results.py, %s, Phil Wilmarth, OHSU ' % VERSION)
print('=====================================================')

# let user browse to folder with filtered sqt files
decoy = decoy_string
default = default_location
if default == '' or not os.path.exists(default):
    default = os.getcwd()
print('\n...select folder of filtered files')
filtered_folder = PAW_lib.get_folder(default, 'Select folder of filtered files')
if not filtered_folder: sys.exit()  # cancel button response

# put protein and peptide files in a "results_files" folder
results_folder = os.path.join(os.path.dirname(filtered_folder), 'results_files')
if not os.path.exists(results_folder):
    os.mkdir(results_folder)

# set up the log file in results folder
log_obj = open(os.path.join(results_folder, 'PAW_results.log'), 'wt')
write = [None, log_obj]
print('\n>>> starting PAW_results (%s) at: %s' % (VERSION, time.ctime()), file=log_obj)

# get the list of PAW txt files
os.chdir(filtered_folder)
txt_list = glob.glob('*_filtered.txt')
if len(txt_list) == 0:
    txt_list = glob.glob('*_filtered.txt.gz')
if len(txt_list) == 0:
    for obj in write:
        print('...WARNING: no sqt or txt files found!', file=obj)
    sys.exit()
txt_file_list = [os.path.join(filtered_folder, x) for x in txt_list]

# interactive shell to determine sample name mappings
paw_loop = PAW_lib.PAWShell(filtered_folder)
sample_dict, txt_to_sample = paw_loop.cmdloop('\nPAW command line loop. Type Help to see commands.')
samples = sorted(sample_dict.keys())

# if we do not have a samples list, something went wrong
if len(samples) == 0 or samples[0] is None:
    print('\n...PAW_results ERROR: parsing did not result in any samples!')
    sys.exit()

# print list of samples to console (should be what user expects)    
for obj in write:
    print('...sample keys have been determined:', file=obj)
    for i, sample in enumerate(samples):
        print('   (%s) %s' % (i+1, sample), file=obj)
    print(file=obj)

# get results from PAW txt files
proteins = load_results_from_txt_files(txt_file_list, txt_to_sample, write)

# get FASTA database file name and read in the identified protein sequences
database = get_database(filtered_folder)
(DBProteins, DB_map, DB_name, DB_total) = get_database_proteins(database, proteins, write)

# lookup peptide positions for first proteins in each group (for peptide reports)
for prot in proteins:
    db_prot = DB_map[prot.locus_list[0].ID]
    for hit in prot.tophit_list:
        match_list = PAW_lib.find_peptide(hit.Sequence, DBProteins[db_prot])
        try:
            hit.beg = match_list[0][1]
            hit.end = match_list[0][2]
        except IndexError:
            hit.beg = 0
            hit.end = 0

# build map of protein, tophit coordinates for MS2 filenames and peptide sequences
# used in intensity calculations and for identifying shared peptides
all_ms2 = {}
all_seq = {}
all_ms2_sample = [{} for s in samples]
for i in range(len(proteins)):
    for j in range(len(proteins[i].tophit_list)):
        
        # populate the entire dataset ms2 dictionary
        if proteins[i].tophit_list[j].Filename in all_ms2:
            L = all_ms2[proteins[i].tophit_list[j].Filename]
            L.append((i,j))
            all_ms2[proteins[i].tophit_list[j].Filename] = L
        else:            
            all_ms2[proteins[i].tophit_list[j].Filename] = [(i,j)]
        
        # populate the entire dataset peptide sequence dictionary
        Seq = proteins[i].tophit_list[j].Sequence.split('.')[1]
        Seq = re.sub(r'[IL]', r'j', Seq)    # take care of I/L ambiguity
        if Seq in all_seq:
            L = all_seq[Seq]
            L.append((i,j))
            all_seq[Seq] = L
        else:            
            all_seq[Seq] = [(i,j)]
        
        # populate the per sample dictionaries
        try:
            s = samples.index(proteins[i].tophit_list[j].sample)
        except ValueError:
            s = -1
        if s != -1:
            if proteins[i].tophit_list[j].Filename in all_ms2_sample[s]:
                L = all_ms2_sample[s][proteins[i].tophit_list[j].Filename]
                L.append((i,j))
                all_ms2_sample[s][proteins[i].tophit_list[j].Filename] = L
            else:            
                all_ms2_sample[s][proteins[i].tophit_list[j].Filename] = [(i,j)]
for obj in write:
    print('...ms2 filename maps created...', file=obj)        

# compute fragment ion intensities if desired
if calc_ms2_int:
    compute_ms2_intensity(proteins, all_ms2, filtered_folder, write)

# keep count of several quantities (spectral count, MS2 intensities).
# Create n dictionaries to hold peptides for each of the n samples
peptide = []    # will hold list of all dta filenames for a given protein group
distinct = []   # will hold number of distinct peptides for a given protein group
unique = []     # will hold number of unique peptides for a given protein group
intensity = []  # will hold list of all MS2 intensities for a given protein group
unique_int = [] # will hold list of unique MS2 intensities for a given protein group
normalization = []  # will hold total number of MS2 filenames (or intensities) for each sample
contaminant = []    # count contaminant and decoy matches per sample
for s in samples:
    peptide.append({})
    distinct.append({})
    unique.append({})
    intensity.append({})
    unique_int.append({})
    normalization.append({}) # this dictionary is NOT reset for each protein
    contaminant.append({}) # this dictionary is NOT reset for each protein

# fill up the peptides dictionaries with dta filenames for spectral counts
#   Note: filenames may occur more than once if peptide is in a protein more than once
#         using dictionaries makes sure each spectrum is counted once
# fill up the distinct dictionaries with Sequences for unique peptide counts
for p in range(len(proteins)):
    for s in range(len(samples)):   # reset dictionaries for each protein group
        peptide[s], distinct[s], unique[s] = {}, {}, {}
        intensity[s], unique_int[s] = {}, {}
        proteins[p].spc_list.append([])
    
    # loop over all top hits, assign spectra to respective samples, add to dictionaries
    for hit in proteins[p].tophit_list:
        try:
            s = samples.index(hit.sample)
        except ValueError:
            s = -1
        if s != -1:
            peptide[s][hit.Filename] = 1
            distinct[s][hit.distinct] = 1
            intensity[s][hit.Filename] = float(hit.TotalIntensity)
            if calc_ms2_int:
                value = float(hit.TotalIntensity)
            else:
                value = 1
            normalization[s][hit.Filename] = value
            if (proteins[p].locus_list[0].ID.startswith('CONT_') or
                proteins[p].locus_list[0].ID.startswith('CONT|') or
                proteins[p].locus_list[0].ID.startswith(decoy)):
                contaminant[s][hit.Filename] = value
            
            # I think this would be the same if per sample ms2 dictionaries were used
            if len(all_ms2[hit.Filename]) == 1:
                i, j = all_ms2[hit.Filename][0]
                proteins[i].tophit_list[j].Unique = 'TRUE'
                unique[s][hit.Filename] = 1
                unique_int[s][hit.Filename] = float(hit.TotalIntensity)
            else:
                same_locus = {}
                for i, j in all_ms2[hit.Filename]:
                    same_locus[i] = None
                if len(same_locus) == 1:
                    for i, j in all_ms2[hit.Filename]:
                        proteins[i].tophit_list[j].Unique = 'TRUE'
                    unique[s][hit.Filename] = 1
                    unique_int[s][hit.Filename] = float(hit.TotalIntensity)
    
    # need to compute total spectral (intensity) count for each protein
    for s in range(len(samples)):   
        spc_tot, uniq_tot = 0, 0
        int_tot, uniq_int_tot = 0.0, 0.0
        spc_tot = sum(peptide[s].values())     # sum of spectral counts
        uniq_tot = sum(unique[s].values())      # sum of unique spectra
        if calc_ms2_int:
            int_tot = sum(intensity[s].values())
            uniq_int_tot = sum(unique_int[s].values())                              
        proteins[p].spc_list[s] = [spc_tot, uniq_tot, int_tot, uniq_int_tot]
        
        # test for per sample protein evidence, criteria set in function
        is_protein_valid(proteins[p], distinct[s], unique[s])
        
for obj in write:
    print('...proteins and peptides sorted and counted...', file=obj)

# print out the valid protein results tables.  This is kind of messy because
# we have to make sure proteins not meeting criteria per sample are not included.
# We have multiple samples (usually), each of which has multiple protein, each
# protein has multiple peptides, each peptide may have been sequenced more than
# once, we have an option for either spectral counts or intensities, etc.

# Print status and details to console and log file
if calc_ms2_int:
    string = '...MS2 intensity-weighted protein summary:'
else:
    string = '...spectral count protein summary:'
for obj in write:
    print('\nRESULTS FOLDER: %s' % results_folder, file=obj)
    print(string, time.asctime(time.localtime()), file=obj)
    print('...writing protein_summary_9.txt...', file=obj)
protein_filename = os.path.join(results_folder, 'protein_summary_9.txt')
protein_results = open(protein_filename, 'w')

# first figure out valid protein spectral count numbers
valid_protein_count, valid_reversed = 0, 0
reject_protein_count, reject_reversed = 0, 0

table_rows = []
for prot in proteins:
    valid = False
    counts, uniq_counts, split_counts, items = [], [], [], []
    
    # see if protein is valid and get the list of spectral counts for each of the samples
    for s in range(len(samples)):
        valid = valid or prot.valid[s]
    if valid:
        for s in range(len(samples)):
            split_count = calc_split_count(prot, proteins, all_ms2_sample[s], s)                                           
            split_counts.append(split_count)                                                 
            if calc_ms2_int:
                counts.append('%0.03f' % (prot.spc_list[s][2],))
                uniq_counts.append('%0.03f' % (prot.spc_list[s][3],))
                counts_total = sum([prot.spc_list[s][2] for s in range(len(samples))])
                unique_total = sum([prot.spc_list[s][3] for s in range(len(samples))])
                try:
                    unique_frac = unique_total / counts_total
                except ZeroDivisionError:
                    print('...WARNING: total counts were zero for %s' % (prot.locus_list[0].ID,))
                    unique_frac = 0.0
                count_unique = '%0.03f\t%0.03f\t%0.03f' % (counts_total, unique_total, unique_frac)
            else:
                counts.append(str(prot.spc_list[s][0]))
                uniq_counts.append(str(prot.spc_list[s][1]))
                counts_total = sum([prot.spc_list[s][0] for s in range(len(samples))])
                unique_total = sum([prot.spc_list[s][1] for s in range(len(samples))])
                unique_frac = float(unique_total) / float(counts_total)
                count_unique = '%s\t%s\t%0.03f' % (counts_total, unique_total, unique_frac)
        counts = counts + uniq_counts + split_counts
        
        # add other loci for shared peptides, if any
        others = {}
        other_loci = []
        for hit in prot.tophit_list:
            Seq = hit.Sequence.split('.')[1]
            Seq = re.sub(r'[IL]', r'j', Seq) # take care of I/L masking
            for (p, x) in all_seq[Seq]:
                others[proteins[p].locus_list[0].ID] = None
        other = [x for x in others.keys()]
        other.sort()
        if len(other) > 1:
            other_loci = [', '.join(other)]
        else:
            other_loci = [' ']
        
        # generate the output line for valid proteins
        valid_protein_count += 1
        if reversed_hit(prot.locus_list, decoy): valid_reversed += 1
        k = 0
        for loci in prot.locus_list:
            
            # build up the output line's items
            filter_tag = get_filter_tag(loci.ID, decoy, k)
            items = [str(valid_protein_count+0.0001*k), '1', loci.ID,
                     filter_tag, loci.Coverage, loci.SeqLength,
                     loci.MW, loci.Description, count_unique] + counts + other_loci
            string = '\t'.join([str(x) if str(x) else ' ' for x in items])  # add tabs and print          
            table_rows.append(string)
            k += 1
    else:
        # keep count of rejected proteins
        reject_protein_count += 1
        if reversed_hit(prot.locus_list, decoy): reject_reversed += 1

# build the header line and print
header = ['ProtGroup', 'Counter', 'Accession', 'Filter', 'Coverage', 'SeqLength', 'MW', 'Description',
          'CountsTot', 'UniqueTot', 'UniqFrac']
##preheader = (['', '=subtotal(109,B3:B3)'] + (len(header)-2)*[''] + len(samples)*['Total'] +
##             len(samples)*['Unique'] + len(samples)*['Corrected'] + [''])
pre_lines = 4
print('PAW protein_summary_9 results file', file=protein_results)
print('written:', time.ctime(), file=protein_results)
print(file=protein_results)
formula = '=subtotal(109,B%d:B%d)' % (pre_lines+2, pre_lines+1+len(table_rows))
preheader = ([' ', formula] + (len(header)-2)*[' '] + 3*len(samples)*[' '] + ['x'])
print('\t'.join(preheader), file=protein_results)
header = (header + ['Total_'+x for x in samples] + ['Unique_'+x for x in samples] +
          ['Corrected_'+x for x in samples] + ['OtherLoci'])
print('\t'.join([x if x else ' ' for x in header]), file=protein_results)

for row in table_rows:
    print(row, file=protein_results)

"""Things after the main table are not easy to parse with R"""
# print total spectral counts per sample (normalization factors)
##if calc_ms2_int:
##    string1 = 'Per Sample Total MS2 Intensity:'
##    string2 = 'Per Sample Contaminant MS2 Intensity:'
##else:
##    string1 = 'Per Sample Total Spectral Count:'
##    string2 = 'Per Sample Total Contaminant Count:'
##for s in range(len(samples)):
##    string1 += '\t' + str(sum(normalization[s].values()))
##    string2 += '\t' + str(sum(contaminant[s].values()))
##
##print('Protein summary 9 companion file:', file=protein_companion)
##print('written on:', time.ctime(), file=protein_companion)
##print(file=protein_companion)
##print(string1, file=protein_companion)
##print(string2, file=protein_companion)

# finally, print summary stats to log file and console
for obj in write:
    print('...valid proteins: %s(%s) [rejected proteins: %s(%s)]' %
          (valid_protein_count, valid_reversed, reject_protein_count, reject_reversed), file=obj)
    print('...database: %s (%s entries)' % (DB_name, DB_total), file=obj)
    print('...minimum_ntt_per_peptide: %s' % (minimum_ntt_per_peptide,), file=obj)
    print('...minimum_peptide_per_protein: %s' % (minimum_peptide_per_protein,), file=obj)
    print('...minimum_unique_per_protein: %s' % (minimum_unique_per_protein,), file=obj)
    # less common parameters
    if calc_ms2_int:
        print('...calc_ms2_int flag: %s' % (calc_ms2_int,), file=obj)
        print('......max_num_peaks: %s' % (max_num_peaks,), file=obj)
    if not full_peptide_list:
        print('...full_peptide_list flag: %s' % (full_peptide_list,), file=obj)
    if not remove_subsets:
        print('...remove_subsets flag: %s' % (remove_subsets,), file=obj)
    if multiple_charge_states_ok:
        print('...multiple_charge_states_ok flag: %s' % (multiple_charge_states_ok,), file=obj)
    if not modifications_ok:
        print('...modifications_ok flag: %s' % (modifications_ok,), file=obj)
    if allow_prot_nterm_acytl:
        print('...allow_prot_nterm_acytl flag: %s' % (allow_prot_nterm_acytl,), file=obj)
    print(file=obj)
#
protein_results.close()

# now do the peptide summary report by counts across samples: build and print header line
for obj in write:
    print('...writing peptide_summary_9.txt...', file=obj)
peptide_results = open(os.path.join(results_folder, 'peptide_summary_9.txt'), 'w')
if calc_ms2_int:
    string = 'Program "peptide_summary_9.py" report with MS2 intensites performed on'
else:
    string = 'Program "peptide_summary_9.py" report with spectral counts performed on'
print(string, time.asctime(time.localtime()), file=peptide_results)

header = ['\n\n\nProtGroup', 'Accession', 'Sequence', 'Begin', 'End', 'Unique', \
          'NTT', 'Z', 'TotCount'] + samples
header.append('OtherLoci')
print('\t'.join([x if x else ' ' for x in header]), file=peptide_results)

valid_count = 0
for p in proteins:
    valid = False
    peptide_matches = {}
    
    # count valid proteins to keep peptides in sync with protein report
    for s in range(len(samples)):
        valid = valid or p.valid[s]
    if valid:
        valid_count += 1    # aka protein group number
        
        # for each peptide, print the spectral counts across the samples
        # find all of the proteins associated with non-unique peptides
        for hit in p.tophit_list:
            charge = hit.Filename[-1]
            other_loci = ' '
            if len(all_ms2[hit.Filename]) == 1:
                unique_ms2 = 'TRUE'
            else:
                unique_ms2 = 'FALSE'
                
                # trap repeated motifs in the same protein
                same_locus = {}
                for i, j in all_ms2[hit.Filename]:
                    same_locus[str(i)] = None
                if len(same_locus) == 1:
                    unique_ms = 'TRUE'
                others = {}
                Seq = hit.Sequence.split('.')[1]
                Seq = re.sub(r'[IL]', r'j', Seq) # take care of I/L masking
                for (prot, x) in all_seq[Seq]:
                    others[proteins[prot].locus_list[0].ID] = None
                other = [x for x in others.keys()]
                other.sort()
                if len(other) > 1:
                    other_loci = ', '.join(other)
                else:                    
                    # if we get here, we have different sequences with same score
                    unique_ms2 = 'TRUE'
            
            # get all values for current peptide
            sort_key = hit.beg + (hit.end/100000.)
            line = [sort_key, hit.distinct, valid_count, p.locus_list[0].ID,
                    hit.Sequence, hit.beg, hit.end, unique_ms2, hit.ntt, charge, other_loci]
            
            # figure out which sample count to increment
            try:
                sample_index = samples.index(hit.sample)
            except ValueError:
                sample_index = -1
            if sample_index != -1:                
                # if peptide has already been seen keep the best scoring peptide
                if hit.distinct in peptide_matches:
                    line = peptide_matches[hit.distinct]
                    if calc_ms2_int:
                        line[11][sample_index] += float(hit.TotalIntensity)
                    else:
                        line[11][sample_index] += 1
                    peptide_matches[hit.distinct] = line
                
                # peptide not yet seen, so save values and initialize counters
                else:
                    count = [0 for s in samples]
                    if calc_ms2_int:
                        count[sample_index] += float(hit.TotalIntensity)
                    else:
                        count[sample_index] += 1
                    line.append(count)
                    peptide_matches[hit.distinct] = line        
        
        # get a list of the information for each peptide, sort, and print.
        print_list = list(peptide_matches.values())
        print_list.sort()
        for i, line in enumerate(print_list):
            print_line = line[2:10]
##            count = line[11]    # nested list of counts
            print_line.append(sum(line[11]))
            print_line += line[11]
            print_line.append(line[10])
            string = '\t'.join([str(x) if str(x) else ' ' for x in print_line])
            print(string, file=peptide_results)

# close peptide results file
peptide_results.close()

# now do the peptide report for each sample with scoring details, masses, etc.
for s in range(len(samples)):
    for obj in write:
        print('...writing %s...' % (samples[s] + '_peptide_results_9.txt',), file=obj)
    peptide_results = open(os.path.join(results_folder, samples[s] + '_peptide_results_9.txt'), 'w')
    print('Peptide report for sample "%s" performed on %s' %
          (samples[s], time.asctime(time.localtime())), file=peptide_results)
    header = ['\n\n\nProtGroup', 'Accession', 'Sequence', 'Unique', \
              'TotCount', 'NTT', 'XCorr', 'DeltaCN', 'SpRank', 'NewDisc', 'Z', 'Delta_Mass', \
              'Exp_Mass', 'Calc_Mass', 'DTA_filename']
    print('\t'.join([x if x else ' ' for x in header]), file=peptide_results)
    #
    valid_count = 0
    for p in proteins:
        
        # count valid proteins to keep peptides in sync with protein report
        valid = False
        for sample in range(len(samples)):
            valid = valid or p.valid[sample]
        if valid:
            valid_count += 1
            print_list = []
            best_hit = {}
            for hit in p.tophit_list:
                try: 
                    if s == samples.index(hit.sample):
                        sort_key = hit.beg + (hit.end/100000.)
                        diff = float(hit.PrecursorMass) - float(hit.CalculatedMass) # use SEQUEST or Comet calculated masses
                        charge = hit.Filename[-1]
                        if len(all_ms2[hit.Filename]) == 1:
                            unique_ms2 = 'TRUE'
                        else:
                            unique_ms2 = 'FALSE'
                        
                        # get all values for current peptide
                        if calc_ms2_int:
                            count = hit.TotalIntensity
                        else:
                            count = 1
                        line = [sort_key, hit.distinct, valid_count, p.locus_list[0].ID,
                                hit.Sequence, unique_ms2, count, hit.ntt, hit.XCorr, hit.DeltCN,
                                hit.SpRank, hit.NewDisc, charge, diff, hit.PrecursorMass,
                                hit.CalculatedMass, hit.Filename]
                        
                        # if full report
                        if full_peptide_list:
                            print_list.append(line)
                        
                        # or keep peptide with the best score
                        else:
                            if hit.distinct in best_hit:
                                if calc_ms2_int:
                                    best_hit[hit.distinct][6] = str( float(best_hit[hit.distinct][6]) + float(line[6]))
                                else:
                                    best_hit[hit.distinct][6] += 1
                                if float(line[8]) >= float(best_hit[hit.distinct][8]):
                                    line[6] = best_hit[hit.distinct][6]
                                    best_hit[hit.distinct] = line
                            
                            # peptide not yet seen, so save values
                            else:
                                best_hit[hit.distinct] = line
                except ValueError:
                    pass
            
            # get a list of the best score values for each peptide, sort, and print.
            # sort by increasing amino acid residue position in sequence currently
            if not full_peptide_list:
                print_list = list(best_hit.values())
            print_list.sort()
            for line in print_list:
                print('\t'.join([str(x) if x else ' ' for x in line[2:]]), file=peptide_results)
    peptide_results.close()

# run protein grouping
PAW_protein_grouper.main(protein_filename)

# echo final stats and ending message, close log file
for obj in write:
    print('\nValid proteins: %s(%s) [rejected proteins: %s(%s)]' %
          (valid_protein_count, valid_reversed, reject_protein_count,
           reject_reversed), file=obj)
    print('Ending PAW_results (%s) at: %s' % (VERSION, time.ctime()), file=obj)
    
log_obj.close()

# write a tables description file
tables_filename = os.path.join(results_folder, 'PAW_table_descriptions_9.txt')
with open(tables_filename, 'w') as tables_obj:
    for line in PAW_lib.column_keys():
        print(line, file=tables_obj)
    
# Fini (Yeah!!!)

