"""program "make_PAW_TXT_from_PD2.x.py"
Extracts information from PD 2.x PSM export files and makes PAW pipeline compatible text files.
The text fles also contain the reporter ion quantities.
Written by Phil Wilmarth, OHSU, fall 2017.

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

import os
import sys
import time
import copy
import re
from io import StringIO
from itertools import groupby
import pandas as pd
import numpy as np
import PAW_lib

VERSION = 'v1.0.1'

# globals
SEPARATOR = '\t'  # column separator character
QVALUE = 0.05     # q-Value cutoff
PPM = 20.0        # plus/minus width for PPM deltamass window
INTERF = 50.0     # maximum % interference in precursor isolation

"""Discussion:
It is tempting to do some filterieng of reporter ion information in this script.
That is better to do after protein inference. We want as many confident PSMs
as possible to build the best list of proteins. Removing PSMs before protein
inference will decrease the number of identified proteins.

v1.0.1: Removed options to test reporter ion intensities or input zero replacements. -PW 10/17/2017
"""


class Fasta(object):
    """Class for FASTA related things: sequences, indexes."""
    def __init__(self, db_filename, write):
        self.filename = db_filename     # FASTA file path
        self.write = write              # file objects for logging prints
        self.proteins = []              # list of Protein objects
        self.prot_index = {}            # key: accession, value: protein list index
        self.prot_desc = {}             # key: accession, value: description
        self.seen = {}                  # used in accession/description repair

        self.read_proteins()
        self.make_protein_index()
        self.make_prot_desc()

    def read_proteins(self):
        """Read in the protein sequences from FASTA file."""
        for obj in self.write:
            print('Reading proteins from %s' % os.path.basename(self.filename), file=obj)
        for prot in fasta_iter(self.filename):
            self.proteins.append(prot)
        for obj in self.write:
            print('%s proteins sucessfully read in' % len(self.proteins), file=obj)

    def make_protein_index(self):
        """Indexes proteins by accession parts."""
        self.skip = set(['sp', 'tr', 'gi', 'ref', ' ', ''])
        for i, p in enumerate(self.proteins):
            self.prot_index[p.accession] = i    # index by full accession string
            accs = p.accession.split('|')
            for acc in accs:
                if acc in self.skip:
                    continue
                self.prot_index[acc] = i    # also index by parts of compound accessions

    def make_prot_desc(self):
        """Makes accesssion to description dictionary."""
        for p in self.proteins:
            self.prot_desc[p.accession] = p.description

class Counter(object):
    """Generic container for counter attributes."""
    def __init__(self):
        return

def fasta_iter(fasta_name):
    """Yields Protein objects from fasta files.
    Adapted from "https://www.biostars.org/p/710/" post.
    See also: https://drj11.wordpress.com/2010/02/22/python-getting-fasta-with-itertools-groupby/
    """
    with open(fasta_name, mode='rt') as fasta_handle:
        # skip x[0] (boolean something from groupby) and keep the alternating header, sequences
        # more on groupby here: https://docs.python.org/3/library/itertools.html
        fasta_iter = (x[1] for x in groupby(fasta_handle, lambda line: line[0] == '>'))
        for header in fasta_iter:
            # make a container
            p = PAW_lib.Protein()
            
            # process the header line
            p.header = header.__next__()[1:].rstrip()
            p.accession = p.header.split()[0]
            p.new_acc = p.accession 
            p.description = p.header[len(p.accession)+1:]
            p.new_desc = p.description
            
            # join all sequence lines
            p.sequence = "".join(s.strip() for s in fasta_iter.__next__())
            yield p
        
class PSM(object):
    """Container for PSM information."""
    def __init__(self, row_series, channels, separator='\t'):
        """Creates container for each row, and parses row elements."""
        self.channels = channels    # set of reporter ion labels
        self.separator = separator  # separator character
        
        # row_series is a pandas series object
        self.confidence = row_series['Confidence']
        self.psm_ambiguity = row_series['PSM Ambiguity']
        self.annotated_sequence = row_series['Annotated Sequence']
        try:
            self.sequence = row_series['Sequence']
        except KeyError:
            self.sequence = row_series['Annotated Sequence'].upper()
        self.modifications = row_series['Modifications']
        try:
            self.number_groups = self._int(row_series['Number of Protein Groups'])
        except KeyError:
            self.number_groups = self._int(row_series['# Protein Groups'])
        try:
            self.number_proteins = self._int(row_series['Number of Proteins'])
        except KeyError:
            self.number_proteins = self._int(row_series['# Proteins'])
        self.master_accessions = row_series['Master Protein Accessions']
        try:
            self.master_descriptions = row_series['Master Protein Descriptions']
        except KeyError:
            self.master_descriptions = None
        self.accessions = row_series['Protein Accessions']            
        try:
            self.descriptions = row_series['Protein Descriptions']
        except KeyError:
            self.descriptions = None
        try:
            self.missed_cleavages = self._int(row_series['Number of Missed Cleavages'])
        except KeyError:
            self.missed_cleavages = self._int(row_series['# Missed Cleavages'])        
        self.charge = self._int(row_series['Charge'])
        try:
            self.original_charge = self._int(row_series['Original Precursor Charge'])
        except KeyError:
            self.original_charge = None
        try:
            self.delta_score = self._float(row_series['Delta Score'])
        except KeyError:
            self.delta_score = self._float(row_series['DeltaScore'])
        if self.delta_score == 1.0:
            self.delta_score = 0.0
        try:
            self.delta_cn = self._float(row_series['Delta Cn'])
        except KeyError:
            self.delta_cn = self._float(row_series['DeltaCn'])
        self.rank = self._int(row_series['Rank'])
        self.search_engine_rank = self._int(row_series['Search Engine Rank'])
        try:
            self.concatenated_rank = self._int(row_series['Concatenated Rank'])
        except KeyError:
            self.concatenated_rank = None
        try:
            self.moverz = self._float(row_series['mz in Da'])
        except KeyError:
            self.moverz = self._float(row_series['m/z [Da]'])
        try:
            self.mhplus = self._float(row_series['MHplus in Da'])
        except KeyError:
            self.mhplus = self._float(row_series['MH+ [Da]'])
        try:
            self.theo_mhplus = self._float(row_series['Theo MHplus in Da'])
        except KeyError:
            self.theo_mhplus = None
        try:
            self.deltamassppm = self._float(row_series['Delta M in ppm'])
        except KeyError:
            self.deltamassppm = self._float(row_series['DeltaM [ppm]'])       
        try:
            self.deltamassda = self._float(row_series['Delta mz in Da'])
        except KeyError:
            self.deltamassda = self._float(row_series['Deltam/z [Da]'])
        self.ions_matched = row_series['Ions Matched']
        try:
            self.matched_ions = row_series['Matched Ions']
        except KeyError:
            self.matched_ions = None
        try:
            self.total_ions = row_series['Total Ions']
        except KeyError:
            self.total_ions = None
        try:
            self.intensity = self._float(row_series['Intensity'])
        except KeyError:
            self.intensity = None            
        self.activation = row_series['Activation Type']
        self.ms_order = row_series['MS Order']
        try:
            self.interference = self._int(row_series['Isolation Interference in Percent'])
        except KeyError:
            self.interference = self._int(row_series['Isolation Interference [%]'])
        try:
            self.ave_reporter_SN = row_series['Average Reporter SN']
        except KeyError:
            self.ave_reporter_SN = row_series['Average Reporter S/N']
        try:
            self.inject_time = self._int(row_series['Ion Inject Time in ms'])
        except KeyError:
            self.inject_time = self._int(row_series['Ion Inject Time [ms]'])
        try:
            self.rt = self._float(row_series['RT in min'])
        except KeyError:
            self.rt = self._float(row_series['RT [min]'])
        self.first = self._int(row_series['First Scan'])
        try:
            self.last = self._int(row_series['Last Scan'])
        except KeyError:
            self.last = self.first
        try:
            self.master_Scans = row_series['Master Scans']
        except KeyError:
            self.master_Scans = None
        self.spectrum_file = row_series['Spectrum File'].replace('.raw', '')
        try:
            self.file_id = row_series['File ID']
        except KeyError:
            self.file_id = None

        # get the reporter ion values
        self.intensities = [0.0 for x in self.channels]
        for i, channel in enumerate(self.channels):
            self.intensities[i]  =  self._float(row_series[channel], 0.0)
            
        self.quant_info = row_series['Quan Info']
        try:
            self.quant_info = row_series['Peptide Quan Usage']
        except KeyError:
            self.quant_info = None
        try:
            self.peptides_matched = row_series['Peptides Matched']
        except KeyError:
            self.peptides_matched = None
        self.xcorr = self._float(row_series['XCorr'])
        try:
            self.contam = row_series['Contaminant']
        except KeyError:
            self.contam = None        
        self.qval = self._float(row_series['Percolator q-Value'])
        self.pep = self._float(row_series['Percolator PEP'])
        try:
            self.SVM_score = self._float(row_series['Percolator SVMScore'])
        except KeyError:
            self.SVM_score = None       

        # some computed/added attributes        
        self.peptide_length = len(self.sequence)    # compute peptide length
        self.total = sum(self.intensities)          # sum of reporter ion values 
        self.meets_all_criteria = False             # Flag for final PSM status
        self.valid_acc = True                       # False if accession could not be found in protein database                      
        self.start = 0                              # locations of peptide in protein sequence
        self.match = []                             # peptide match tuples to portein sequence
        return

    def _float(self, string, default=0.0):
        """Converts string to a float, set to "default" if ValueError (missing)."""
        try:
            val = float(string)
        except ValueError:
            val = default
        if str(val) == 'nan':
            val = default
        return val

    def _int(self, string):
        """Converts string to an integer, zero if ValueError."""
        try:
            val = int(string)
        except ValueError:
            val = 0
        return val

    def make_header(self):
        """Makes a header line for PSM data."""
        header_line = (['Counter', 'Confidence', 'Sequence', 'PSM Ambiguity', 'Protein Descriptions',
                        'Number of Protein Groups', 'Protein Group Accessions', 'Modifications', 'Activation Type',
                        'Delta Score', 'Rank']
                       + self.channels +
                       ['Total Int', 'Quan Info', 'Percolator q-Value', 'Percolator PEP',
                        'Peptides Matched', 'XCorr', 'Number of Missed Cleavages',
                        'Isolation Inteference in percent', 'Ion Inject Time in ms', 'Intensity', 'Charge', 'mz in Da',
                        'MHPlus in Da', 'Delta Mass in Da', 'Delta Mass in ppm', 'RT in min', 'First Scan',
                        'Last Scan', 'Spectrum File', 'MeetsAllCriteria', 'New Sequence', 'New Modifications',
                        'Peptide Length', 'PeptideMatchInfo', 'SEQUEST', 'length', 'ntt', 'ForR'])
        return '\t'.join(header_line)

    def make_data(self):
        """Makes a data line for PSM data
        """
        data_list = ([1, self.confidence, self.sequence, self.psm_ambiguity, self.descriptions,
                      self.number_groups, self.accessions, self.modifications, self.activation,
                      self.delta_score, self.rank] +
                     self.intensities +
                     [self.total, self.quant_info, self.qval,
                      self.pep, self.peptides_matched, self.xcorr, self.missed_cleavages,
                      self.interference, self.inject_time, self.intensity, self.charge, self.moverz,
                      self.mhplus, self.deltamassda, self.deltamassppm, self.rt, self.first, self.last,
                      self.spectrum_file, self.meets_all_criteria, self.new_sequence, self.new_modifications,
                      self.peptide_length, self.match, self.SEQUEST, self.length, self.ntt, self.ForR])
        return '\t'.join([str(x) for x in data_list])
    
# end classes

def analyze_headers(headers):
    """Gets the list of TMT channel headers."""
    # look for reporter ion channels
    channels = [c for c in headers if c.startswith('Abundance')]
    if not channels:
        channels = [c for c in headers if (c.startswith('12') or c.startswith('13'))]
    if not channels:
        print('...WARNING: no reporter ion channels were found')
    return channels
            
##def test_trimmed_average_intensity(intensities, intensity):
##    """Finds average intensity of channels, excluding the top and bottom values.
##    """
##    int_vector = sorted(intensities)
##
##    # compute the trimmed average
##    average = sum(int_vector[1:-1])/float(len(int_vector[1:-1]))
##        
##    # test threshold
##    if average >= intensity:
##        return True
##    else:
##        return False
    
def amino_acid_count(sequence_string, enzyme='Tryp', return_base_pep=False):
    """Counts amino acids in peptides.  Returns (length, ntt) tuple.

    Usage: (length, ntt) = amino_acid_count(sequence_string),
        where "sequence_string" is a peptide sequence with bounding residues,
        "enzyme" is a string for the specific protease used,
        "length" is the returned number of amino acids, and
        "ntt" is the number of tryptic termini.

    Written by Phil Wilmarth, OHSU, 2008.
    """    
    # find the string between the bounding periods '.'
    parts = len(sequence_string.split('.'))
    if parts == 3: # we have bounding residues (supports PTMs as embedded deltamass values)
        start = sequence_string.index('.') + 1   # start is after first period
        temp = sequence_string[::-1] # reverse string
        end = temp.index('.')+1     # find first period in reversed string
        end = len(sequence_string) - end     # end is past the last period
    elif parts == 1:
        start = 0
        end = len(sequence_string)
    else:
        print('...amino_acid_count WARNING: incorrect number of periods -', sequence_string)
        if return_base_pep:
            return (0, 0, "")
        else:
            return(0, 0)
    sequence = sequence_string[start:end]
    
    # remove any special characters from the sequence string
    newseq = ''
    for c in sequence:
        if c.isalpha() and c.isupper():
            newseq += c
    
    # get the prefix and suffix amino acids
    prefix = sequence_string[start-2:start-1]
    if (prefix == "") or (start == 0):
        prefix = "X"    # no bounding residue info so unknown AA
    cterm = newseq[-1]  # last amino acid in sequence
    nterm = newseq[0]   # first amino acid in sequence
    suffix = sequence_string[end+1:end+2]
    if suffix == "":
        suffix = "X"    # no bounding residue info so unknown AA
    
    # determine number of enzymatic termini, ntt
    ntt = 0
    if enzyme.upper() == 'TRYP':  # cleaves at c-side of K, R
        if (prefix in 'KR-'):
            ntt += 1
        if (cterm in 'KR') or (suffix == '-'):
            ntt += 1
    elif enzyme.upper() == 'GLUC':  # cleaves at c-side of D, E
        if prefix in 'DE-':
            ntt += 1
        if (cterm in 'DE') or (suffix == '-'):
            ntt += 1
    elif enzyme.upper() == 'ASPN': # cleaves at n-side of D
        if (prefix == '-') or (nterm == 'D'):
            ntt += 1
        if suffix in 'D-':
            ntt += 1
    else:
        print('...amino_acid_count WARNING: unknown enzyme -', enzyme)
    
    # return length, number of tryptic termini, and (optional) base peptide sequence
    if return_base_pep:
        return (len(newseq), ntt, newseq)
    else:
        return (len(newseq), ntt)    
    # end

def repair_accessions(psm, db):
    """Repairs protein group accessions from PD 1.4."""
    # check the accession or accession list
    try:
        if psm.accessions.endswith('"'):        # not sure if this is needed
            print('accessions had ending quote')
            psm.accessions = psm.accessions[:-1]
    except AttributeError:
        # pandas puts in "nan" for missing data, which is a float not a string
        psm.accessions = ''

    # lookup any missing protein accession group lists
    temp_list = [x.strip() for x in psm.accessions.split(';')]   # clean list of accessions
    if not [x for x in temp_list if x]:   # can have peptides with no accessions in PD 2.1
        matches = []
        for prot in db.proteins:
            m = prot.findPeptide(psm.sequence.upper())
            if m:
                matches.append(prot.accession)
        psm.accessions = '; '.join(matches) # gather up the accessions into a string
        if not psm.accessions:
            print('...WARNING: peptide not found in DB:', psm.sequence)
        
    if psm.accessions in db.seen:   # we will see psm accession strings multiple times
        psm.descriptions = db.seen[psm.accessions][1]
        psm.accessions = db.seen[psm.accessions][0]
    else:
        accessions = []
        descriptions = []
        # PD separates with ';' but ';' can be in description strings, too.
        # this should get actual accessions and they should be in db.prot_index
        for item in [x.split()[0] for x in psm.accessions.split(';')]:
            if (item not in db.skip) and (item in db.prot_index):
                pidx = db.prot_index[item]
                accessions.append(db.proteins[pidx].accession)
                descriptions.append(db.proteins[pidx].description.replace(';', ':')) # replace any semicolons in descs
            else:
                print('...WARNING: %s not in DB index' % item)
        # make repaired accession and description strings
        acc_str = '; '.join(accessions)
        desc_str = '; '.join(descriptions)
        # save in db.seen dictionary
        db.seen[psm.accessions] = (acc_str, desc_str)
        # replace strings
        psm.descriptions = desc_str
        psm.accessions = acc_str                    

def parse_psm_lines(table, db, max_qvalue=0.05, max_ppm=20.0, max_interfere=50.0, separator='\t'):
    """Parses PSM exports and adds information to PSM objects.
    Returns a list of all PSMs passing q-value, deltamass, and intensity cutoffs.
    """
    # define some fixed ranges of valid PSM attributes
    min_charge = 2
    max_charge = 4
    min_length = 7
    max_length = 45

    # initialize counters, etc.
    c = Counter()
    c.qval_good = 0         # PSMs with q-value less than cutoff (good)
    c.qval_bad = 0          # PSMs with higher q-values (bad)
    c.in_ppm = 0            # PSMs with accurate deltamass inside PPM window (good)
    c.out_ppm = 0           # PSMs with accurate deltamass outside PPM window (bad)
    c.interfere_low = 0     # low % co-isolation inteference (good)
    c.interfere_high = 0    # higher % co-isolation interence (bad)
    c.total = 0             # total lines in PSM file
    c.top = 0               # number of PSMs that were top ranked hits
    c.valid = 0             # overall number of PSMs surviving tests
    c.reject = 0            # overall number of PSMs failing the tests
    c.misc_reject = 0       # numbers rejected by length, charge
    psm_list = []

    # get the channel list
    channels = analyze_headers(table.columns)
    
    start = False    
    for i, row in table.iterrows():
        
        c.total += 1
        psm = PSM(row, channels, separator)
        if not psm:
            c.total += -1
            c.reject += -1
            continue

        # PSM tables have extra lines for non-top hits to get deltaCN values, skip those
        if psm.rank == 1:
            c.top += 1
        else:
            continue

        # test various criteria
        if psm.qval <= max_qvalue:
            c.qval_good += 1
        else:
            c.qval_bad += 1
            c.reject += 1
            continue
            
        if abs(psm.deltamassppm) <= max_ppm:
            c.in_ppm += 1
        else:
            c.out_ppm += 1
            c.reject += 1
            continue

        if psm.interference <= max_interfere:     # smaller is better
            c.interfere_low += 1
        else:
            c.interfere_high += 1
            psm.intensities = [0.0 for x in psm.channels]

        if min_charge <= psm.charge <= max_charge:
            pass
        else:
            c.reject += 1
            c.misc_reject += 1
            continue

        if min_length <= psm.peptide_length <= max_length:
            pass
        else:
            c.reject += 1
            c.misc_reject += 1
            continue

        psm.meets_all_criteria = True
        c.valid += 1
        
        # add valid PSMs to list (expand matches to multiple proteins)
        repair_accessions(psm, db)  # repair possible messed up accessions and descriptions
        for acc in psm.accessions.split(';'):
            acc = acc.strip()
            new_psm = copy.deepcopy(psm)
            new_psm.accessions = acc
            if acc == '':
                new_psm.descriptions = ''
            else:
                new_psm.descriptions = db.prot_desc[acc]
            psm_list.append(new_psm)
                    
    return c, psm_list

def analyze_modifications(psm_list):
    """gets freqeuncies of amino acids in sequences, gets frequencies of modifications
    """
    aa_freq = {}
    all_mods = {}
    
    aa_freq['N-Term'] = len(psm_list)
    aa_freq['C-Term'] = len(psm_list)
    for psm in psm_list:
        amino_acid_frequency(psm.sequence, aa_freq)
        mods = parse_mods(psm.modifications)
        update_dictionary(all_mods, mods)

    return aa_freq, all_mods

def amino_acid_frequency(pepstring, aa_freq):
    """Counts amino acids frequencies of peptide sequences from PD,
    sequence and frequency dictionary are passed as arguments.
    """
    pepstring = pepstring.upper()
    for aa in pepstring:
        if aa in aa_freq:
            aa_freq[aa] +=1
        else:
            aa_freq[aa] = 1
    return

def parse_mods(modstring):
    """Parses PD modification descriptions to get
    modification types and count of affected residues.
    
    modifications returned as dictionary of modtypes and dictionary of residues and counts.
    empty mod descriptions should return empty structures.
    """
    mods = {}

    # split modification description string
    modlist = modstring.split(';')
    for mod in modlist:
        mod = mod.strip()   # get rid of whitespace
        mod = mod.replace(')(', '_')    # protein terminal mods look like: "N-Term(Prot)(Acetyl)"
        temp = mod[:-1].split('(') # separate part inside (); mod strings end in ")"
        residue = ''.join([c for c in temp[0] if not c.isdigit()]) # ignore positions
        modtype = temp[1]

        # mods is a dictionary of modtype where each modtype value is a count dictionary {residue: count}
        if modtype in mods:
            if residue in mods[modtype]:
                mods[modtype][residue] += 1
            else:
                mods[modtype][residue] = 1
        else:
            mods[modtype] = {residue: 1}
    return mods

def update_dictionary(big, little):
    """Big and little are dictionaries with dictionaries of count values.
    Big probably has more keys than little.
    """
    for key in little:
        if key in big:
            for secondkey in little[key]:
                if secondkey in big[key]:
                    big[key][secondkey] += little[key][secondkey]
                else:
                    big[key][secondkey] = little[key][secondkey]
        else:
            big[key] = little[key]
    return

def fixed_or_variable(all_mods, aa_freq):
    """Assigns special symbols to variable mods. None assigned to fixed mods."""
    # define some symbols and such
    reg_symbols = ['*', '#', '@', '^', '~', '$', '%', '!', '+']     # as of last Comet iteration (first 6 match SEQUEST)
    nt_symbols = [']', ')', '}']     # old style terminal mods for SEQUEST (SEQUEST only had one)
    ct_symbols = ['[', '(', '{']     # old style terminal mods for SEQUEST (SEQUEST only had one)
    labels = ['regular mods', 'N-term mods', 'C-term mods']
    symbols = [reg_symbols, nt_symbols, ct_symbols]

    mod_type = {}   # double duty, denotes which mods are static or variable and the variable mod symbol
    global_variable_freq = {}   # frequency counts of the different variable modifications
    # if mod frequency is the same as total aa frequency, mod was static
    for mod in all_mods:
        for residue in all_mods[mod]:
            if all_mods[mod][residue] == aa_freq[residue]:
                mod_type[mod] = None
            else:
                mod_type[mod] = residue
                
    # get overall variable mod frequencies                
    for mod in all_mods:
        if mod_type[mod]:
            global_variable_freq[mod] = (mod_type[mod], sum(all_mods[mod].values()))

    # separate the mods by location
    reg_list = [x for x in global_variable_freq.items() if '-Term' not in x[1][0]]
    nt_list = [x for x in global_variable_freq.items() if x[1][0] == 'N-Term']
    ct_list = [x for x in global_variable_freq.items() if x[1][0] == 'C-Term']

    # iterate over the lists and set the mod symbols by decreasing frequency
    for i, mod_list in enumerate([reg_list, nt_list, ct_list]):
        variable_freq = sorted(mod_list, key=lambda x: x[1][1], reverse=True)
        for j, (mod, (residue, frequency)) in enumerate(variable_freq):
            try:
                mod_type[mod] = symbols[i][j]
            except IndexError:
                print('...WARNING: maximum number of %s exceeded (reusing last symbol: %s)' %
                      (labels[i], symbols[i][-1]))
                mod_type[mod] = symbols[i][-1]
    return mod_type

def get_variable_positions(psm, mod_type):
    """Parses PD modification descriptions to get modification positions and symbols."""
    modmask = {}

    # split modification description string
    modlist = psm.modifications.split(';')
    for mod in modlist:
        mod = mod.strip()   # get rid of whitespace
        mod = mod.replace(')(', '_')    # protein terminal mods look like: "N-Term(Prot)(Acetyl)"
        temp = mod[:-1].split('(')  # separate part inside ()
        position = ''.join([c for c in temp[0] if c.isdigit()]) # get positions
        try:
            position = int(position)
        except ValueError:
            if temp[0] == 'N-Term':
                position = 1
            elif temp[0] == 'C-Term':
                position = len(psm.sequence)
            else:
                position = 0
        modtype = temp[1]

        if mod_type[modtype]:   # just variable mods
            modmask[position-1] = mod_type[modtype]
    return modmask

def print_modification_report(all_mods, mod_type, write):
    """Prints a summary of variable and static modifications
    """
    variable = []
    for obj in write:
        print('\nVariable modifications:', file=obj)   # print mod type, symbol, affected residues
    for mod in mod_type:
        if mod_type[mod]:
            variable.append([mod, mod_type[mod], sorted(all_mods[mod].keys()), sum(all_mods[mod].values())])
    variable = sorted(variable, reverse=True, key=lambda x: x[3])   # order mods by decreasing frequency
    for mod in variable:
        for obj in write:
            print('  ', mod[0], mod[1], mod[2], file=obj)

    static = []
    for obj in write:
        print('Static modifications:', file=obj)   # print mod type and afected residues
    for mod in mod_type:
        if not mod_type[mod]:
            static.append([mod, sorted(all_mods[mod].keys()), sum(all_mods[mod].values())])
    static = sorted(static, reverse=True, key=lambda x: x[2])   # order mods by decreasing freqiuency
    for mod in static:
        for obj in write:
            print('  ', mod[0], mod[1], file=obj)
    for obj in write:
        print(file=obj)
    return

def fix_PTM_info(psm, mod_type):
    """Makes SEQUEST-style sequences and removes static mods from modifications strings."""
    new_seq = list(psm.sequence.upper())
    new_symbols = ['' for x in new_seq]
    modmask = get_variable_positions(psm, mod_type)
    for index in modmask:
        new_symbols[index] = modmask[index]
    psm.new_sequence = ''.join([j for i in zip(new_seq, new_symbols) for j in i])

    new_mod_list = [x.strip() for x in psm.modifications.split(';') if mod_type[x[:-1].replace(')(', '_').split('(')[1]]]
    psm.new_modifications = '; '.join(new_mod_list)

def make_protein_index(proteins):
    """Indexes proteins."""
    prot_index = {}
    skip = set(['sp', 'tr', 'gi', 'ref'])
    for i, p in enumerate(proteins):
        accs = p.accession.split('|')
        for acc in accs:
            if acc in skip:
                continue
            prot_index[acc] = i
    return prot_index        

def lookup_peptides(psm_list, proteins, prot_index, write):
    """Finds starting residue number for peptide sequence in protein sequence.
    """
    for psm in psm_list:
        try:
            # eventually add lookup of all accessions, just first to test
            acc = psm.accessions.split(';')[0].strip() # usual parsing
            try:
                psm.match = proteins[prot_index[acc]].findPeptide(psm.new_sequence)
                psm.start = psm.match[0][0]
            except KeyError:
                for obj in write:
                    print('...WARNING: acc not in prot_index', acc, file=obj)
                try:
                    acc = acc.split()[0]
                    psm.match = proteins[prot_index[acc]].findPeptide(psm.new_sequence) # phrog DB
                    psm.start = psm.match[0][0]
                except KeyError:
                    psm.valid_acc = False
        except IndexError:
            for obj in write:
                print('\n...peptide lookup issue:', file=obj)
                print('psm accessions:', psm.accessions, file=obj)
                print('acc:', acc, file=obj)
                print('index:', prot_index[acc], file=obj)
                print('full acc:', proteins[prot_index[acc]].accession, file=obj)
                print('peptide:', psm.new_sequence, file=obj)
                print('peptide in sequence?', psm.new_sequence in proteins[prot_index[acc]].sequence, '\n', file=obj)
            psm.start = 0
        
def make_group_index(psm_list):
    """Makes an index of PSM indexes that have the same grouper string
    """
    group_index = {}
    for i, psm in enumerate(psm_list):
        if psm.meets_all_criteria:
            if psm.grouper in group_index:
                group_index[psm.grouper].append(i)
            else:
                group_index[psm.grouper] = [i]
    return group_index
 
def add_SEQUEST_sequences(psm_list):
    for psm in psm_list:
        try:
            psm.SEQUEST = psm.match[0][3]
        except IndexError:
##            print('\nindex error:', psm.match)
##            psm._snoop()
            psm.SEQUEST = 'X.X.X'   # index error should mean a failed peptide lookup
        psm.length, psm.ntt = amino_acid_count(psm.SEQUEST)
        psm.ForR = 'F'
        

###########################################
######## main program starts here #########
###########################################

# get the PSM results PD export file information
default_location = r'F:\PSR_Core_Analysis'
if not os.path.exists(default_location):
    default_location = os.getcwd()
print('Select the PD 2.x PSM export file')
psm_filename = PAW_lib.get_file(default_location,
                                [('Text files', '*.txt'), ('CSV files', '*.csv'), ('All files', '*.*')],
                                'Select a PD2.x PSM export file')
if not psm_filename: sys.exit()     # exit if not file selected

# get the path so we can open a log file there
path, basename = os.path.split(psm_filename)
log_obj = open(os.path.join(path, 'make_PAW_TXT_from_PD2.X.log'), 'wt')
write = [None, log_obj]

for obj in write:
    print(file=obj)
    print('====================================================================', file=obj)
    print(' make_PAW_TXT_from_PD2.x.py, %s, written by Phil Wilmarth, OHSU ' % VERSION, file=obj)
    print('====================================================================', file=obj)
    print('Ran on:', time.ctime(), file=obj)
    print('QVALUE = %s, PPM (+/-) = %s, INTERF = %s%%' % (QVALUE, PPM, INTERF), file=obj)

# get the FASTA database file
print('Select the FASTA protein database file')
db_filename = PAW_lib.get_file(default_location,
                               [('Fasta files', '*.fasta'), ('All files', '*.*')],
                               'Select a FASTA database')
if not db_filename: sys.exit()

for obj in write:   # echo selected files
    print('PSM export file:', basename, file=obj)
    print('FASTA file:', os.path.basename(db_filename), file=obj)

# read in the FASTA file information and make some indexes
db = Fasta(db_filename, write)

# load PD 2.x export into a pandas data frame
if psm_filename.endswith('.csv'):
    table = pd.read_csv(psm_filename, low_memory=False)
elif psm_filename.endswith('.txt'):
    table = pd.read_table(psm_filename, low_memory=False)
else:
    print('...Did you select the correct PSM file?')
    sys.exit()

# get the filtered psm information from the PSM export file
counts, psm_list = parse_psm_lines(table, db, QVALUE, PPM, INTERF, SEPARATOR)

# analyze PTMs: fixed or variable? what residues? what names?
aa_freq, all_mods = analyze_modifications(psm_list)

# analyze mods and print summary
mod_type = fixed_or_variable(all_mods, aa_freq)
print_modification_report(all_mods, mod_type, write)

# add alternatively formatted sequence string (SEQUEST style), reformat modifications string (remove static mods)
for psm in psm_list:
    fix_PTM_info(psm, mod_type)

# find peptide starting residue numbers
lookup_peptides(psm_list, db.proteins, db.prot_index, write)

# add true SEQUEST peptide string as attribute
add_SEQUEST_sequences(psm_list)

# when everything is done, sort and write to new files
# open results file for all PSMs, print header lines
psmout = StringIO()
##dump = open(os.path.join(path, 'psm_dump.txt'), 'wt')
print(psm_list[0].make_header(), file=psmout)
##print(psm_list[0].make_header(), file=dump)

# print PSM data to file sorted by decreasing total intensity
psm_list = sorted(psm_list, key=lambda x: x.total, reverse=True)
for psm in psm_list:
    if psm.valid_acc:
        print(psm.make_data(), file=psmout)
##        print(psm.make_data(), file=dump)
##dump.close()

# read the table into a pandas data frame
psmout.seek(0)
df = pd.read_table(psmout, low_memory=False)

# close StrinIO "file"
psmout.close()
where = os.path.dirname(psm_filename)

# copy the essentials to TXT data frame
txt_table = df[['First Scan', 'Last Scan', 'Charge', 'MHPlus in Da']].copy()
txt_table.columns = ['start', 'end', 'Z', 'expM']
txt_table['SpRank'] = df['Rank']
txt_table['theoM'] = df['MHPlus in Da'] + df['Delta Mass in Da']
txt_table['deltaCN'] = df['Delta Score']
txt_table['Xcorr'] = df['XCorr']
txt_table['Sequence'] = df['SEQUEST']
txt_table['Loci'] = df['Protein Group Accessions']
txt_table['NewDeltaCN'] = df['Delta Score']
txt_table['NewDisc'] = df['Percolator q-Value']
txt_table['ntt'] = df['ntt']
txt_table['ForR'] = df['ForR']

# add the reporter ion peak heights
cols = psm_list[0].channels
new_cols = [x.replace('Abundance ', '') for x in cols]
new_cols = ['height_' + x for x in cols]
txt_table[new_cols] = df[cols]

# add LC run names for filtering
txt_table['LCRun'] = df['Spectrum File']

# write a file for each LC run (also a dummy SQT file for each)
txt_table.sort_values(by='start', inplace=True)
folder = os.path.join(os.path.dirname(psm_filename), 'filtered_files')
if not os.path.exists(folder):
    os.mkdir(folder)
for base_name in txt_table['LCRun'].unique():
    mask = txt_table['LCRun'] == base_name
    txt_table[mask].drop('LCRun', axis=1).to_csv(os.path.join(folder, base_name+'_filtered.txt'), sep='\t', index=False)
    fsqt = open(os.path.join(os.path.join(folder, base_name + '_filtered.sqt')), 'w')
    print('H\tDatabase\t%s' % db_filename, file=fsqt)
    fsqt.close()

# print the parsing statistics
for obj in write:
    print('There were', counts.total, 'Total rows in PSM table export', file=obj)
    print('There were', counts.top, 'Top ranked PSMs', file=obj)
    print('There were', counts.valid, 'PSMs that passed all criteria', file=obj)
    print('There were %d total PSMs that were rejected (%0.2f%% rejected)' %
          (counts.reject, 100.0 * counts.reject / (counts.valid + counts.reject)), file=obj)
    print('There were %d PSMs that were rejected by len, Z (%0.2f%% rejected)' %
          (counts.misc_reject, 100.0 * counts.misc_reject / (counts.valid + counts.misc_reject)), file=obj)
    print('There were', counts.qval_good, 'PSMs passing q-value cutoff', file=obj)
    print('There were %d PSMs that failed q-value cutoff (%0.2f%% rejected)' %
          (counts.qval_bad, 100.0 * counts.qval_bad / (counts.qval_good + counts.qval_bad)), file=obj)
    print('There were', counts.in_ppm, 'PSMs passing deltamass cutoff', file=obj)
    print('There were %d PSMs that failed deltamass cutoff (%0.2f%% rejected)' %
          (counts.out_ppm, 100.0 * counts.out_ppm / (counts.in_ppm + counts.out_ppm)), file=obj)
    print('There were', counts.interfere_low, 'PSMs with low interfering signal in MS1', file=obj)
    print('There were %d PSMs with too much interference in MS1 (%0.2f%% of PSM intensities zeroed)' %
          (counts.interfere_high, 100.0 * counts.interfere_high / (counts.interfere_low + counts.interfere_high)), file=obj)

log_obj.close()

# end



