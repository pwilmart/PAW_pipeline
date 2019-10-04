""" program "make_PAW_TXT_from_PD_export.py"
Makes PAW TXT files from PAW PSM dump files created during the
processing of TMT export files from Proteome Discoverer 1.4.

added median reporter ion intensity cutoff, spring 2014 -PW
added deltamass filter, spring 2015 -PW
double checked things and added some additional comments, Nov. 2015, -PW
added reformatted PSM data dump, Nov. 2015, -PW
parses extra data from phospho searches, Feb. 2016, -PW
added protein description strings, July 2016, -PW
changed for Python 3, Aug. 2017, -PW

written by Phil Wilmarth 8/2016.
    
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
# import common modules here
import os
import sys
import time
import copy
from io import StringIO

import numpy as np
import pandas as pd

import PAW_lib

"""Since this processing comes before the protein inference in the PAW pipeline,
we need to retain as many PSMs as possible. We do not have access to the parent ion
interference in the PAW TXT files so the reporter ion intensities for PSMs that
fail the test are zeroed out rather than the PSM dropped. Any reporter ion minimum
intensity tesing moves to the add_TMT_intensities step. The input for an reporter
ion channels with zero should be done at the protein level after all of the PSMs
have been summed. -PW 10/18/2017
"""
# globals
SEPARATOR = '\t'  # column separator character
QVALUE = 0.05     # q-Value cutoff
PPM = 20.0        # plus/minus width for PPM deltamass window
INTERF = 50.0    # maximum % interference in precursor isolation

VERSION = 'v1.0.2'  # October 2017 -PW


class Fasta(object):
    """Class for FASTA related things: sequences, indexes."""
    def __init__(self, db_filename, write):
        self.filename = db_filename     # FASTA file path
        self.write = write              # list of file objects for writing
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
        f = PAW_lib.FastaReader(self.filename)
        p = PAW_lib.Protein()
        while f.readNextProtein(p, False):
            prot = copy.deepcopy(p)
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
    
class PSM(object):
    """Container for PSM information. """
    def __init__(self, line, header_map, channels, separator='\t'):
        """Parse line and populate fields
        """
        items = line.split(separator)
        self.channels = channels

        # note that PD exports fields enclosed in double quotes
        self.confidence_level = items[header_map['Confidence Level']][1:-1]
        self.sequence = items[header_map['Sequence']][1:-1]
        self.peptide_length = len(self.sequence)    # compute peptide length
        self.ambiguity = items[header_map['PSM Ambiguity']][1:-1]
        self.descriptions = items[header_map['Protein Descriptions']][1:-1]
        self.number_proteins = self._int(items[header_map['# Proteins']][1:-1])
        self.number_groups = self._int(items[header_map['# Protein Groups']][1:-1])
        self.accessions = items[header_map['Protein Group Accessions']][1:-1]
        self.modifications = items[header_map['Modifications']][1:-1]
        self.number_phospho = self.modifications.count('Phospho')   # add count of phospho sites
        self.activation = items[header_map['Activation Type']][1:-1]
        self.delta_score = self._float(items[header_map['DeltaScore']][1:-1])
        if self.delta_score == 1.0:
            self.delta_score = 0.0
        self.deltacn = self._float(items[header_map['DeltaCn']][1:-1])
        self.rank = self._int(items[header_map['Rank']][1:-1])
        self.se_rank = self._int(items[header_map['Search Engine Rank']][1:-1])
        self.intensities = [0.0 for x in self.channels]
        for i, channel in enumerate(self.channels):
            self.intensities[i]  =  self._float(items[header_map[channel]][1:-1], 0.0)
        self.total = sum(self.intensities)
        self.quant_info = items[header_map['Quan Info']][1:-1]
        self.quant_usage = items[header_map['Quan Usage']][1:-1]
        try:
            self.binomial_score = self._int(items[header_map['phosphoRS Binomial Peptide Score']][1:-1])
        except:
            self.binomial_score = 0
        self.qval = self._float(items[header_map['q-Value']][1:-1])
        self.pep = self._float(items[header_map['PEP']][1:-1])
        self.decoy_candidates = self._int(items[header_map['Decoy Peptides Matched']][1:-1])
        self.target_candidates = self._int(items[header_map['Peptides Matched']][1:-1])
        self.xcorr = self._float(items[header_map['XCorr']][1:-1])
        self.missed = self._int(items[header_map['# Missed Cleavages']][1:-1])
        self.interference = self._int(items[header_map['Isolation Interference [%]']][1:-1])
        self.inject_time = self._int(items[header_map['Ion Inject Time [ms]']][1:-1])
        self.intensity = self._float(items[header_map['Intensity']][1:-1])
        self.charge = self._int(items[header_map['Charge']][1:-1])
        self.moverz = self._float(items[header_map['m/z [Da]']][1:-1])
        self.mhplus = self._float(items[header_map['MH+ [Da]']][1:-1])
        self.grouper = self.sequence.upper() + '_' + str(int(round(self.mhplus, 0)))    # key for grouping similar PSMs
        self.deltamassda = self._float(items[header_map['Delta Mass [Da]']][1:-1])
        self.deltamassppm = self._float(items[header_map['Delta Mass [PPM]']][1:-1])
        self.rt = self._float(items[header_map['RT [min]']][1:-1])
        self.first = self._int(items[header_map['First Scan']][1:-1])
        self.last = self._int(items[header_map['Last Scan']][1:-1])
        self.ms_order = items[header_map['MS Order']][1:-1]
        self.ions_matched = items[header_map['Ions Matched']][1:-1]
        self.matched_ions = items[header_map['Matched Ions']][1:-1]
        self.total_ions = items[header_map['Total Ions']][1:-1]
        self.spectrum_file = items[header_map['Spectrum File']][1:-1].replace('.raw', '')
        self.meets_all_criteria = False
        self.valid_acc = True
        self.match = []
        return

    def _float(self, string, default=0.0):
        """Converts string to a float, set to "default" if ValueError (missing)
        """
        try:
            val = float(string)
        except ValueError:
            val = default
        return val

    def _int(self, string):
        """Converts string to an integer, zero if ValueError
        """
        try:
            val = int(string)
        except ValueError:
            val = 0
        return val

    def _snoop(self):
        """Diagnostic dump to console
        """
        print('Accessions:', self.accessions)
        print('Sequence:', self.sequence)
        print('Ambiguity:', self.ambiguity)
        print('Modifications:', self.modifications)
        for i, channel in enumerate(self.channels):
            print('%s: %s' % (channel, self.intensities[i]))
        print('Total Intensity:', self.total)
        print('Quan Info:', self.quant_info)
        print('Quan Usage:', self.quant_usage)
        print('q-Value:', self.qval)
        print('XCorr:', self.xcorr)
        print('DeltaCn:', self.delta_score)
        print('# Missed Cleavages:', self.missed)
        print('Spectrum File:', self.spectrum_file)
        print('First Scan:', self.first)
        print('Last Scan:', self.last)
        print('Charge:', self.charge)
        print('RT [min]:', self.rt)
        print('MH+ [Da]:', self.mhplus)
        print('Delta Mass [Da]:', self.deltamassda)
        print('Delta Mass [PPM]:', self.deltamassppm)
        print('MeetsAllCriteria:', self.meets_all_criteria)
        return

    def make_header(self):
        """Makes a header line for PSM data
        """
        header_line = (['Counter', 'Confidence Level', 'Sequence', 'PSM Ambiguity', 'Protein Descriptions',
                        '# Protein Groups', 'Protein Group Accessions', 'Modifications', 'Activation Type',
                        'DeltaScore', 'Rank']
                       + self.channels +
                       ['Total Int', 'Quan Info', 'Quan Usage', 'q-Value', 'PEP',
                        'Decoy Peptides Matched', 'Peptides Matched', 'XCorr', '# Missed Cleavages',
                        'Isolation Inteference [%]', 'Ion Inject Time [ms]', 'Intensity', 'Charge', 'm/z [Da]',
                        'MH+ [Da]', 'Delta Mass [Da]', 'Delta Mass [PPM]', 'RT [min]', 'First Scan',
                        'Last Scan', 'Spectrum File', 'MeetsAllCriteria', 'New Sequence', 'New Modifications',
                        'Peptide Length', 'PeptideMatchInfo', 'SEQUEST', 'length', 'ntt', 'ForR'])
        return '\t'.join(header_line)

    def make_data(self):
        """Makes a data line for PSM data
        """
        data_list = ([1, self.confidence_level, self.sequence, self.ambiguity, self.descriptions,
                      self.number_groups, self.accessions, self.modifications, self.activation,
                      self.delta_score, self.rank] +
                     self.intensities +
                     [self.total, self.quant_info, self.quant_usage, self.qval,
                      self.pep, self.decoy_candidates, self.target_candidates, self.xcorr, self.missed,
                      self.interference, self.inject_time, self.intensity, self.charge, self.moverz,
                      self.mhplus, self.deltamassda, self.deltamassppm, self.rt, self.first, self.last,
                      self.spectrum_file, self.meets_all_criteria, self.new_sequence, self.new_modifications,
                      self.peptide_length, self.match, self.SEQUEST, self.length, self.ntt, self.ForR])
        return '\t'.join([str(x) for x in data_list])

# end classes

def make_header_map(line, write, separator='\t'):
    """Returns a dictionary of "column header":index."""
    header_map = {}
    if not line:
        return header_map

    # split line at separator and populate dictionary (remove double quotes around header items)
    channels = []
    headers = line.split(separator)
    for i, header in enumerate(headers):
        header_map[header[1:-1]] = i
        if (header[1:-1].startswith('12') or header[1:-1].startswith('13')) and '/' not in header:
            channels.append(header[1:-1])

    # find the reporter ion channel headers
    for obj in write:
        print('channels are:', channels, file=obj)
           
    return header_map, channels
    # end
    
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
    if psm.accessions.endswith('"'):        # not sure if this is needed
        print('...WARNING: accessions had ending quote')
        psm.accessions = psm.accessions[:-1]

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
            elif item in db.skip:
                continue
            else:
                print('...WARNING: %s not in DB index' % item)
                db.skip.add(item)
        # make repaired accession and description strings
        acc_str = '; '.join(accessions)
        desc_str = '; '.join(descriptions)
        # save in db.seen dictionary
        db.seen[psm.accessions] = (acc_str, desc_str)
        # replace strings
        psm.descriptions = desc_str
        psm.accessions = acc_str                    

def parse_psm_lines(psm_file, db, write, max_qvalue=0.05, max_ppm=20.0,
                    max_interfere=50.0, separator='\t'):
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

    # start reading lines    
    start = False    
    for line in open(psm_file, 'r'):
        line = line.strip()
        if not line:
            continue    # skip blank lines        
        c.total += 1
        
        if not start:   # skip lines until header
            if line.startswith('"Confidence Level"'):   # look for header line
                psm_map, channels = make_header_map(line, write)
                start = True
                continue
            
        else:   # parse table lines
            psm = PSM(line, psm_map, channels, separator)
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
                psm.intensities = [0.0 for x in channels]

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
        temp = mod[:-1].split('(') # separate part inside ()
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
        temp = mod[:-1].split('(') # separate part inside ()
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
    """Makes SEQUEST-style sequence stings and removes static mods from modifications strings
    """
    new_seq = list(psm.sequence.upper())
    new_symbols = ['' for x in new_seq]
    modmask = get_variable_positions(psm, mod_type)
    for index in modmask:
        new_symbols[index] = modmask[index]
    psm.new_sequence = ''.join([j for i in zip(new_seq, new_symbols) for j in i])

    new_mod_list = [x.strip() for x in psm.modifications.split(';') if mod_type[x[:-1].split('(')[1]]]
    psm.new_modifications = '; '.join(new_mod_list)

def lookup_peptides(psm_list, proteins, prot_index):
    """Finds starting residue number for peptide sequence in protein sequence.
    """
    for psm in psm_list:
        try:
            # eventually add lookup of all accessions, just first to test
            acc = psm.accessions.split(';')[0].strip() # usual parsing
            try:
                psm.match = proteins[prot_index[acc]].findPeptide(psm.new_sequence)
                psm.start = psm.match[0][1]
            except KeyError:
                print('acc not in prot_index', acc)
                try:
                    acc = acc.split()[0]
                    psm.match = proteins[prot_index[acc]].findPeptide(psm.new_sequence) # phrog DB
                    psm.start = psm.match[0][1]
                except KeyError:
                    psm.valid_acc = False
        except IndexError:
            print('lookup failed:', psm.new_sequence, psm.accessions)
            print('\n...peptide lookup issue:')
            print('acc:', acc)
            print('index:', prot_index[acc])
            print('full acc:', proteins[prot_index[acc]].accession)
            print('peptide:', psm.new_sequence)
            print('peptide in sequence?', psm.new_sequence in proteins[prot_index[acc]].sequence, '\n')
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
print('Select the PSM export file')
psm_filename = PAW_lib.get_file(default_location,
                                [('Text files', '*.txt'), ('All files', '*.*')],
                                'Select a PD1.4 PSM export file')
if not psm_filename: sys.exit()     # exit if not file selected

# get the path so we can open a log file there
path, basename = os.path.split(psm_filename)
log_obj = open(os.path.join(path, 'make_PAW_TXT_from_PD1.4.log'), 'wt')
write = [None, log_obj]

for obj in write:
    print('====================================================================', file=obj)
    print(' make_PAW_TXT_from_PD1.4.py, %s, written by Phil Wilmarth, OHSU ' % VERSION, file=obj)
    print('====================================================================', file=obj)
    print('Ran on:', time.ctime(), file=obj)
    print('QVALUE = %s, PPM (+/-) = %s, INTERF = %s%%' % (QVALUE, PPM, INTERF), file=obj)

# get the FASTA database file
print('Select the FASTA protein database file')
db_filename = PAW_lib.get_file(default_location,
                                [('Fasta files', '*.fasta'), ('All files', '*.*')],
                                'Select a FASTA database')
if not db_filename: sys.exit()

for obj in write:
    print('PSM export file:', os.path.basename(psm_filename), file=obj)
    print('FASTA file:', os.path.basename(db_filename), file=obj)

# read in the FASTA file information and make some indexes
db = Fasta(db_filename, write)

# get the filtered psm information from the PSM export file
print('Reading PSMs from', os.path.split(psm_filename)[1])
counts, psm_list = parse_psm_lines(psm_filename, db, write, QVALUE, PPM, INTERF, SEPARATOR)

# analyze PTMs: fixed or variable? what residues? what names?
aa_freq, all_mods = analyze_modifications(psm_list)

# analyze mods and print summary
mod_type = fixed_or_variable(all_mods, aa_freq)
print_modification_report(all_mods, mod_type, write)

### repair accessions and descriptions to match what is in the FASTA file
##repair_accessions(psm_list, db)

# add alternatively formatted sequence string (SEQUEST style), reformat modifications string (remove static mods)
for psm in psm_list:
    fix_PTM_info(psm, mod_type)

# find peptide starting residue numbers
lookup_peptides(psm_list, db.proteins, db.prot_index)

# add true SEQUEST peptide string as attribute
add_SEQUEST_sequences(psm_list)

######################################################
# when everything is done, sort and write to new files

# open results file for all PSMs, print header lines
psmout = StringIO()
##dump = open(os.path.join(path, 'psm_dump.txt'), 'wt')
##print(psm_list[0].make_header(), file=dump)
print(psm_list[0].make_header(), file=psmout)

# print PSM data to file sorted by decreasing total intensity
psm_list = sorted(psm_list, key=lambda x: x.total, reverse=True)
for psm in psm_list:
    if psm.valid_acc:
        print(psm.make_data(), file=psmout)
##        print(psm.make_data(), file=dump)
##dump.close()

# read the table into a pandas data frame
psmout.seek(0)
df = pd.read_table(psmout)

# close file
psmout.close()

# copy the essentials to TXT data frame
txt_table = df[['First Scan', 'Last Scan', 'Charge', 'MH+ [Da]']].copy()
txt_table.columns = ['start', 'end', 'Z', 'expM']
txt_table['SpRank'] = df['Rank']
txt_table['theoM'] = df['MH+ [Da]'] + df['Delta Mass [Da]']
txt_table['deltaCN'] = df['DeltaScore']
txt_table['Xcorr'] = df['XCorr']
txt_table['Sequence'] = df['SEQUEST']
txt_table['Loci'] = df['Protein Group Accessions']
txt_table['NewDeltaCN'] = df['DeltaScore']
txt_table['NewDisc'] = df['q-Value']
txt_table['ntt'] = df['ntt']
txt_table['ForR'] = df['ForR']

# add the reporter ion peak heights
cols = [col for col in df.columns.values if (col.startswith('12') or col.startswith('13'))]
cols = [x.replace('-N', '_N') for x in sorted([x.replace('_N', '-N') for x in cols])] # want tags in mass order (N before C)
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
    print('\nThere were', counts.total, 'Total rows in PSM table export', file=obj)
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
