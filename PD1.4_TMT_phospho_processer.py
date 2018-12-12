""" program "PD1.4_TMT_phospho_processer_v2.py"
    Process TMT export files from Proteome Discoverer. For proteins
    with a minimum number of unique peptides, total reporter ion
    intensities are computed.

    added median reporter ion intensity cutoff, spring 2015 -PW
    added deltamass filter, spring 2015 -PW
    double checked things and added some additional comments, Nov. 2015, -PW
    added reformatted PSM data dump, Nov. 2015, -PW
    parses extra data from phospho searches, Feb. 2016, -PW
    site positions were off by one (too large), Aug. 2016, -PW
    
    written by Phil Wilmarth 2015.

    3/20/2017 (PW): added flag to skip or not skip unmodified peptides
    6/15/2017 (PW): moved missing data inputation to the combined PSMs only; added DB to results files; 

    TODO:
    rewrite intensity functions to return test value and do the
    testing in the calling function/method.
    
"""
# import common modules here
import os
import sys
import time
import operator
import copy
import PAW_lib

# column separator character
SEPARATOR = '\t'

# missing value replacement (for combined PSMs)
MISSING = 150.0

# minimum trimmed average intensity
INTENSITY = 500.0

# q-Value cutoff
QVALUE = 0.01

# plus/minus PPM deltamass window
PPM = 20.0

# skip unmodified peptides?
SKIP_UNMOD = True

def make_header_map(line, separator='\t'):
    """Returns a dictionary of "column header":index
    """
    header_map = {}
    if not line:
        return header_map

    # split line at separator and populate dictionary (remove double quotes around header)
    headers = line.split(separator)
    for i, header in enumerate(headers):
        header_map[header[1:-1]] = i
        
    return header_map
    # end

def test_median_intensity(intensities, intensity, zero_input):
    """Find median intensity of channels, excluding "empty" channels.
    Return True if median above intensity, False otherwise.
    """
    int_vector = [x for x in intensities if x > zero_input]

    # exclude PSM if too few channels with minimal intensity
    if len(int_vector) < 3:
        return False

    # compute the median
    int_vector.sort()
    length = len(int_vector)
    if int_vector:
        if not length % 2:
            median = (int_vector[length/2] + int_vector[length/2 - 1]) / 2.0
        else:
            median = int_vector[length/2]
    else:
        median = 0.0
        
    # test threshold
    if median >= intensity:
        return True
    else:
        return False
    
def test_average_intensity(intensities, intensity, zero_input):
    """Find average intensity of channels, excluding "empty" channels.
    Return True if average above "intensity", False otherwise.
    """
    int_vector = [x for x in intensities if x > zero_input]

    # exclude PSM if too few channels with minimal intensity
    if len(int_vector) < 3:
        return False

    # compute the average
    average = sum(int_vector)/float(len(int_vector))
        
    # test threshold
    if average >= intensity:
        return True
    else:
        return False
    
def test_trimmed_average_intensity(intensities, intensity):
    """Parse line and get intensities. Find average intensity of channels,
    excluding the top and bottom values. It does not skip any low values.
    Return True if average above "intensity", False otherwise.
    """
    int_vector = intensities

    # compute the trimmed average
    average = sum(int_vector[1:-1])/float(len(int_vector[1:-1]))
        
    # test threshold
    if average >= intensity:
        return True
    else:
        return False
    
class PSM():
    """Container for PSM information
    """
    def __init__(self, line, header_map, missing=0.0, separator='\t'):
        """Parse line and populate fields
        """
        items = line.split(separator)
        self.channels = ['126', '127_N', '127_C', '128_N', '128_C',
                         '129_N', '129_C', '130_N', '130_C', '131']

        # note that PD exports fields enclosed in double quotes
        self.confidence_level = items[header_map['Confidence Level']][1:-1]
#        self.search_id = items[header_map['Search ID']][1:-1]
#        self.processing_node = items[header_map['Processing Node No']][1:-1]
        self.sequence = items[header_map['Sequence']][1:-1]
        self.peptide_length = len(self.sequence)    # compute peptide length
#        self.sequence_id = self._int(items[header_map['Unique Sequence ID']][1:-1])
        self.ambiguity = items[header_map['PSM Ambiguity']][1:-1]
        self.descriptions = items[header_map['Protein Descriptions']][1:-1]
#        self.number_proteins = self._int(items[header_map['# Proteins']][1:-1])
        self.number_groups = self._int(items[header_map['# Protein Groups']][1:-1])
        self.accessions = items[header_map['Protein Group Accessions']][1:-1]
        self.modifications = items[header_map['Modifications']][1:-1]
        self.number_phospho = self.modifications.count('Phospho')   # add count of phospho sites
        self.activation = items[header_map['Activation Type']][1:-1]
        self.delta_score = self._float(items[header_map['DeltaScore']][1:-1])
        if self.delta_score == 1.0:
            self.delta_score = 0.0
#        self.deltacn = self._float(items[header_map['DeltaCn']][1:-1])
        self.rank = self._int(items[header_map['Rank']][1:-1])
#        self.se_rank = self._int(items[header_map['Search Engine Rank']][1:-1])
        try:
            self.isoform_prob = items[header_map['phosphoRS Isoform Probability']][1:-1]
            self.isoform_prob = float(self.isoform_prob.replace('%', ''))/100.0
        except:
            self.isoform_prob = 0.0
        try:
            self.site_prob = items[header_map['phosphoRS Site Probabilities']][1:-1]
        except:
            self.site_prob = ''
        self.intensities = [0.0 for x in self.channels]
        for i, channel in enumerate(self.channels):
            self.intensities[i]  =  self._float(items[header_map[channel]][1:-1], missing)
        self.total = sum(self.intensities)
        self.quant_info = items[header_map['Quan Info']][1:-1]
        self.quant_usage = items[header_map['Quan Usage']][1:-1]
#        self.quant_result_id = self._int(items[header_map['QuanResultID']][1:-1])
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
#        self.ms_order = items[header_map['MS Order']][1:-1]
#        self.ions_matched = items[header_map['Ions Matched']][1:-1]
#        self.matched_ions = items[header_map['Matched Ions']][1:-1]
#        self.total_ions = items[header_map['Total Ions']][1:-1]
        self.spectrum_file = items[header_map['Spectrum File']][1:-1]
#        self.annotation = items[header_map['Annotation']][1:-1]    # seems empty
        self.meets_all_criteria = False
        self.grouped_psms = ''
        self.psm_count = 0
        self.localization = 'NA'
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
                        'DeltaScore', 'Rank', 'phosphoRS Isoform Probability', 'phosphoRS Site Probabilities']
                       + self.channels +
                       ['Total Int', 'Quan Info', 'Quan Usage', 'phosphoRS Binomial Peptide Score', 'q-Value', 'PEP',
                        'Decoy Peptides Matched', 'Peptides Matched', 'XCorr', '# Missed Cleavages',
                        'Isolation Inteference [%]', 'Ion Inject Time [ms]', 'Intensity', 'Charge', 'm/z [Da]',
                        'MH+ [Da]', 'Delta Mass [Da]', 'Delta Mass [PPM]', 'RT [min]', 'First Scan',
                        'Last Scan', 'Spectrum File', 'PSM Number', 'MeetsAllCriteria', 'New Sequence', 'New Modifications',
                        'Number Phospho Sites', 'Peptide Length', 'New Site Prob Peptide', 'New Site Prob Protein',
                        'Max Prob', 'Min Prob', 'Site List', 'Localization Status', 'Group String', 'Grouped PSMs',
                        'Used PSM Count', 'PeptideMatchInfo'])
        return '\t'.join(header_line)

    def make_data(self):
        """Makes a data line for PSM data
        """
        data_list = ([1, self.confidence_level, self.sequence, self.ambiguity, self.descriptions,
                      self.number_groups, self.accessions, self.modifications, self.activation,
                      self.delta_score, self.rank, self.isoform_prob, self.site_prob] +
                     self.intensities +
                     [self.total, self.quant_info, self.quant_usage, self.binomial_score, self.qval,
                      self.pep, self.decoy_candidates, self.target_candidates, self.xcorr, self.missed,
                      self.interference, self.inject_time, self.intensity, self.charge, self.moverz,
                      self.mhplus, self.deltamassda, self.deltamassppm, self.rt, self.first, self.last,
                      self.spectrum_file, self.psm_number, self.meets_all_criteria, self.new_sequence, self.new_modifications,
                      self.number_phospho, self.peptide_length, self.new_site_prob_peptide, self.new_site_prob,
                      self.max_prob, self.min_prob, self.sites, self.localization, self.grouper, self.grouped_psms,
                      self.psm_count, self.match])
        return '\t'.join([str(x) for x in data_list])

    def make_key(self):
        """Makes a column header key
        """
        return """\nColumn Header\tDescription
Counter\tColumn of 'ones' for counting rows
Confidence Level\tPD confidence level flag
Sequence\tPepetide sequence (PD 1.4 format)
PSM Ambiguity\tWhether peptide is shared or unique
Protein Descriptions\tProtein description strings
# Protein Groups\tNumber of proteins in group
Protein Group Accessions\tProtein accession(s)
Modifications\tModifications (PD 1.4 format)
Activation Type\tActivation (fragmentation) type
DeltaScore\tDelta CN vlaue from SEQUEST
Rank\tPSM rank from search engine
phosphoRS Isoform Probability\tPhospoRS localization probability (%)
phosphoRS Site Probabilities\tPhospoRS site localization probabilities
TotInt 126\tIntensity for respective TMT channel
TotInt 127_N\tIntensity for respective TMT channel
TotInt 127_C\tTotal protein intensity (sum over PSMs) for respective TMT channel
TotInt 128_N\tIntensity for respective TMT channel
TotInt 128_C\tIntensity for respective TMT channel
TotInt 129_N\tIntensity for respective TMT channel
TotInt 129_C\tIntensity for respective TMT channel
TotInt 130_N\tIntensity for respective TMT channel
TotInt 130_C\tIntensity for respective TMT channel
TotInt 131\tIntensity for respective TMT channel
Total Intensity\tSum across reporter ion intensities
Quan Info\tQuantitation Information tag
Quan Usage\tWhether or not peptide is useable for TMT quantitation
q-Value\tPSM q-value computer by Percolator
PEP\tPosterior error probability
Decoy Peptides Matched\tNumber of decoy candidates
Peptides Matched\tNumber of target peptide candidates
XCorr\tPSM XCorr value from SEQUEST
# Missed Cleavages\tNumber of missed enzymatic cleavages
Isolation Interference [%]\tEstimated precursor isolation  purity
Ion Inject Time [ms]\tIon inject time in milliseconds
Intensity\tPrecursor intensity
Charge\tPeptide charge
MH+ [Da]\tMeasure MH+ value (Da)
Delta Mass [Da]\tDiffernece betwen measured and calcuated masses (Da)
Delta Mass [PPM]\tDiffernece betwen measured and calcuated masses (ppm)
RT [min]\tRetention time for MS2 scan
First Scan\tStarting scan number
Last Scan\tEnding scan number
Spectrum File\tRAW file base name
Annotation\tWho knows
PSM Number\tArbitrary PSM number to correlate grouped and ungrouped PSMs
MeetsAllCriteira\tTrue if PSM meets a set of validation criteria (q-score, length, charge range, deltamass, etc.)
New Sequence\tNew peptide sequence string in SEQUEST format
New Modifications\tModifications string without static mods
Number Phospho Sites\tNumber of phospho sites in peptide
Peptide Length\tNumber of amino acids in peptide sequence
New Site Prob Peptide\tSimplified phosphoRS site probabilities (reduced to number of phosphosites from high to low p-val), relative positions
New Site Prob Protein\tSimplified phosphoRS site probabilities (reduced to number of phosphosites from high to low p-val), aa positions in protein (probability ties are summed)
Max Prob\tMaximum site probability
Min Prob\tMinimum site probability (cut at number of sites)
Site List\tList of localized sites without any probabilities
Localization Status\tConsistant string means that all PSM localizations have similar sites, Varies denotes otherwise
Group String\tString used to group similar PSM (sequence plus nominal MH+ mass)
Used PSM Count\tThe number of PSMs contributing to the intensity totals
Match\tMatch tuple from peptide lookup in protein sequence (multiple tuples if peptide is present more than once in protein)

"""
    # end class

def parse_psm_lines(psm_file, max_qvalue=0.05, max_ppm=20.0, min_intensity=500, missing=0.0, separator='\t'):
    """Parses PSM exports and adds information to PSM objects.
    Returns a list of all PSMs passing q-value, deltamass, and intensity cutoffs.
    """
    # define some fixed ranges of valid PSM attributes
    min_charge = 2
    max_charge = 4
    min_length = 7
    max_length = 40

    # initialize counters, etc.    
    qval_good = 0
    qval_bad = 0
    in_ppm = 0
    out_ppm = 0
    above_int = 0
    below_int = 0
    total = 0
    top = 0
    valid = 0
    reject = 0
    psm_list = []
    
    start = False    
    for line in open(psm_file, 'r'):
        line = line.strip()
        if not line:
            continue    # skip blank lines
        
        total += 1        
        if not start:   # skip lines until header
            if line.startswith('"Confidence Level"'):   # look for header line
                psm_map = make_header_map(line)
                start = True
                continue
        else:   # parse table lines
            psm = PSM(line, psm_map, missing, separator)
##            if not psm:
            if not psm.accessions: # this may filter out low q-value matches
                total += -1
                reject += -1
                continue

            # PSM tables have extra lines for non-top hits to get deltaCN values, skip those
            if psm.rank == 1:
                top += 1
            else:
                continue

            # add real PSMs to list
            psm_list.append(psm)

            # test various criteria
            if psm.qval <= max_qvalue:
                qval_good += 1
            else:
                qval_bad += 1
                reject += 1
                continue
                
            if abs(psm.deltamassppm) <= max_ppm:
                in_ppm += 1
            else:
                out_ppm += 1
                reject += 1
                continue

            if test_trimmed_average_intensity(psm.intensities, min_intensity):
                above_int += 1
            else:
                below_int += 1
                reject += 1
                continue

            if min_charge <= psm.charge <= max_charge:
                pass
            else:
                reject += 1
                continue

            if min_length <= psm.peptide_length <= max_length:
                pass
            else:
                reject += 1
                continue

            psm.meets_all_criteria = True
            valid += 1
                    
    return [total, top, valid, reject, qval_good, qval_bad, in_ppm, out_ppm, above_int, below_int], psm_list

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
    """Assigns special symbols to variable mods. None assigned to fixed mods.
    """
    symbols = ['*', '#', '@', '^', '~', '$', '[', ']']

    mod_type = {}
    global_mod_freq = {}
    for mod in all_mods:
        for residue in all_mods[mod]:
            if all_mods[mod][residue] == aa_freq[residue]:
                mod_type[mod] = None
            else:
                mod_type[mod] = True
    for mod in all_mods:
        if mod_type[mod]:
            global_mod_freq[mod] = sum(all_mods[mod].values())
    variable_freq = sorted(global_mod_freq.items(), key=lambda x: x[1], reverse=True)
    for i, (mod, count) in enumerate(variable_freq):
        mod_type[mod] = symbols[i]
    return mod_type

def get_variable_positions(psm, mod_type):
    """Parses PD modification descriptions to get modification positions and symbols
    """
    modmask = {}

    # split modification description string
    modlist = psm.modifications.split(';')
    for mod in modlist:
        mod = mod.strip()   # get rid of whitespace
        temp = mod[:-1].split('(') # separate part inside ()
        position = ''.join([c for c in temp[0] if c.isdigit()]) # get positions
        try:
            position = int(position)
        except ValueError:
            if temp[0] == 'N-Term':
                position = 0
            elif temp[0] == 'C-Term':
                position = len(psm.sequence)
            else:
                position = -1
        modtype = temp[1]

        if mod_type[modtype]:   # just variable mods
            modmask[position-1] = mod_type[modtype]

    return modmask

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

def sequence_length_analysis(psm_list):
    length_count = [0 for x in range(60)]
    print()
    for i in range(60):
        length_count[i] = len([x for x in psm_list if x.peptide_length == i])

    total = float(sum(length_count))
    running_tot = 0.0
    for i in range(60):
        running_tot += length_count[i]
        print(i, length_count[i], round(100.0*running_tot/total, 2))
    print()
    return

class RSProb:
    """Container for parsed PhosphoRS data
    """
    def __init__(self, one_site):
        try:
            self.residue = [one_site.split('(')[0]]   # amino acid (S, T, or Y)
            self.position = [int(one_site.split(':')[0][:-1].split('(')[1])]   # relative position in peptide string
            self.probability = float(one_site.split(':')[1])  # assigned probability in percent
        except:
            print('...WARNING: site description:', one_site)
##            self.residue = ''
##            self.position = 0
##            self.probability = 0.0

def parse_RSProb(psm):
    """Parse phosphoRS probability string and returns N sites of the highest probability.
    """
    RS_full_list = []
    if ((not psm.site_prob) or
        (psm.site_prob == 'Too many NL-allowing PTMs') or
        (psm.site_prob == 'Too many isoforms')):
        psm.new_site_prob = ''
        psm.new_site_prob_peptide = ''
        psm.max_prob = 0.0
        psm.min_prob = 0.0
        psm.sites = ''
        return
    for prob in psm.site_prob.split(';'):
        prob = prob.strip()
        prob = RSProb(prob)     # this can return empty object for odd PhosphoRS site descriptions
        prob.position[0] += (psm.start-1)
        RS_full_list.append(prob)        

    # sort descending and get the top N sites
    RS_full_list.sort(reverse=True, key=lambda x: x.probability)
    RS_list = RS_full_list[:psm.number_phospho]
    psm.max_prob = RS_list[0].probability
    psm.min_prob = RS_list[-1].probability

    # check for prob ties
    for rs in RS_full_list[psm.number_phospho:]:
        if rs.probability == RS_list[-1].probability:
            RS_list[-1].residue += rs.residue
            RS_list[-1].position += rs.position
    
    # put together output string           
    site_list_prob = []
    for rs in RS_list:
        residue_position = ', '.join(rs.residue) + ' (' + ', '.join([str(x) for x in rs.position]) + ')'
        # maybe combine probabliities for ties?
        site_list_prob.append('%s: %0.1f' % (residue_position, len(rs.position)*rs.probability))
    psm.new_site_prob = '; '.join(site_list_prob)
#    psm.sites = '; '.join(site_list)
    
    # put together output string with relative positions          
    site_list_prob = []
    for rs in RS_list:
        positions = [x-psm.start+1 for x in rs.position]
        residue_position = ', '.join(rs.residue) + ' (' + ', '.join([str(x) for x in positions]) + ')'
        site_list_prob.append('%s: %0.1f' % (residue_position, rs.probability))
    psm.new_site_prob_peptide = '; '.join(site_list_prob)

    # make site list (no ties), sorted by residue position
    site_set = set()
    for rs in RS_list:
        [site_set.add(x) for x in zip(rs.residue, rs.position)]
    site_list = sorted(list(site_set), key=lambda x: x[1])
    site_list = [str(x[0])+str(x[1]) for x in site_list]
    psm.sites = '; '.join(site_list)

    return

def make_protein_index(proteins):
    """Indexes proteins
    """
    prot_index = {}
    skip = set(['sp', 'tr', 'gi', 'ref', ''])
    for i, p in enumerate(proteins):
        accs = p.accession.split('|')
        for acc in accs:
            if acc in skip:
                continue
            prot_index[acc] = i
    return prot_index        

def lookup_peptides(psm_list, proteins, prot_index):
    """Finds starting residue number for peptide sequence in protein sequence.
    """
    for psm in psm_list:
        try:
            # eventually add lookup of all accessions, just first to test
            acc = psm.accessions.split(';')[0].strip()
            psm.match = proteins[prot_index[acc]].findPeptide(psm.new_sequence, pad_count=3)
            psm.start = psm.match[0][1]
        except IndexError:
            print()
            print('...peptide lookup issue:')
            print('pre-acc:', psm.accessions)
            print('acc:', acc)
            print('index:', prot_index[acc])
            print('full acc:', proteins[prot_index[acc]].accession)
            print('peptide:', psm.new_sequence)
            print('peptide in sequence?', psm.new_sequence in proteins[prot_index[acc]].sequence)
            print()
            psm.start = 0
        
def make_group_index(psm_list):
    """Makes a dictionary of PSM index lists (values) that have the same grouper string (key)
    """
    group_index = {}
    for i, psm in enumerate(psm_list):
        if psm.meets_all_criteria:
            if psm.grouper in group_index:
                group_index[psm.grouper].append(i)
            else:
                group_index[psm.grouper] = [i]
    return group_index

def analyze_modifications(psm_list):
    """gets freqeuncies of amino acids in sequences, gets frequencies of modifications
    """
    aa_freq = {}
    all_mods = {}
    
    aa_freq['N-Term'] = aa_freq['C-Term'] = len(psm_list)
    for psm in psm_list:
        amino_acid_frequency(psm.sequence, aa_freq)
        mods = parse_mods(psm.modifications)
        update_dictionary(all_mods, mods)

    return aa_freq, all_mods

def print_modification_report(all_mods, mod_type):
    """Prints a summary of variable and static modifications
    """
    variable = []
    print('\nVariable modifications:')   # print mod type, symbol, affected residues
    for mod in mod_type:
        if mod_type[mod]:
            variable.append([mod, mod_type[mod], sorted(all_mods[mod].keys()), sum(all_mods[mod].values())])
    variable = sorted(variable, reverse=True, key=lambda x: x[3])   # order mods by decreasing frequency
    for mod in variable:
        print('  ', mod[0], mod[1], mod[2])

    static = []
    print('Static modifications:')   # print mod type and afected residues
    for mod in mod_type:
        if not mod_type[mod]:
            static.append([mod, sorted(all_mods[mod].keys()), sum(all_mods[mod].values())])
    static = sorted(static, reverse=True, key=lambda x: x[2])   # order mods by decreasing freqiuency
    for mod in static:
        print('  ', mod[0], mod[1])
    print()
    return

def test_phosphoRS_localization(psm_list):
    """Need some localization probability cutoff testing here.
    Should depend on number of sites in peptide.
    """
    for psm in psm_list:
        if psm.number_phospho == 0:
            if SKIP_UNMOD:
                psm.meets_all_criteria = False  # this rejects unmodified peptides
        elif psm.number_phospho == 1 and psm.min_prob < 45.0:
            psm.meets_all_criteria = False
        elif psm.number_phospho == 2 and psm.min_prob < 0.30:
            psm.meets_all_criteria = False
        elif psm.number_phospho == 3 and psm.min_prob < 22.5:
            psm.meets_all_criteria = False

def make_grouped_PSMs(psm_list, group_index, missing=150.0):
    """Groups PSMs by grouper key, keeps psm with best q-value, sums intensities
    """
    new_psm_list = []
    for key in group_index:
        if len(group_index[key]) == 1:
            # nothing to group
            best_psm = psm_list[group_index[key][0]]
            best_psm.localization = 'consistent'
            best_psm.grouped_psms = str(best_psm.psm_number)
            best_psm.psm_count = 1
            new_psm_list.append(best_psm)
        else:
            # pick group member with best q-value to represent group
            psm_group_list = sorted([psm_list[i] for i in group_index[key]], key=lambda x: x.qval)
            best_psm = copy.deepcopy(psm_group_list[0])
            best_psm.psm_count = len(psm_group_list)
            best_psm.grouped_psms = '; '.join([str(x.psm_number) for x in psm_group_list])
            # sum the intensities for the group
            summed = [0.0 for i in best_psm.channels]
            for i in range(len(best_psm.channels)):
                for j in range(len(psm_group_list)):
                    summed[i] += psm_group_list[j].intensities[i]
            best_psm.intensities = copy.deepcopy(summed)
            # do zero replacement
            best_psm.intensities = [(x if x > 0.0 else missing) for x in best_psm.intensities]
            # sum reporter ions
            best_psm.total = sum(best_psm.intensities)
            # see if site localization is the same for all PSMs
            site_set = set([x.sites for x in psm_group_list])
            if len(site_set) == 1:
                best_psm.localization = 'consistent'
            else:
                best_psm.localization = 'varies'
            new_psm_list.append(best_psm)
    return new_psm_list
            
    

###########################################
######## main program starts here #########
###########################################

print('\n====================================================')
print(' program "PD_TMT_phospho_processer.py", version 2.0 ')
print('====================================================')
print('\nRan on:', time.ctime())
print('INTENSITY = %s, QVALUE = %s, MISSING = %s, PPM = %s' % (INTENSITY, QVALUE, MISSING, PPM))

# get the PSM results PD export file information
default_location = r'F:\PSR_Core_Analysis'
if not os.path.exists(default_location):
    default_location = os.getcwd()
print('Select the PSM export file')
psm_filename = PAW_lib.get_file(default_location,
                                [('Text files', '*.txt'), ('All files', '*.*')],
                                'Select a phospho PSM export file')
if not psm_filename: sys.exit()     # exit if not file selected

# get the FASTA database file
print('Select the FASTA protein database file')
db_filename = PAW_lib.get_file(default_location,
                               [('Fasta files', '*.fasta'), ('All files', '*.*')],
                               'Select a FASTA database')
if not db_filename: sys.exit()

# read in the protein sequences
print('Reading proteins...')
f = PAW_lib.FastaReader(db_filename)
p = PAW_lib.Protein()
proteins = []
while f.readNextProtein(p, False):
    prot = copy.deepcopy(p)
    proteins.append(prot)
prot_index = make_protein_index(proteins)

# get the psm information from the PSM export file
counts, psm_list = parse_psm_lines(psm_filename, QVALUE, PPM, INTENSITY, 0.0, SEPARATOR)

# analyze PTMs: fixed or variable? what residues? what names?
aa_freq, all_mods = analyze_modifications(psm_list)

# analyze mods and print summary
mod_type = fixed_or_variable(all_mods, aa_freq)
print_modification_report(all_mods, mod_type)    

# add alternatively formatted sequence string (SEQUEST style), reformat modifications string (remove static mods)
for psm in psm_list:
    fix_PTM_info(psm, mod_type)

# find peptide starting residue numbers (do before RS probabilities)
lookup_peptides(psm_list, proteins, prot_index)

# add a simplified RS probability string (top N sites, check for ties)
for psm in psm_list:
    parse_RSProb(psm)

# test for some minimal phosphoRS localization and update "meets_all_criteria"
# need to set "meets_all_criteria" since it controls combining of similar PSMs
test_phosphoRS_localization(psm_list)

# add PSM number for indexing purposes (DO NOT sort psms until after grouping is done)
for i, psm in enumerate(psm_list):
    psm.psm_number = i

# figure out which PSMs to group together (sum intensities)
group_index = make_group_index(psm_list)    # this is a filtered list ("meets_all_criteria")
new_psm_list = make_grouped_PSMs(psm_list, group_index, MISSING)    # do the zero replacement here

######################################################
# when everything is done, sort and write to new files

# open results file for all PSMs, print header lines
psmout = open(psm_filename.replace('_psms', '_psm_filtered_all'), 'w')
print(psm_list[0].make_header(), file=psmout)

# print PSM data to file sorted by decreasing total intensity
psm_list = sorted(psm_list, key=lambda x: x.total, reverse=True)
for psm in psm_list:
    print(psm.make_data(), file=psmout)

# print keys at end of table (psm should still be last psm in psm_list)
print(psm.make_key(), file=psmout)

for out in [None, psmout]:
    # print the parsing statistics
    print('The FASTA database was:', db_filename, file=out)
    print('There were', counts[0], 'Total rows in PSM table export', file=out)
    print('There were', counts[1], 'Top ranked PSMs', file=out)
    print('There were', counts[2], 'PSMs that passed all criteria', file=out)
    print('There were', counts[3], 'PSMs that were rejected', file=out)
    print('There were', counts[4], 'PSMs passing q-value cutoff', file=out)
    print('There were', counts[5], 'PSMs that failed q-value cutoff', file=out)
    print('There were', counts[6], 'PSMs passing deltamass cutoff', file=out)
    print('There were', counts[7], 'PSMs that failed deltamass cutoff', file=out)
    print('There were', counts[8], 'PSMs with TMT ions above intensity threshold', file=out) 
    print('There were', counts[9], 'PSMs with TMT ions below intensity threshold', file=out)
    print('(%0.2f%% of scans rejected)' % (100.0 * counts[9] / (counts[8]+counts[9]),), file=out)

# close file
psmout.close()

# open results file for combined PSMs, print header lines
psmout = open(psm_filename.replace('_psms', '_psm_filtered_combined'), 'w')
print(new_psm_list[0].make_header(), file=psmout)

# print PSM data to file sorted by decreasing total intensity
new_psm_list = sorted(new_psm_list, key=lambda x: x.total, reverse=True)
for psm in new_psm_list:
    print(psm.make_data(), file=psmout)

# print keys at end of tables
print('The FASTA database was:', db_filename, file=psmout)
print(psm.make_key(), file=psmout)
    
# close files
psmout.close()

# end
