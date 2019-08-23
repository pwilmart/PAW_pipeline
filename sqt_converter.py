"""sqt_converter.py written by Phil Wilmarth, OHSU.

    original version 2008 by PW supported SEQUEST only.
    support for Mascot and X!Tandem added in 2012 by PW.
    support for Comet added 2014 by PW.
    rewritten as callable function by Delan Huang, 2014.
    added support for Lys-C enzyme, Oct. 27, 2014 by PW.
    added support for non-standard Comet SQT files, Oct. 31, 2014 by PW.
        (Comet does not list multiple L lines)
    updated calls to use new library files, Oct. 28, 2015 by PW.
    peptide lookup redone to use theretical peptide digests, April 2016 -PW.
    removed any older SEQUEST-related code, April 2016 -PW
    
    updated for python 3 -PW 9/16/2017
    adds reporter ion quantities to TXT files -PW 10/12/2017
"""
import os
import glob
import gzip
import sys
import time
import math
import copy
import re

from collections import OrderedDict

import PAW_lib
try:
    from configuration import *
except:
    default_location = r'F:\PSR_Core_Analysis'
    decoy_string = 'REV'
    min_pep_len = 7

DECOY_STRING = decoy_string
MIN_PEP_LEN = min_pep_len
DEFAULT_LOCATION = default_location
MASK = False    # I think Comet list all sequences explicitly
VERSION = 'v2.0.3'

##############################################################################################
# routines for discriminant scores, alternate DeltaCNs, and data structures for search results
##############################################################################################

def new_discrim_score(Z, Xcorr, length, newDCN, SpRank, MH1exp, MH1calc):
    """Calculate modified discriminant score using newDCN.  Returns a float.

    Usage: newdiscrim =
        discriminant_score(Z, Xcorr, length, deltaCN, SpRank, MH1exp, MH1calc),
        where "Z" is peptide charge, "Xcorr" is SEQUEST XCorr score,
        "length" is peptide length, "deltaCN" is traditional DCN score,
        "SpRank" is the preliminary score ranking, "MH1exp" is the measured
        MH+ mass, "MH1calc" is the calculated MH+ mass, and "discrimval" is
        the returned discriminant function value.

    Coefficients optimized for newDCN (based on average of hits 4 to 12
       instead of hit #2) using in-house LTQ control mix, 10/28/08 -PW.
    """
    # trap low values and set minimums
    if length < 2:
        length = 2
    if Xcorr < 0.01:
        Xcorr = 0.01
    
    # 1+ peptides
    if Z == 1:
        if length >= 100: length = 100
        numions = 2 * length
        xcorrprime = math.log(Xcorr) / math.log(numions)
        newdiscrim = (10.0 * xcorrprime +
                      7.0 * newDCN +
                      (-0.25) * math.log(SpRank) +
                      (-0.45) * math.fabs(MH1exp - MH1calc) + (-.05))
    
    # 2+ peptides
    elif Z == 2:
        if length >= 15: length = 15
        numions = 2 * length
        xcorrprime = math.log(Xcorr) / math.log(numions)
        newdiscrim = (8.0 * xcorrprime +
                      7.5 * newDCN +
                      (-0.20) * math.log(SpRank) +
                      (-0.55) * math.fabs(MH1exp - MH1calc) + (-0.80))
    
    # 3+ and 4+ peptides
    elif Z == 3 or Z == 4:
        if length >= 25: length = 25
        numions = 4 * length
        xcorrprime =  math.log(Xcorr) / math.log(numions)
        newdiscrim = (10.0 * xcorrprime +
                      8.0 * newDCN +
                      (-0.20) * math.log(SpRank) +
                      (-0.40) * math.fabs(MH1exp - MH1calc) + (-1.5))
        
    # set higher charge states to large negative score
    else:
        newdiscrim = -10.0
    
    # return the new, modified discriminant score
    return newdiscrim

def get_new_old_dcn(matches, top):
    """Computes a new deltaCN by averaging scores for hits 4 to 12 and
    using that Xcorr instead of Xcorr of second best match. More robust
    when PTMs are specified (less sensitive to positional isomers).
    """
    # comment this out to get "newdcn = 0" with no exit when there are too few matches reported
    if int(top) < 4:
        print('get_new_old_dcn ERROR: not enough top hits to compute new deltaCN values!')
        sys.exit()
        
    dcn = 0.0       # traditional deltaCN
    ave_xcorr = 0.0 # average of scores from hits 4 to "top"
    newdcn = 0.0    # new deltaCN value
    scores = 0      # counter
    oldmatch = 0
    for match in matches:
        newmatch = match.match
        
        # traditional deltaCN is from first match where Xcorr != top Xcorr
        if (dcn == 0.0) and (matches[0].xcorr != match.xcorr):
            dcn = float(match.deltaCN)
        
        # average hits 4 to "top" (if present)
        if (int(match.XcRank) > 3 and newmatch <= int(top)): 
            if newmatch != oldmatch:
                ave_xcorr = ave_xcorr + float(match.xcorr)
                scores += 1
        oldmatch = newmatch
    
    # there were some hits between 4 and "top" so compute newDCN
    if scores != 0: 
        ave_xcorr = ave_xcorr/scores
        newdcn = (float(matches[0].xcorr) - ave_xcorr)/float(matches[0].xcorr)
    else:
        newdcn = dcn
    #            
    return newdcn, dcn  # return new and old DCN values

class SequestMatch:
    """Object to hold the main data present in a SEQUEST OUT file.
    written by Phil Wilmarth, OHSU, 2008.
    """
    def __init__(self):
        """Set up the basic attributes in the constructor, no parameters passed in.
        """
        self.match = '1'            # match number
        self.xcorr = '2.5000'       # Xcorr score
        self.deltaCN = '0.0000'     # DeltaCN score
        self.XcRank = '1'           # Rank by Xcorr score
        self.SpRank = '1'           # Rank by Sp score
        self.Sp = '250.0'           # Sp score or expectation value
        self.MH1calc = '1000.0'     # Calculated mass
        self.ionmatch = '10'        # Number of matched ions
        self.ionpredict = '10'      # Number of predicted ions
        self.accession = ''         # Protein accession of first match
        self.sequence = ''          # Peptide sequence (SEQUEST format)
        self.match_tuple_list = []  # List of protein match tuples
        self.additional = 0      # most attributes are STRINGS, this is integer

    def _snoop(self):
        """Prints some useful diagnostic information.
        """
        print('...', self.xcorr, self.deltaCN, self.accession, self.sequence, self.additional)
        return

    # end class

class OutData:
    """Object to hold all of the relevant info from SEQUEST search (SQT) files.
    Written by Phil Wilmarth, OHSU, 2008, 2016.
    """
    def __init__(self):
        """Basic constuctor sets several attributes.  No parameters.
        """
        self.beg = '100'            # beginning scan number
        self.end = '100'            # ending scan number
        self.z = '2'                # peptide charge state 
        self.MH1exp = '1000.0'      # experimental precursor MH+ value
        self.top = '12'             # no. of peptide matches in SQT file
        self.matches = []           # list of "SequestMatch" objects
        return

    def _snoop(self):
        """
        Prints out some info for diagnostics.
        """
        print(self.base_name, self.beg, self.end, self.z)
        for i in self.matches:
            print(i.xcorr, i.deltaCN, i.sequence, i.accession)
        return
        
    # end class
    
def make_PAW_txt_file(sqt_file, outs, params):
    """
    Makes a PAW text file (top hit summary) from list of Out objects
    (loaded from one SQT file).
    """
    for obj in params.log_obj:
        print('...starting TXT file creation at:', time.ctime(), file=obj)
    print('length of outs:', len(outs))
    
    # open file or zipped file
    txt_name = sqt_file.replace('.sqt', '.txt')
    if params.zipflag:
        txt_gz = open(txt_name, 'wb')
        txt = gzip.GzipFile(fileobj=txt_gz, mode='wt')
    else:
        txt = open(txt_name, 'w')

    # get the TMT intensities (if they exist)
    heights, tmt_intensity_dict = load_tmt_data(sqt_file, params)

    # get the raw file basename
    lc_name = os.path.split(sqt_file)[1]
    lc_name = os.path.splitext(lc_name)[0]
        
    # build and write header line
    header = ['start','end', 'Z', 'expM', 'SpRank', 'theoM', 'deltaCN',
              'Xcorr', 'Sequence', 'Loci', 'NewDeltaCN', 'ISBDisc',
              'NewDisc', 'ntt', 'ForR', 'Sp']
    if heights:
        header += heights
    txt.write('\t'.join(header)+'\n')
    
    # loop over scans and build output lines
    empty_scans = 0
    short_peptides = 0
    for out in outs:

        # double check peptide lookups
        for match in out.matches:
            for match_tuple in match.match_tuple_list:
                if len(match_tuple) == 0:
                    for obj in params.log_obj:
                        print('......WARNING: protein_matches not populated - Scan:', out.beg, out.matches[0].sequence, file=obj)

        # skip empty OUT files
        if len(out.matches) == 0:
            for obj in params.log_obj:
                print('......WARNING: empty OUT file encountered - Scan:', out.beg, file=obj)
            empty_scans += 1
            continue
        
        statline = [out.beg, out.end, out.z, out.MH1exp] # the "static" data        
        newdcn, olddcn = get_new_old_dcn(out.matches, out.top)
        max_xcorr = out.matches[0].xcorr
        
        for match in out.matches:
            
            # write lines for top hits only (match number equal to one)
            if match.match == 1:
                matchline = statline + [match.SpRank, match.MH1calc, str(olddcn), match.xcorr] # Match data
                try:
                    for match_tuple in match.match_tuple_list:                  
                        line = matchline + [match_tuple[3], match_tuple[0]]                  
                        length, ntt = PAW_lib.amino_acid_count(match_tuple[3], params.enzyme)
                        NewDisc = new_discrim_score(int(out.z), float(match.xcorr), length, newdcn,
                                                    int(match.SpRank), float(out.MH1exp), float(match.MH1calc))
                        line += [str(round(newdcn,5)), '', str(round(NewDisc,5)), str(ntt)]
                        if params.decoy_string in match_tuple[0]:
                            line += ['R']
                        else:
                            line += ['F']
                        line += [match.Sp]
                        if tmt_intensity_dict:
                            try:
                                line += tmt_intensity_dict[lc_name + '.' + out.beg]
                            except KeyError:
                                print('...WARNING: no key:', lc_name + '.' + out.beg)
                                line += ['0.0' for x in heights]
                        txt.write('\t'.join(line) + '\n')
                except IndexError:
                    short_peptides += 1
                    pass    # empty match tuple (short peptide)
    txt.close()
    if params.zipflag:
        txt_gz.close()

    # print out number of empty scans
    if empty_scans:
        for obj in params.log_obj:
            print('...there were %s/%s empty scans' % (empty_scans, len(outs)), file=obj)
    if short_peptides:
        for obj in params.log_obj:
            print('...there were %s short peptides' % (short_peptides,), file=obj)
        
    return    

def lookup_peptide_sequences(out_list, already_seen, params):
    """Function to lookup peptide sequences in protein sequences when peptides match to
    more than one protein.

    out_list is list of all OUT objects, many things in params object
    """
    for obj in params.log_obj:
        print('...starting lookups for %d results at %s' % (len(out_list), time.ctime()), file=obj)
    agree, disagree = 0, 0
    if params.enzyme_termini in [8, 9]:
        minimum_termini = 1
    else:
        minimum_termini = params.enzyme_termini

    out_count = 0
    for out in out_list:
        out_count += 1              
        if (out_count % 5000) == 0.0:              
            print('..%s scans looked up (%s)..' % (out_count, time.ctime()))
        try:
            for match in out.matches:
                length, ntt = PAW_lib.amino_acid_count(match.sequence)
                if length < params.min_pep_len or match.match != 1:     # skip short peptides, non-top hits
                    continue
                
                # try to lookup the peptide (faster)
                base_pep_masked = PAW_lib.get_base_peptide_sequence(match.sequence, params.mask)
                if base_pep_masked in already_seen:
                    match.match_tuple_list = already_seen[base_pep_masked]
                    continue

                if match.additional == 0:
                    prot_list = [params.proteins[params.prot_map[match.accession]]]

                else:
                    prot_list = []
                    lookup = True
                    if base_pep_masked in params.peptide_to_accessions:
                        # make protein list from digest dictionaries - look for full peptides first
                        for acc in params.peptide_to_accessions[base_pep_masked]:
                            prot_list.append(params.proteins[params.prot_map[acc]])
                        if len(prot_list) != match.additional + 1: # if we have all matches, skip semis
                            prot_list = params.proteins
                            lookup = False
                    if lookup:
                        # look for semi-trpytics next
                        acc_list = []
                        if base_pep_masked[:6] in params.Nterm_to_accessions:
                            for acc in params.Nterm_to_accessions[base_pep_masked[:6]]:
                                acc_list.append(acc)
                        if base_pep_masked[-6:] in params.Cterm_to_accessions:
                            for acc in params.Cterm_to_accessions[base_pep_masked[-6:]]:
                                acc_list.append(acc)
                        # simplify indexes of accession list (it might be redundant)
                        idx_list = sorted(set([params.prot_map[acc] for acc in acc_list]))
                        for idx in idx_list:
                            prot_list.append(params.proteins[idx])
                            
                    if not prot_list:
                        for obj in params.log_obj:
                            print('...WARNING: out count:', out_count, file=obj)
                            print('...WARNING: %s not in peptide dictionaries' % match.sequence, file=obj)
                        # if all else fails, search all proteins
                        prot_list = params.proteins

                # lookup matches, save in ordered dict, add matches to already_seen                        
                all_matches = PAW_lib.find_peptide(match.sequence, prot_list, params.mask, verbose=True)
                if not all_matches:
                    all_matches = PAW_lib.find_peptide(match.sequence, params.proteins, params.mask)

                # make sure matches are non-redundant and preserve order
                all_matches = list(OrderedDict([(key, None) for key in all_matches]).keys())
                
                """This is where premature stop codon peptides get rejected - modifiy amino_acid_count function?"""
                filtered_matches = [t for t in all_matches if PAW_lib.amino_acid_count(t[3])[1] >= minimum_termini]

                number_tryptic = len(filtered_matches)
                match.match_tuple_list = filtered_matches
                already_seen[base_pep_masked] = filtered_matches
                
                # see if number of protein matches agrees with Comet number
                if (match.additional+1) != number_tryptic:
                    print('Length proteins list:', len(prot_list))
                    if len(all_matches) != len(filtered_matches):
                        print('All matches:')
                        for _match in all_matches:
                            print(_match)
                    print('Filtered matches:')
                    for _match in filtered_matches:
                        print(_match)

                    disagree += 1
                    if length >= params.min_pep_len:
                        for obj in params.log_obj:
                            print('......WARNING: scan #%s, seq: %s , Comet had %s, lookup was %s' %
                                  (out.beg, match.sequence, match.additional+1, number_tryptic), file=obj)
                            if len(filtered_matches) < 5:                                                                            
                                print(filtered_matches, file=obj)
                else:
                    agree +=1
                            
        except IndexError:
            continue    # seems we get some empty scans with Comet
        
    for obj in params.log_obj:
        print('...match counts agreed: %s, disagreed: %s' % (agree, disagree), file=obj)
    return already_seen

def make_lookup_tables(params):
    """Reads a FASTA file and makes a dictionary of peptide sequence to
    protein accession list for theoretically digested proteins. Protein
    entries are retained and an accession to index map created.
    """
    # read the sequences from the FASTA file
    f = PAW_lib.FastaReader(params.data_base)
    p = PAW_lib.Protein()
    params.proteins = []
    while f.readNextProtein(p, check_for_errs=False):
        params.proteins.append(copy.deepcopy(p))
    start_digest = time.ctime()
    for obj in params.log_obj:
        print('...%s had %s proteins' % (os.path.basename(params.data_base),
                                             len(params.proteins)), file=obj)
        
    # make dictionary where key is accession string, value is protein index number
    params.prot_map = {}
    for i, protein in enumerate(params.proteins):
        params.prot_map[protein.accession] = i
    
    # ENZYME is one of "Tryp", "GluC", "AspN", or "LysC"
    if params.enzyme == "Tryp":
        regex = re.compile(r".(?:(?<![KR](?!P)).)*")    # trypsin strict
    elif params.enzyme == "GluC":
        regex = re.compile(r".(?:(?<![DE](?!P)).)*")    # Glu-C
    elif params.enzyme == "AspN":
        regex = re.compile(r".(?:(?![D]).)*")           # Asp-N
    elif params.enzyme == "LysC":
        regex = re.compile(r".(?:(?<![K](?!P)).)*")     # Lys-C strict
    else:
        for obj in params.log_obj:
            print('...WARNING: enzyme will default to trypsin', file=obj)
        regex = None
        
    # make the peptide lookup dictionary            
    params.peptide_to_accessions = {}
    params.Nterm_to_accessions = {}
    params.Cterm_to_accessions = {}
    
    for prot in params.proteins:
        # NOTE: need a large high peptide mass cutoff to get all semi-tryptics in dictionary below
        peptides = prot.enzymaticDigest(enzyme_regex=regex, high=1000000)    # use the digestion method of Protein class
        
        for pep in peptides:
            
            # see what to do with I and L residues
            if params.mask:
                mass_spec_seq = re.sub(r'[IL]', 'j', pep.seq)
            else:
                mass_spec_seq = pep.seq
                
            # add full peptides to their dictionary
            if mass_spec_seq in params.peptide_to_accessions:
                if prot.accession not in params.peptide_to_accessions[mass_spec_seq]:
                    params.peptide_to_accessions[mass_spec_seq].append(prot.accession)
            else:
                params.peptide_to_accessions[mass_spec_seq] = [prot.accession]

            # N-term tags (if doing semi-tryptic searches)
            if params.enzyme_termini == 1 or params.enzyme_termini == 8:
                if mass_spec_seq[:6] in params.Nterm_to_accessions:
                    if prot.accession not in params.Nterm_to_accessions[mass_spec_seq[:6]]:
                        params.Nterm_to_accessions[mass_spec_seq[:6]].append(prot.accession)
                else:
                    params.Nterm_to_accessions[mass_spec_seq[:6]] = [prot.accession]
                    
            # C-term tags (if doing semi-tryptic searches)
            if params.enzyme_termini == 1 or params.enzyme_termini == 9:
                if mass_spec_seq[-6:] in params.Cterm_to_accessions:
                    if prot.accession not in params.Cterm_to_accessions[mass_spec_seq[-6:]]:
                        params.Cterm_to_accessions[mass_spec_seq[-6:]].append(prot.accession)
                else:
                    params.Cterm_to_accessions[mass_spec_seq[-6:]] = [prot.accession]

    end_digest = time.ctime()
    for obj in params.log_obj:
        print('...peptide digest has %s peptides' % (len(params.peptide_to_accessions),), file=obj)
        if params.Nterm_to_accessions:
            print('...N-term lookup has %s tags' % (len(params.Nterm_to_accessions),), file=obj)
        if params.Nterm_to_accessions:
            print('...C-term lookup has %s tags' % (len(params.Cterm_to_accessions),), file=obj)
        print('...start: %s, end: %s' % (start_digest, end_digest), file=obj)
    
    return params

def process_one_match(match_number, static_match_data, match_to_dynamic):
    """Makes a list of SequestMatch objects from matches having the same Xcorr
    """
    match_list = []

    for key in match_to_dynamic:
        match = SequestMatch()
        match.match = match_number
        match.xcorr, match.deltaCN, match.XcRank, match.SpRank, match.Sp = static_match_data
        match.sequence, match.MH1calc, match.ionmatch, match.ionpredict, match.accession = match_to_dynamic[key][0]
        match.match_tuple_list = list(match.accession)
        if len(match_to_dynamic[key]) > 1:
            match.additional = len(match_to_dynamic[key]) - 1
        match_list.append(match)
    
    return match_list

def process_one_scan(buff, out):
    """Processes the equivalent of a single OUT file (one SQT scan block)
    """
    if not buff:
        return

    # initializations
    last_xcorr = 100.0
    match_number = 0
    static_match_data = []

    for line in buff:
        temp = line.split('\t')
        
        # match lines
        if temp[0] == 'M':
            
            # see if Xcorr is different
            if float(temp[5]) < last_xcorr:
                if static_match_data:
                    out.matches += process_one_match(match_number, static_match_data, match_to_dynamic)
                match_number += 1
                match_to_dynamic = {}
                static_match_data = []
                last_xcorr = float(temp[5])

            if not static_match_data:
                xcorr = temp[5]       # Xcorr score
                deltaCN = temp[4]     # deltaCN score (actually runs one behind in matches)
                XcRank = temp[1]      # Rank by Xcorr
                SpRank = temp[2]      # Rank by Sp score
                Sp = temp[6]          # Sp score (or expectation value)
                static_match_data = [xcorr, deltaCN, XcRank, SpRank, Sp] 
            
                sequence = temp[9]    # Peptide sequence (Comet format: bounding residues, special symbol PTMs)
                MH1calc = temp[3]     # Calculated MH+ mass
                ionmatch = temp[7]    # Number of matched ions
                ionpredict = temp[8]  # Number of predicted ions
                dynamic_match_data = [sequence, MH1calc, ionmatch, ionpredict]

            else:
                sequence = temp[9]    # Peptide sequence (Comet format: bounding residues, special symbol PTMs)
                MH1calc = temp[3]     # Calculated MH+ mass
                ionmatch = temp[7]    # Number of matched ions
                ionpredict = temp[8]  # Number of predicted ions
                dynamic_match_data = [sequence, MH1calc, ionmatch, ionpredict]

        # loci lines
        elif temp[0] == 'L':
            dynamic_match_data.append(temp[1]) # append accession to dynamic match list
            try:
                number_matches = int(temp[2]) + 1
            except IndexError:
                number_matches = 1
            key = sequence.split('.')[1]
            for dummy in range(number_matches):
                if key in match_to_dynamic:
                    match_to_dynamic[key].append(dynamic_match_data)
                else:
                    match_to_dynamic[key] = [dynamic_match_data]

    out.matches += process_one_match(match_number, static_match_data, match_to_dynamic) # need to get last one

def convert_one_sqt_file(sqt_file, params):
    """Parses SQT file and makes a list of OutData objects.
    """
    # check file type and open accordingly
    if params.zipflag:
        sqt_file_obj = gzip.open(sqt_file, 'rb').readlines()
    else:
        sqt_file_obj = open(sqt_file,'r').readlines()

    # initializations
    out = OutData()        
    outs = []
    buff = []

    # read scan, match, and loci lines. Skip everything else    
    for line in sqt_file_obj:
        line = line.rstrip()
        # skip header lines
        if line.startswith('H'):
            continue

        # scan lines
        if line.startswith('S'):

            # process last scan
            process_one_scan(buff, out)
            buff = []
            
            # append "last" scan to "outs" (runs one behind) and create new object for next scan
            outs.append(out)
            out = OutData()

            # get static data from S line
            temp = line.split('\t')
            out.beg, out.end, out.z, out.MH1exp = temp[1], temp[2], temp[3], temp[6]
            out.top = params.output_lines

        # save M and L lines in buff
        else:
            buff.append(line)
            
    # process the last scan
    process_one_scan(buff, out)
    outs.append(out)
    
    return outs[1:] # the first object is blank, so skip

def set_up_conversions(sqt_file, params):
    """Reads header lines from SQT file and sets some global parameters.
    Assumes all SQT files in a folder have same set of header lines,
    i.e. all from same Comet search.

    Does an in silico digest of the FASTA file for peptide lookups. Adds
    semi-tryptic N- and C-term 6 aa tags if semi-tryptic searches were used.
    """        
    # get path to location of SQT files and open log file there
    path = os.path.dirname(sqt_file)
    params.default_location = path
    log_name = 'sqt_conversion_log.txt'
##    log_name = 'sqt_conversion_%s_log.txt' % time.time() # this add a time stamp - will have multiple log files
    log_obj = open(os.path.join(path, log_name), 'wt')
    params.log_obj = [None, log_obj]

    # we have log file setup so print some info
    for obj in params.log_obj:
        print('==========================================================', file=obj)
        print(' sqt_converter.py, %s, written by Phil Wilmarth, 2019  ' % VERSION, file=obj)
        print('==========================================================', file=obj)
        print('processing: %s: %s' % (os.path.dirname(sqt_file), time.ctime()), file=obj)
    
    # check file type and open accordingly
    if params.zipflag:
        sqt_file_obj = gzip.open(sqt_file, 'rb').readlines()
    else:
        sqt_file_obj = open(sqt_file,'r').readlines()

    # read header lines, then close
    for line in sqt_file_obj:
        line = line.strip()
        
        if line.startswith('H\tCometParams'): # CometParams header lines
            parameters = line.split('\t')
            if len(parameters) <= 2:
                continue
            if parameters[2].startswith('database_name'):
                params.data_base = parameters[2].split('= ')[1]
                if not os.path.exists(params.data_base):    # browse to database if path not found
                    ext_list = [('FASTA files', '*.fasta'), ('Zipped files', '*.gz'), ('All files', '*.*')]
                    title = 'Please locate %s' % os.path.basename(params.data_base)
                    params.data_base = PAW_lib.get_file(params.default_location, ext_list, title)
                    if not params.data_base: sys.exit()  # cancel button response is exit()
            if parameters[2].startswith('num_output_lines'):
                params.output_lines = int(parameters[2].split('= ')[1]) # number of lines (M/L pairs) written to sqt file
            if parameters[2].startswith('num_enzyme_termini'):
                params.enzyme_termini = int(parameters[2].split('= ')[1]) # minimum number enzymatic termini
            if parameters[2].startswith('allowed_missed_cleavage'):
                params.missed_cleavages = int(parameters[2].split('= ')[1]) # max number missed cleavages
            if parameters[2].startswith('print_expect_score'): # whether Sp score is that (0) or expect value (1)
                params.expect_score = int(parameters[2].split('= ')[1])
            if parameters[2].startswith('search_enzyme_number'): # enzyme number used in search
                params.enzyme_number = int(parameters[2].split('= ')[1])
                """Need to add support for all enzymes"""
                if params.enzyme_number in [0, 1]: # trypsin for No_enzyme, too
                    params.enzyme = 'Tryp'
                elif params.enzyme_number == 3:
                    params.enzyme = 'LysC'
                elif params.enzyme_number == 6:
                    params.enzyme = 'AspN'
                elif params.enzyme_number == 8:
                    params.enzyme = 'GluC'
                else:
                    for obj in params.log_obj:
                        print('...WARNING: enzyme not supported!', file=obj)
                    params.enzyme = None
                    
            """Should get static mods from CometParams and build amino acid table for digestion routine
            """
        if line.startswith('S\t'):
            break
        
    # read the FASTA file and make lookup tables
    params = make_lookup_tables(params)
            
    return params

def load_tmt_data(sqt_file, params):
    """Loads TMT reporter ion data (if present) from RAW files folder."""
    # get container folder path
    sqt_folder, sqt_filename = os.path.split(sqt_file)
    if sqt_filename.endswith('.sqt.gz'):
        tmt_filename = sqt_filename[:-7] + '.PAW_tmt.txt'
    elif sqt_filename.endswith('.sqt'):
        tmt_filename = sqt_filename[:-4] + '.PAW_tmt.txt'
    else:
        for obj in params.log_obj:
            print('...WARNING: sqt filename does not end in sqt', file=obj)
        
    project_folder = os.path.dirname(sqt_folder)

    # walk project folder and look for ".PAW_tmt.txt" files
    found = False
    for root, dirs, files in os.walk(project_folder):
        if tmt_filename in files:
            tmt_filename = os.path.join(root, tmt_filename)
            found = True
            for obj in params.log_obj:
                print('...found:', tmt_filename, file=obj)

    # warning if no files found
    if not found:
        for obj in params.log_obj:
            print('...WARNING: the PAW_tmt.txt file was not found!', file=obj)
        return None, None

    # read and parse the TMT information
    heights = []
    tmt_intensity_dict = {}
    for obj in params.log_obj:
        print('...getting tmt data from:', os.path.basename(tmt_filename), file=obj)
    bad_keys = set()
    with open(tmt_filename, 'rt') as tmt_obj:
        for line in tmt_obj:
            line = line.rstrip()
            if line.startswith('lc_name\tms2_scan'):    # find the header line (assumes first line)
                old_heights = heights
                col_map = {v: i for i, v in enumerate(line.split('\t'))}
                heights = [x for x in line.split('\t') if x.startswith('height_')]
                if old_heights and (heights != old_heights):
                    for obj in params.log_obj:
                        print('...WARNING: TMT reporter ions are not consistent!', file=obj)
                continue
            
            # this does the rest of the lines
            intensities = []
            items = line.split('\t')
            key = items[col_map['lc_name']] + '.' + items[col_map['ms2_scan']] # key is LCName.MS2Scan
            for height in heights:
                intensities.append(items[col_map[height]])
                
            # update dictionary
            if key in tmt_intensity_dict:
                for obj in params.log_obj:
                    print('...WARNING: non unique key: (will be zeroed out)', key, file=obj)
                bad_keys.add(key)
            else:
                tmt_intensity_dict[key] = intensities
    for key in bad_keys:
        tmt_intensity_dict[key] = ['0.0' for x in heights]

    return heights, tmt_intensity_dict        
    
def set_default_params(params):
    """Sets up some default values in params object.
    """
    params.min_pep_len = MIN_PEP_LEN - 1 # one less than actual desired minimum peptide length
    if params.min_pep_len < 6: # but not less than 6
        params.min_pep_len = 6
    params.mask = MASK  # set to global value (from top)
    params.decoy_string = DECOY_STRING
    params.default_location = DEFAULT_LOCATION
    return

def convert_sqt_files(sqt_list, params):
    """Converts all SQT files in list into PAW text files.
    """    
    # read header lines from one sqt file, set up parameters, make lookup tables
    params = set_up_conversions(sqt_list[0], params)

    # loop over SQT files and unpack data into OutData objects
    scan_total = 0
    already_seen = {} # save looked up peptides
    for sqt_file in sqt_list:
        outs = convert_one_sqt_file(sqt_file, params)
        scan_total += len(outs)
        
        # print status line after file parsing has finished
        for obj in params.log_obj:
            print('\n...done reading: %s (%s scans)' % (os.path.basename(sqt_file), len(outs)), file=obj)

        # lookup peptides that match to more than one protein, and get match tuples for all                
        lookup_peptide_sequences(outs, already_seen, params)
            
        # we have all of the data in OutData objects so write the TXT file
        make_PAW_txt_file(sqt_file, outs, params)
        
    for obj in params.log_obj:
        print('\nConversion of %s scans has completed!\n' % scan_total, file=obj)

    # try and close log files
    for obj in params.log_obj:
        try:
            obj.close()
        except:
            pass
        
    return

#############################
# MAIN program starts here:
#############################

def main(folder, overwrite=True):
    """ Converts Comet search results (SQT format) to PAW TXT files.    
    written by Phil Wilmarth, David Lab, OHSU, 2016.    
    """    
    # check for valid path
    if not folder:
        print('......WARNING: no file folder was selected!')
        return
    if not os.path.exists(folder):
        print('...WARNING: folder could not be found')
        return

    # create container for parameters and lookup tables (easier to pass around)
    params = PAW_lib.CometParams()  # this might be mostly a namespace mechanism
    set_default_params(params)
   
    # see if there are SQT files or compressed files
    params.zipflag = False
    os.chdir(folder)
    sqt_list = glob.glob('*.sqt')
    txt_list = [x for x in sqt_list if os.path.exists(os.path.join(folder, x.replace('.sqt', '.txt')))]
    gzip_sqt_list = glob.glob('*.sqt.gz')
    gzip_txt_list = [x for x in gzip_sqt_list if os.path.exists(os.path.join(folder, x.replace('.sqt.gz', '.txt.gx')))]
    if len(gzip_sqt_list) > len(sqt_list):
        sqt_list = gzip_sqt_list
        txt_list = gzip_txt_list
        params.zipflag = True

    # skip folder if no SQT files were found
    if len(sqt_list) == 0:
        print('......WARNING - folder does not contain any SQT files.')
        return

    # skip folder if TXT files already exist and overwrite flag is false
    if (len(txt_list) == len(sqt_list)) and not overwrite:
        print('......WARNING - TXT files already exist so skipping folder.')
        return

    # convert all of the SQT files
    sqt_list = [os.path.join(folder, x) for x in sqt_list]
    convert_sqt_files(sqt_list, params)
    return

if __name__ == '__main__':
    try:
        if os.path.exists(sys.argv[1]):
            sqt_folder = sys.argv[1]
            
            # stand alone execution will overwrite TXT files (safer)
            main(sqt_folder, overwrite=True)            
        else:
            print('FATAL: invalid file path')
            sys.exit()
            
    except IndexError:

        # no arguments passed so browse and select SQT files
        ext_list = [('SQT files', '*.sqt'), ('Compressed SQT files', '*.sqt.gz')]
        title = 'Select one or more SQT files'
        sqt_list = PAW_lib.get_files(DEFAULT_LOCATION, ext_list, title)
        if not sqt_list: sys.exit()    # cancel button response

        # see if compressed files or not
        params = PAW_lib.CometParams()
        set_default_params(params)
        sqt_count = len([x for x in sqt_list if x.endswith('.sqt')])
        gz_count = len([x for x in sqt_list if x.endswith('.sqt.gz')])
        if gz_count > sqt_count:
            params.zipflag = True
        else:
            params.zipflag = False

        # convert all of the SQT files
        convert_sqt_files(sqt_list, params)

# end
        
    


            
