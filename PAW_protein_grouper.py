"""program "PAW_protein_grouper.py"
Produces quantitative protein and peptide results tables.

This is a protein grouping step in the PAW analysis pipeline for Comet
searches of large data sets where spectral count or TMT quantitation is being
performed. The program combines similar proteins into "families" (homology
grouping), redetermines shared and unique peptide status, and re-splits
shared peptide counts.

Written by Phil Wilmarth, 2011-2017, OHSU.

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

# global imports statements
import os
import sys
#
VERSION = 'v1.0.4'  # updated 8/30/2017

# updated for Python 3, 8/29/17 -PW

"""More extensive documentation added April, 2016 -PW

Protein and peptide summary files from the PAW pipeline are read in for one or more
biological samples in a proteomics experiment.

Grouping involves testing proteins for similarity. Comparisons are always done between pairs of
proteins. Peptide sequence sets are constructed for each protein from its associated peptides.
Peptides sequences are "masp spec" sequences where "I" and "L" redidues are replaced with a "j".
When comparing two proteins' peptide sets, three sets are determined: the intersection (the
set of peptides the two proteins share), and the left and right difference sets (the peptides
that are unique to each protein, respectively). A quantitative value is computed for each set.
These values are the experiment-wide total spectral count of all the peptides contained within
each set.

Extenstions to the basic parsimony concepts of redundant proteins and subset proteins are tested
first. Instead of exaclty identical peptide sets for redundant proteins and true formal set
subsets for subsets, we define "nearly" identical and "almost" a subset using the "pseudo"
parameter. This is based on the long standing minimum protein identification requirement of
two peptides per protein.
"""
# variable parameters to control grouping
VERBOSE = True      # controls how much information is printed to console (global parameter)
PSEUDO = 2.0        # how much independent total count evidence is needed for pseudo-redundants, -subsets
LOW = 2.5           # true unique total count evidence threshold used in combination test
SHARE = 40.0        # absolute total shared count threshold used in combination test
MASS = 20.0         # relative fraction of shared-to-unique counts threshold used in combination test
#
IL = True           # make I and L indistinguisable (True)
diffZ = True        # treat different charge states as distinct (True)


def get_file(default_location, extension_list, title_string=""):
    """Dialog box to browse to a file.  Returns full file name.

    Usage: full_file_name = get_file(default_location, extension_list, [title]),
        where "default_location" is a starting folder location,
        extension_list is a list of (label, pattern) tuples,
        e.g. extension_list = [('Text files', '*.txt')],
        "title" is an optional message to list in the dialog box, and
        "full_file_name" is the complete path name of the selected file.

    Written by Phil Wilmarth, OHSU, 2008.
    """
    import tkinter
    from tkinter import filedialog

    # set up GUI elements
    root = tkinter.Tk()
    root.withdraw()

    # set default title string if not passed
    if title_string == "":
        title_string = 'Please select a single FILE'

    # create dialog box and return file selection
    root.update()
    return filedialog.askopenfilename(parent=root, initialdir=default_location,
                                      filetypes=extension_list, title=title_string)
    # end


def time_stamp_logfile(message, file_obj):
    """Prints "message" and time stamp to a log file.

    Written by Phil Wilmarth, OHSU, 2009.
    """
    import time
    print('%s on %s' % (message, time.ctime()), file=file_obj)
    return


class ProteinTuplet:
    """Object to hold protein information from the protein_summary
    file rows. This is a simple data container object. No methods.
    All data from each row are collected so results table can be
    very similar to the input table.

    Redundant proteins occupy multuiple rows in the standard PAW
    summary file so that acceessions and descriptions are easier
    to read. Here we need single protein objects that correspond
    to their respective sets of peptides. We will not create separate
    protein objects for redundant proteins. We will make an object for
    the first member and add some information about any additional
    redundant members to that object. This data structure will also
    hold information about the larger protein families that result
    from the protein grouping tests. Since "protein groups" already
    have a defined meaning in proteomics, we will use the term
    "families". It is really more specific than that since we are
    really talking about very similar family members such as triplets.
    That is why the object is called a ProteinTuplet instead of a
    ProteinGroup.

    written by Phil Wilmarth, 2011, OHSU.
    """

    def __init__(self):
        """Basic constructor, no parameters. Attributes are
        listed below along with descriptions.
        """
        self.number = '0'           # ProtGroup (group number string)
        self.counter = ''           # Counter column of ones
        self.name = ''              # Accession (or derived from)
        self.redundants = []           # List of other group/redundant member accessions
        self.filter = ''            # Filter (flag column)
        self.coverge = ''           # Coverage (sequence coverage)
        self.length = ''            # SeqLength (sequence length)
        self.mw = ''                # MW (calculated MW)
        self.description = ''       # Description (protein description string)
        self.tot_count = ''         # CountTot column
        self.tot_unique = ''        # UniqueTot column
        self.uniqfrac = ''          # UniqFrac column
        self.totals = []            # total count list (strings)
        self.uniques = []           # unique count list (strings)
        self.correcteds = []        # corrected count list (strings)
        self.otherloci = ''         # OtherLoci (list of accessions)
        self.all_loci = set()       # expanded set of other loci
        self.new_totcounts = []     # new total count list (list is number of samples long)
        self.new_uniqcounts = []    # new unique count list
        self.new_corrcounts = []    # new corrected count list
        self.done = False           # set to True if OtherLoci combining has been completed
        self.keep = True            # set to False if protein has been combined with another
        self.combined = []          # constituent proteins for combinations

        return
    # end ProteinTuplet class


class ProteinSummary:
    """Container for proteins and protein tuplets (indistinguishable by peptides).
    In addition to the traditional definition of redundant proteins, protein
    tuplets may include pseudo-redundant protein groupings, pseudo-subset protein
    groupings, and highly similar protein groupings.

    Several summary mappings (dictionaries) are also created to facilitate
    processing.

    Methods:
        load_results: reads in a protein results file (with redundant protein
            rollup).
        get_samples: parses the sample names from column headers.
        parse_row: parses the fields from each row of the results file.
        fix_otherloci: updates OtherLoci to include all proteins having any shared peptides.
        print_summary: prints summary table to new results file

    written by Phil Wilmarth, OHSU, 2011
    """

    def __init__(self):
        """Basic constructor, no parameters. Attributes and definitions are listed.
        """
        self.samples = []   # list of sample names (from Total count column headings)
        self.tuplets = []   # list of protein results objects (ProteinTuplets)
        self.col_map = {}   # map of column headings to column position
        self.header = ''    # save a copy of the header line
        self.acc_map = {}   # key is accession (str), value is tuplet index (int)
        self.acc_to_desc = {} # dictionary of accessions and associated descriptions (for reports)
        self.tuplet_map = {} # key is group number (str), value is accession (str)
        self.row_count = 0  # count of number of protein rows read in
        self.protein_count = 0  # number of protein tuplets
        self.distinguishable_count = 0   # number of proteins having some distinguisable peptides
        self.indistinguishable_count = 0    # number of protein groups having indistinguishable peptides
        return

    def load_results(self, file_name):
        """Read protein results file and parse main table row fields
        into ProteinTuplet objects. Protein_summary files have some lines
        before the main table and some short lines after the main table.

        Lines are skipped until the header line is found. The length of the
        header line (after splitting) is used to determine what lines are
        in the main table body (the last column is optional in the table so
        we have to adjust the length). The header line is split into its
        column headers and the column header text is used to index the column
        position using a "col_map" dictionary (saved as an attribute).
        """
        try:
            open_file = open(file_name, 'r')
        except IOError:
            print('...WARNING: could not open protein results file.')
        #
        columns = 10000 # so that lines before the header get skipped
        for line in open_file:
            line = '\t'.join([x.strip() if x.strip() else ' ' for x in line.split('\t')]) # remove EOL and clean up headers
            col = line.split('\t')    # split into column headers
            if col[0] == 'ProtGroup' and col[1] == 'Counter': # look for header line
                self.header = line  # save a copy of the header line
                columns = len(col) - 2 # minimum length of actual protein rows
                for i, col_header in enumerate(col):
                    self.col_map[col_header] = i # map column headings to positions
                self.get_samples(col) # get the sample names
            elif len(col) < columns: # skip any short lines (before or after the protein table)
                continue
            else:
                self.parse_row(col) # parse the protein row
            #
        open_file.close()
        """We "walk" the OtherLoci proteins to determine the full set of all
        proteins that are linked via any shared peptides. Proteins removed as
        subsets during previous PAW processing may still be present in OtherLoci
        lists but missing from the results table. We remove any references to those.
        """
        self.fix_otherloci() # create full sets of OtherLoci and remove any left overs
        self.protein_count = len(self.tuplets)
        for prot_tuplet in self.tuplets:
            if len(prot_tuplet.redundants) > 1:
                self.indistinguishable_count += 1
            else:
                self.distinguishable_count += 1
        return

    def get_samples(self, col):
        """Parses sample names from header line.

        Counts from the samples are in three repeating blocks of columns:
        Total counts, unique counts, and corrected counts (fractional
        counting of shared peptides, split based on other protein information).
        The columns start after the "UniqFrac" column and end before the
        "OtherLoci" column. The type of counts are labeled in the row above
        the main header row. The column headers in the main header row are
        just the sample names.
        """

        # get the sample names
        start = self.col_map['UniqFrac'] + 1
        end = int(start + (self.col_map['OtherLoci'] - start) / 3)
        for sample in col[start:end]:
            self.samples.append(sample)
        return

    def parse_row(self, col):
        """Parses protein results rows with redundant protein rollup.
        col is the list result after splitting the row line on tab character.

        "col_map" is the dictionary to get index into "col" from column header
        text string (more readible and independent of column order). As we combine
        redundant proteins and/or create groups, we keep one protein description
        and make lists of additional protein accessions. We will make a bigger
        dictionary to remember all the protein descriptions so we can make reports
        at the end. We add a new data object when we get a new group number,
        otherwise we "rollup" additional redundant protein information into previous
        data object.
        """

        # add accession (key) and description (value) to acc_to_desc dictionary
        self.acc_to_desc[col[self.col_map['Accession']]] = col[self.col_map['Description']]

        # see if group number is new (not in dictionary)
        group = col[self.col_map['ProtGroup']].split('.')[0] # NOTE: group is a string not an integer
        if group not in self.tuplet_map:
            # add dictionary entry and instantiate new ProteinTuplet object (for group)
            self.tuplet_map[group] = col[self.col_map['Accession']]
            self.acc_map[col[self.col_map['Accession']]] = len(self.tuplets) # add accession to dictionary
            new = ProteinTuplet()

            # set tuplet object attributes
            new.number = group  # NOTE: a string not an integer
            new.counter = col[self.col_map['Counter']]
            new.name = col[self.col_map['Accession']]
            new.redundants.append(new.name) # add self to list of accessions
            new.filter = col[self.col_map['Filter']]
            new.coverage = col[self.col_map['Coverage']]
            new.length = col[self.col_map['SeqLength']]
            new.mw = col[self.col_map['MW']]
            new.description = col[self.col_map['Description']]
            new.tot_count = col[self.col_map['CountsTot']]
            new.tot_unique = col[self.col_map['UniqueTot']]
            new.uniqfrac = col[self.col_map['UniqFrac']]

            # save spectral count numbers (no need to fix any non-grouped proteins)
            start = self.col_map['UniqFrac'] + 1
            samples = len(self.samples)
            new.totals = col[start: start+samples]  # saved as strings
            new.uniques = col[start+samples: start+2*samples]   # saved as strings
            new.correcteds = col[start+2*samples: start+3*samples]  # saved as strings
            try:
                new.otherloci = col[self.col_map['OtherLoci']] # rows w/o OtherLoci shorter so catch index error
                if new.otherloci == ' ':
                    new.otherloci = ''
            except IndexError:
                pass

            # add to list of tuplets (groups)
            self.tuplets.append(new)

        else:
            # redundant row so update name and add to the redundant members list of last tuplet
            tuplet = self.tuplets[-1]
            i = len(tuplet.redundants)     # number of additional redundant proteins so far
            if i > 1:
                last_suffix = '(+%s)' % (i-1,)
                next_suffix = '(+%s)' % (i,)
                tuplet.name = tuplet.name.replace(last_suffix, next_suffix) # save name with updated number
            else:
                tuplet.name = tuplet.name + ' (+1)'
            tuplet.redundants.append(col[self.col_map['Accession']]) # add to accession list

        self.row_count += 1
        return

    def fix_otherloci(self):
        """Fixes accessions in OtherLoci that are not in results file.
        Also builds more complete OtherLoci sets.

        If A shares some peptides with B and B shares some peptides with both A and C,
        then OtherLoci for A will be A, B. OtherLoci for B will be A, B, C. For group
        testing we need to consider A, B, and C. So we want the OtherLoci for A to be
        A, B, C. We visit all other proteins pointed to by OtherLoci and add loci. We
        keep following OtherLoci until the set stops growing.
        """

        # get rid of loci that did not make the protein report (because of Parsimony filter)
        for i, prot in enumerate(self.tuplets):
            if prot.otherloci == '': # skip unique proteins
                continue
            acc_list = [x.strip() for x in prot.otherloci.split(',')] # splits and removes whitespace
            new_acc_list = []
            for acc in acc_list:
                acc = acc.replace('"','')   # could have quotes enclosing cells (CSV-ish from Excel)
                if acc in self.acc_map: # skip accessions that do not correspond to rows in results file
                    new_acc_list.append(acc)
            if len(new_acc_list) > 1: # update OtherLoci strings with only relevant accessions
                prot.otherloci = ', '.join(new_acc_list)
            else:
                prot.otherloci = ''

        # build more comprehensive OtherLoci-like lists
        for i, prot in enumerate(self.tuplets):
            if prot.otherloci:
                new = set([x.strip() for x in prot.otherloci.split(',') if x != '']) # starting set of other loci
                previous = len(new) # starting total number of other loci
                current = 0
                while current != previous: # stop when set stops growing
                    old = set(new) # reset old for next iteration (set call makes copy)
                    previous = len(old)
                    for loci in old: # make union of all OtherLoci sets for loci in old
                        otherloci = set([x.strip() for x in
                                         self.tuplets[self.acc_map[loci]].otherloci.split(',') if x != ''])
                        new.update(otherloci)
                    current = len(new) # length of (possibly) more complete set
                prot.all_loci = set(new) # save the more complete set of other loci
        return

    def print_summary(self, fileobj=None, log_list=None):
        """Prints a new protein summary file.

        We create a similar format to the original protein_summary file.
        OtherLoci and grouping information gets moved to the beginning of the
        rows instead of at the ends. Lists of accessions get concatenated into
        strings.
        """
        import time

        # header stuff first
        print('Program "PAW_protein_grouper.py" run on', time.asctime(), file=fileobj)
        print(file=fileobj)
        print(file=fileobj)
        header = self.header.split('\t')
        header.insert(self.col_map['Filter'], header[-1]) # add OtherLoci to before "Filter" col
        header = header[:-1] # remove original OtherLoci header at end
        header.insert(self.col_map['Filter'], 'Similar')  # column for family groups
        header.insert(self.col_map['Filter'], 'Identical')    # column for fully-redundant proteins
        pre_header = ['' for x in header]   # row above main header
        rows = len([x for x in self.tuplets if x.keep == True])
        prefix = 5
        pre_header[1] = '=SUBTOTAL(109,B%s:B%s)' % (prefix+1, rows+prefix)  # make counter subtotal function
        print('\t'.join(pre_header), file=fileobj)
        print('\t'.join(header), file=fileobj)

        # protein rows next (skip any rows not marked as keepers; these were added to family groups)
        for p in self.tuplets:
            if p.keep:

                # most quantities from original protein_summary were saved as strings
                row = [p.number, p.counter, p.name]
                if len(p.redundants) > 1:   # Identical (redundants) column
                    row.append('&'.join(p.redundants))
                else:
                    row.append(' ')
                if len(p.combined) > 0:     # Similar (grouped families) column
                    row.append('&'.join(p.combined))
                else:
                    row.append(' ')
                if len(p.otherloci) > 0: # OtherLoci column (original not expanded?)
                    row.append(p.otherloci)
                else:
                    row.append(' ')
                row += [p.filter]
                try:
                    row.append('%0.1f' % (float(p.coverage),)) # do some formatting if possible
                except ValueError:
                    row.append(p.coverage)
                row += [p.length, p.mw, p.description, p.tot_count, p.tot_unique]
                row.append('%0.3f' % (float(p.uniqfrac),)) # do some formatting
                row = row + p.totals + p.uniques + p.correcteds
                print('\t'.join(row), file=fileobj)

    # end ProteinSummary class


class PeptideSummary:
    """Container for Peptide Summary file contents.

    Methods:
        __init__: basic constructor, defined attributes
        load_results: opens and parses peptide summary file
        parse_row: parses peptide rows into IdentifiedPeptides objects
        print_summary: prints a new peptide summary table


    Each protein has one or more rows of peptide information in the
    peptide_summary file. We will collect all rows for a given protein
    accession as an IdentifiedPeptides object. We will store those
    objects in a list, one item for each protein group/tuplet accession.
    We will store the tuple indices in a dictionary keyed by accessions.

    We will save a copy of the header line and make a dictionary to
    access columns by column header text strings.

    written by Phil Wilmarth, OHSU, 2011.
    """
    def __init__(self, IL=True, diffZ=True):
        """Simple constructor, no parameters. Attributes and definitions:
        """
        self.tuplets = []   # list of IdentifiedPeptides objects, one per protein tuplet
        self.col_map = {}   # maps column headings to column positions
        self.acc_map = {}   # key is accession (str), value is tuplets index (int)
        self.line_count = 0 # counts number of lines read
        self.header = ''    # saves a copy of the header line
        self.samples = []   # list of sample names
        self.IL = IL        # True makes I and L indistinguishable
        self.diffZ = diffZ  # True keeps different charge states as distinct sequences
        return

    def load_results(self, file_name):
        """Reads a peptide summary file to list and indexes accessions.
        We allow for optional "short" lines before and after the main table,
        whose width is determined by the length of the split header row.
        """
        try:
            open_file = open(file_name, 'r')
        except IOError:
            print('...WARNING: could not open peptide results file.')
        #
        columns = 100000 # so that lines before the column headers get skipped
        for line in open_file:
            line = '\t'.join([x.strip() if x.strip() else ' ' for x in line.split('\t')])    # remove end-of-line characters
            col = line.split('\t')  # split into pieces
            if col[0] == 'ProtGroup' and col[-1] == 'OtherLoci':    # header line
                self.header = line  # so save
                columns = len(col) - 1  # the length of data rows in the table (OtherLoci may be empty)
                for i, col_header in enumerate(col):    # map column headings to column positions
                    self.col_map[col_header] = i
                beg = self.col_map['TotCount']+1
                end = self.col_map['OtherLoci']
                self.samples = col[beg:end] # get the sample names (between the above positions)
            elif len(col) < columns: # skip any lines without data
                continue
            else:
                self.parse_row(col, line)
        return

    def parse_row(self, col, line):
        """Saves peptide results rows into an IdentifiedPeptides Object (one object per protein).
        """
        if col[self.col_map['Accession']] not in self.acc_map:   # if new protein, add new object to tuplets
            self.tuplets.append(IdentifiedPeptides(self.IL, self.diffZ))
            self.tuplets[-1].group = col[self.col_map['ProtGroup']]
            self.tuplets[-1].acc = col[self.col_map['Accession']]
            self.acc_map[col[self.col_map['Accession']]] = len(self.tuplets)-1 # map accession to tuplet index
        self.tuplets[-1].rows.append(line) # always append line to the last tuple's row list
        self.line_count += 1
        return

    def print_summary(self, fileobj=None, log_list=None):
        """Prints a new peptide summary file. Any proteins that have been combined are
        flagged as non-keepers and skipped. Combined proteins add new tuplets with combined
        peptide rows.
        Write the header information, then the peptide rows, then the parameter settings.
        """
        import time

        print('Program "PAW_protein_grouper.py" run on', time.asctime(), file=fileobj)
        print('\n\n', file=fileobj)
        print(self.header, file=fileobj)
        for prot in self.tuplets:
            if prot.keep:
                for row in prot.rows:
                    print(row, file=fileobj)
        return

    # end PeptideSummary class


class IdentifiedPeptides:
    """Raw container for each protein's identified peptides.
    No significant parsing is done at initialization.

    Methods:
        get_peptide_sets_counts: makes peptide sequence sets and
            dictionaries of counts associated with each peptide.

    written by Phil Wilmarth, OHSU, 2011.
    """
    def __init__(self, IL=True, diffZ=True):
        """Basic constructor, no parameters.
        """
        self.group = ''     # protein group number
        self.acc = ''       # protein accession string
        self.keep = True    # set to False if protein has been combined with another
        self.rows = []      # rows of identified peptide information from the peptide summary file
        self.IL = IL        # flag to set I and L indistinguishable
        self.diffZ = diffZ  # flag to keep charge states distinguishable

    def get_peptide_sets_counts(self, col_map, samples):
        """Returns a protein's master peptide set for the entire experiment.
        IL flag replaces I/L with j, diffZ flag makes different peptide
        charge states distinct.

        Also returns a list of dictionaries of peptide:count for all
        individual samples and the experiment-wide grand total.
        """
        import re

        pep_set = set()
        pep_count = [{} for s in samples]
        pep_count.append({}) # grand totals
        sample_list = samples + ['TotCount']

        for row in self.rows:
            col = [x.strip() for x in row.split('\t')]

            # get the sequence
            seq = col[col_map['Sequence']].split('.')[1]    # get peptide sequence without bounding residues
            if self.diffZ:
                seq = seq + '.' + col[col_map['Z']]     # add charge state to end of peptide sequence
            if self.IL:
                seq = re.sub(r'[IL]', 'j', seq)     # mask I and L with j (make peptide sequence mass spec-like)
            pep_set.add(seq)    # add to peptide set (only need

            # get the counts
            for i, sample in enumerate(sample_list):
                if col[col_map[sample]] == '':
                    count = 0
                else:
                    try:
                        count = int(col[col_map[sample]])
                    except ValueError:
                        count = float(col[col_map[sample]])
                pep_count[i][seq] = count   # save count in dictionary

        return pep_set, pep_count

    def _snoop(self):
        print('============ identified peptide contents ===============')
        print('group number:', self.group)
        print('accession:', self.acc)
        print('keep flag:', self.keep)
        print('number of rows:', len(self.rows))

    # end IdentifiedPeptides class


class OtherLociCombiner:
    """Tests if any of the proteins having shared peptides are
        similar enough that they should be combined.

        Performs three tests: two are extenstions of redundant protein
        grouping and peptide subset removal, respectively. The third
        looks for cases where the collective shared peptide counts between
        two protein are large compared to any unique peptide counts. The
        test values are total spectral counts of peptides not smaller peptide
        sequence counts. Noise peptide matches will tend to produce small
        counts, real peptide matches can be larger, particularly for
        abundant proteins. Total counts should separate signal from noise
        better than non-redudant peptide sequence counts.

        Some other key points:
        There is extensive use of sets and count dictionaries to keep track
        of everything and calculate the quantities needed for the testing.

        Only sets of proteins that share peptides with each other are tested.
        There can be many such sets in typical experiments. The sets will
        range in size but will typically be relatively small. Sizes will depend
        on the proteins in the sample and the nature of the protein database.

        Comparisons are all pair-wise and shared and unique peptides are often
        defined with respect to the pair being compared.

        Comparisons are ordered from proteins with the greatest unique counts to
        lowest. This helps prevent over grouping due to lack of information.

        Comparison tests are iterative. Test logic is driven by dynamic
        definitions of shared and unique peptides, which change during
        grouping. After any pair of proteins are grouped, everything is
        recomputed and testing starts over. Testing continues until no
        more proteins are grouped.

        Protein grouping usually changes definitions of shared and
        unique peptide status (we define things with respect to the set
        of proteins identified in the sample rather than with respect
        to the entire protein database). That changes how shared protein
        spectral counts are split between proteins. Several things need to
        be updated in the protein and peptide reports. There is quite a
        bit of code to deal with these aspects.

        Methods:
        __init__: basic consructor
        load_data: builds peptide sets and count dictionaries (for each protein, for each sample)
        peptide_to_accession: makes acc lists for each peptide (all proteins that they map to)
        make_unique_share_count_maps: makes peptide sequence to count dictionaries for subsets of peptides (unique or shared)
        total_count_set: computes grand total spectral count of entire passed in set
        combine_proteins: performs the three tests and fixes up the results
        test_redundants: tests for psuedo redundant situations
        test_subsets: test for pseudo subset situations
        test_shares: tests if proteins are similar enough to group together
        sort_accessions: orders proteins for testing from highest unique counts to lowest
        protein_rollup: combines constituent proteins into new family groups
        split_counts: recomputes corrected spectral counts
        make_peptide_rows: makes new peptide_summary output for combined proteins
        get_other_loci: recomputes the OtherLoci lists to reflect family grouping
        trap_combined_accessions: something to do with the one-and-done nature of the subsets test
        sort_combine_list: splits apart compound accessions and sorts all accssions by counts

    written by Phil Wilmarth, OHSU, 2011.
    """

    def __init__(self, pseudo=2.0, low=2.5, share=40.0, mass=20.0, IL=True, diffZ=True):
        """Basic constructor: attributes and definitions listed here:
        """
        self.pseudo = pseudo    # min counts for pseudo-redundants and pseudo-subsets (experiment wide)
        self.low = low          # min ave unique threshold (per sample or per experiment?)
        self.share = share      # ave total count of shared peptides threshold (per sample or per experiment?)
        self.mass = mass        # ratio of shared to unique total counts threshold
        self.IL = IL            # True makes I and L indistinguishable
        self.diffZ = diffZ      # True makes different charge states distinct
        self.number = 0         # passed-in numerical label for console printing (protein group number)
        self.proteins = None    # save pointer to proteins object list (from load_data method)
        self.peptides = None    # save pointer to peptides object list (from load_data method)
        self.samples = []       # list of biological samples in the experiment
        self.pep_to_acc = {}    # map of peptide sequences to accession lists (unique peptides have lists with one acc)
        self.acc_list = []      # list of accessions to be compared (proteins linked by any shared peptides)
        self.sets = []          # list of peptide sets, one for each acc to be compared

        # these are experiment-wide grand totals
        self.tot_count = []     # list of peptide:total count dictionaries for each accession
        self.tot_unique = []    # list of unique peptide:total count dictionaries for each accession
        self.tot_share = []     # list of shared peptide:total count dictionaries for each accession

        # these are counts per sample
        self.tot_samples = []   # nested list peptide:sample count dicts for each acc, sample
        self.uniq_samples = []  # nested list unique peptide:sample count dicts for each acc, sample
        self.share_samples = [] # nested list shared peptide:sample count dicts for each acc, sample

        self.history = ''       # compound accession string capturing combining history
        self.not_to_keep = []   # list of protein accessions that have been combined
        self.log = []           # list of verbose log file lines

        return

    def load_data(self, acc_list, proteins, peptides):
        """Loads data structures for protein combining.
        """
        # make sure there are proteins to compare
        if len(acc_list) <= 1:
            return False

        # set some attributes
        self.acc_list = acc_list
        self.proteins = proteins
        self.peptides = peptides
        self.samples = peptides.samples

        # get peptide set and count dictionary for each protein accession in acc_list
        for acc in acc_list:
            identified_peptides = peptides.tuplets[peptides.acc_map[acc]] # get IDed peptide info for "acc"

            # NOTE: pep_set is created with defaults for IL (True) and diffZ (True)
            """Should test if diffZ does anything. It may well not make any difference..."""
            pep_set, counts = identified_peptides.get_peptide_sets_counts(peptides.col_map, self.samples)
            self.sets.append(pep_set) # save the peptide set for each accession
            self.tot_count.append(counts[-1]) # last dictionary in list is the total counts
            self.tot_samples.append(counts[:-1]) # other dictionaries are the per sample counts

            # set flags to indicate that this protein has been tested
            proteins.tuplets[proteins.acc_map[acc]].done = True

        # build a master dictionary of all peptides and lists of accessions they map to
        self.peptide_to_accession()

        # make unique (to protein collection) and shared peptide count dictionaries
        self.make_unique_share_count_maps()

        return

    def peptide_to_accession(self):
        """Builds a master list of peptide:acc_list for all proteins being compared.
        Peptides with lists of length one are unique with respect to the set of proteins
        being compared. Any peptides with longer lists are shared.
        """
        for i, pep_set in enumerate(self.sets):
            for pep in pep_set:
                if pep in self.pep_to_acc:
                    self.pep_to_acc[pep].append(self.acc_list[i])
                else:
                    self.pep_to_acc[pep] = []
                    self.pep_to_acc[pep].append(self.acc_list[i])
        return

    def make_unique_share_count_maps(self):
        """Makes unique peptides to counts and shared peptides to counts dictionaries.
        Need grand total count (experiment-wide) and per sample count dictionaries.
        """
        for i, acc in enumerate(self.acc_list):

            # create the dictionaries and lists to hold the results, pack into respective lists
            self.tot_unique.append({})
            self.tot_share.append({})
            self.uniq_samples.append([])
            self.share_samples.append([])

            for seq in self.sets[i]:

                # add the peptide keys and counts to the correct dictionaries
                if len(self.pep_to_acc[seq]) == 1:  # unique peptides
                    self.tot_unique[i][seq] = self.tot_count[i][seq] # add total SpC for "unique" peptide
                else:                               # shared peptides
                    self.tot_share[i][seq] = self.tot_count[i][seq] # add total SpC for "shared" peptide

                # do nested per sample stuff next
                for j in range(len(self.samples)):

                    # add the dictionaries to nested lists
                    self.uniq_samples[i].append({})
                    self.share_samples[i].append({})

                    # add the peptide keys and counts to the correct dictionaries
                    if len(self.pep_to_acc[seq]) == 1:  # unique peptides (may not be in every sample)
                        self.uniq_samples[i][j][seq] = self.tot_samples[i][j].get(seq, 0)
                    else:   # shared peptides (may not be present in every sample)
                        self.share_samples[i][j][seq] = self.tot_samples[i][j].get(seq, 0)
        return

    def total_count_set(self, pep_set, count_dict):
        """Returns total SpC for all peptides in pep_set using passed in count dictionary.
        """
        total = 0
        for seq in pep_set:
            total += count_dict.get(seq, 0)
        return total

    def combine_proteins(self, number):
        self.number = number
        """Tests for pseudo-redundants, pseudo-subsets, and finally siblings.
        Each step iterates until the set of proteins is stable.

        There is also considerable book keeping to create new collections
        of information about new family groups for the subsequent reporting.

        The first two tests are similar in concept. In large datasets, there can
        be incorrect matches to correct proteins. It is possible that two fully
        redundant proteins could have small numbers of matches that change
        their identical peptide sets into "nearly" identical sets. We test for
        pairs of proteins where there are mostly shared peptides with small
        numbers of unique peptides. A similar situation can prevent a peptide
        subset from being removed where "nearly" all peptides are a subset of
        larger peptide set. The two situations require different test logic
        but are conceptually more or less the same.

        The final testing is more agressive to find and combine proteins
        with high homologies. Proteins in this category are many
        house keeping genes and well-established familes (actins,
        tubulins, heat shock, HLA, 14-3-3, keratins, etc.). These cases
        become problematic with more complete protein databases but
        are also present even in the most carefully curated databases.
        """
        self.log.append('\n%s. %s' % (number, ', '.join(self.acc_list)))
        redo_split = False

        """The three methods are functionally similar to call. The protein
        group number is passed in (for console printing). The other data
        is kept in object attributes (accession lists and peptide sets).
        What is returned are lists of accesssions to be combined by the
        "protein_rollup" method. The attributes (eg. accession lists) are
        updated during protein rollups. The first and third tests start
        testing a set of accessions or compound accessions (intermediate
        protein families), stop when a grouping has been determined, do the
        grouping, which shortens the accession list attribute. Then the testing
        is started over. Testing stops when all remaining accessions make it
        through the tests and there is nothing to group together.
        """

        # test for pseudo-redundants first
        to_combine = self.test_redundants(number) # need to get the loop started
        while len(to_combine) != 0:
            for acc_list in to_combine: # we only get here if there are acc in to_combine
                self.protein_rollup(acc_list, 'PsRedun')
                redo_split = True
            to_combine = self.test_redundants(number)

        """Subsets is probably the one to modifiy for Scaffold testing.
        """
        # test for pseudo-subsets next
        to_combine = self.test_subsets(number) # this is a one and done test
        for acc_list in to_combine:
            self.trap_combined_accessions(acc_list) # need this check...
            self.protein_rollup(acc_list, 'PsSubset')
            redo_split = True

        # test for similar pairs with significant share overlap (do last)
        to_combine = self.test_shares(number) # need to get the loop started
        while len(to_combine) != 0:
            for acc_list in to_combine: # we only get here if there are acc in to_combine
                self.protein_rollup(acc_list, 'Share')
                redo_split = True
            to_combine = self.test_shares(number)

        self.log.append('.....acc: %s' % (self.acc_list,))
        self.log.append('.....history: %s' % (self.history,))

        # any grouping changes shared and unique definitions and shared peptide splitting
        if redo_split:  # this is only True if something got combined or grouped
            get_counts = self.total_count_set # make short nickname for the function call
            for i, acc in enumerate(self.acc_list):
                if '&' in acc:
                    # need to make new protein entries for the combined groups
                    new = ProteinTuplet()   # create a new data container
                    primary = acc.split('&')[0] # use the first accession as a key
                    new.number = str(len(self.proteins.tuplets)+1) # it will go at the end of the tuplets list
                    new.counter = '1'   # set counter column value
                    new.name = primary + '_family'  # make a proper family accession

                    # new redundant protein list is concatnation of respective lists (if any)
                    redundant_list = []
                    for subacc in acc.split('&'):
                        redundant_list += self.proteins.tuplets[self.proteins.acc_map[subacc]].redundants

                    # several fields get values from the primary protein's data
                    new.filter = self.proteins.tuplets[self.proteins.acc_map[primary]].filter # filter
                    # these values in parentheses to see that they are from primary, extra single quote for Excel
                    new.coverage = float(self.proteins.tuplets[self.proteins.acc_map[primary]].coverage) # coverage
                    new.length = self.proteins.tuplets[self.proteins.acc_map[primary]].length # length
                    new.mw = self.proteins.tuplets[self.proteins.acc_map[primary]].mw # molecular weight
                    new.description = self.proteins.tuplets[self.proteins.acc_map[primary]].description # description

                    # these are recomputed because the grouping changed definitions of shared and unique
                    new.tot_count = str(get_counts(self.sets[i], self.tot_count[i])) # Total counts
                    new.tot_unique = str(get_counts(self.sets[i], self.tot_unique[i])) # total unique counts
                    new.uniqfrac = '%0.3f' % (float(new.tot_unique)/float(new.tot_count)) # unique fraction

                    # recompute the counts per sample and add to their respective lists
                    for j in range(len(self.samples)):
                        new.totals.append(str(get_counts(self.sets[i], self.tot_samples[i][j])))
                        new.uniques.append(str(get_counts(self.sets[i], self.uniq_samples[i][j])))
                        new.correcteds.append(self.split_counts(i, j))
                    new.combined = acc.split('&') # make combined accession list
                    new.otherloci = ', '.join(self.get_other_loci(i)) # make new OtherLoci string
                    self.proteins.tuplets.append(new) # add data object to end of the protein tuplets list

                    # peptides, too
                    newpep = IdentifiedPeptides()   # new data object
                    newpep.group = new.number       # new group number
                    newpep.acc = new.name           # new family accession
                    self.make_peptide_rows(acc, newpep)     # recompute the rows (need to update shared/unique)
                    self.peptides.tuplets.append(newpep)    # add to end of peptide tuplet list

                else:
                    # and recompute corrected counts (per sample) for non-combined proteins
                    new_corr_counts = []
                    tuplet_index = self.proteins.acc_map[acc] # this will be among the original protein tuplets
                    for j in range(len(self.samples)):
                        new_corr_counts.append(self.split_counts(i, j)) # shared splitting done with internal data
                    self.proteins.tuplets[tuplet_index].correcteds = new_corr_counts # replace old corrected counts with new ones
        return

    def test_redundants(self, number):
        """Tests for pseudo-redundant proteins. Breaks after a redundant pair is found.
        This will cause a combination of the pair. The pair's accessions and peptides
        will be merged. This will change shared and unique peptide definitions of peptides
        among the combined pair and the rest of the proteins in the test set. After recomputing
        everthing, we start the testing over again.
        """
        to_combine = []
        self.sort_accessions('descending') # order from most unique to fewest unique

        # upper diagonal loop over decreasing count-sorted proteins
        for i in range(len(self.acc_list)):
            skip = False    # use a flag so we can break out of loop after grouping

            for j in range(i+1, len(self.acc_list)):
                # compare pair-unique to pair-shared total counts ("left" protein is "i" and "right" protein is "j")
                share_set = self.sets[i].intersection(self.sets[j])                     # set of pair-wise shared peptides
                share_count = self.total_count_set(share_set, self.tot_count[j])        # total count of shared peptides (can use either i or j)
                unique_i_set = self.sets[i].difference(self.sets[j])                    # set of "left" unique peptides
                unique_i_count = self.total_count_set(unique_i_set, self.tot_count[i])  # total count of "left" uniques
                unique_j_set = self.sets[j].difference(self.sets[i])                    # set of "right" unique peptides
                unique_j_count = self.total_count_set(unique_j_set, self.tot_count[j])  # total count of "right" peptides
                """left and right unique count totals have to both be less than or equal to
                the PSEUDO parameter AND the total shared count has to be ten times bigger
                than the largest unique total.
                """
                if (unique_i_count <= self.pseudo and
                    share_count >= 10*max(unique_i_count, unique_j_count) and
                    unique_j_count <= self.pseudo):
                    #
                    to_combine.append((self.acc_list[i], self.acc_list[j])) # add the pair to the to_combine list
                    self.log.append('...%s pseudo-redundant with %s [%s (%s) %s]' %
                                    (self.acc_list[i], self.acc_list[j], unique_i_count, share_count, unique_j_count))
                    skip = True
                    break   # breaks out of the inner loop

            if skip:
                break   # breaks out of the outer loop (we need to get out of both after grouping)

        return to_combine

    def test_subsets(self, number):
        """Tests for psuedo-subset proteins. Subsets cascade when proteins are sorted
        by counts so we can do them all in one pass and it is more efficient than
        breaking after each pair to combine.
        """
        to_combine = []
        supersets = {} # this accumulates all test results (values are accession lists)
        self.sort_accessions() # sort from smallest set to largest

        # upper diagonal loop over increasing count-sorted proteins
        for i in range(len(self.acc_list)):
            subsets = {} # this accumulated one "row" (outer loop) of test results
            for j in range(i+1, len(self.acc_list)):

                # compare pair-unique to pair-shared total counts ("left" protein is "i" and "right" protein is "j")
                share_set = self.sets[i].intersection(self.sets[j])                     # set of shared peptides between pair
                share_count = self.total_count_set(share_set, self.tot_count[j])        # total shared counts (could use either i or j)
                unique_i_set = self.sets[i].difference(self.sets[j])                    # "left" unique peptide set
                unique_i_count = self.total_count_set(unique_i_set, self.tot_count[i])  # total counts of "left" uniques
                unique_j_set = self.sets[j].difference(self.sets[i])                    # "right" unique peptide set
                unique_j_count = self.total_count_set(unique_j_set, self.tot_count[j])  # total counts of "right" uniques

                """Pseudo subset if unique total less than or equalt to PSEUDO and shared total
                is ten times larger than unique (or 10 if unique is zero). The first test
                eliminates (maybe) anything that would have been handled with the redundants test
                above.
                """
##                # should we test "j" size or not? Very small unique i or j may fail redundant test due to 10*max(i,j)
##                # but have one pass the subset test here and get removed. I think no test of "j" more true to subset definition
##                if (unique_j_count > self.pseudo and
##                    share_count >= 10*max(unique_i_count, 1) and
##                    unique_i_count <= self.pseudo):

                if (share_count >= 10*max(unique_i_count, 1)) and (unique_i_count <= self.pseudo): # don't care about "right" uniques
                    subsets[self.acc_list[i]] = self.acc_list[j]    # add to the test results dictionary
                                                                    # the dict value will be the largest set in the "row" that contains "i"
                    self.log.append('...%s pseudo-subset of %s [%s (%s) %s]' %
                                    (self.acc_list[i], self.acc_list[j], unique_i_count, share_count, unique_j_count))

            # collect subsets into supersets dictionary (NOTE: this is inside the loop)
            for key, val in subsets.items():
                if val in supersets and key not in supersets.get(val, []):
                    supersets[val].append(key)
                else:
                    supersets[val] = [key]

        # get aggregate accession lists for combining
        for key in supersets:
            supersets[key].append(key)  # this adds "self" to the lists
            to_combine.append(supersets[key])   # NOTE: to_combine will be list of lists

        return to_combine

    def test_shares(self, number):
        """Tests pairs of proteins in acc_list to see if they share enough
        peptides to be combined. Breaks after a pair to be grouped is found.
        This will cause a combination of the pair. The pair's accessions and peptides
        will be merged. This will change shared and unique peptide definitions of peptides
        among the combined pair and the rest of the proteins in the test set. After recomputing
        everthing, we start the testing over again.        """
        to_combine = []
        self.sort_accessions('descending')  # sort by decreasing unique counts

        # upper diagonal loop over decreasing count-sorted proteins
        for i in range(len(self.acc_list)):
            skip = False    # use a flag so we can break out of loop after grouping

            for j in range(i+1, len(self.acc_list)):
                # compare pair-unique to pair-shared total counts ("left" protein is "i" and "right" protein is "j")
                share_set = self.sets[i].intersection(self.sets[j])                     # set of pair-wise shared peptides
                share_count = self.total_count_set(share_set, self.tot_count[j])        # total count of shared peptides (can use either i or j)
                ave_share_count = share_count / float(len(self.samples))                # ave shared count per sample
                unique_i_set = self.sets[i].difference(self.sets[j])                    # set of "left" unique peptides
                unique_i_count = self.total_count_set(unique_i_set, self.tot_count[i])  # total "left" unique counts
                unique_j_set = self.sets[j].difference(self.sets[i])                    # set of "right" unique peptides
                unique_j_count = self.total_count_set(unique_j_set, self.tot_count[j])  # total "right" unique counts
                ave_i_unique = unique_i_count / float(len(self.samples))                # ave unique per sample for "left"
                ave_j_unique = unique_j_count / float(len(self.samples))                # ave unique per sample for "right"

                if unique_i_count > 0:
                    share_i_mass = share_count / float(unique_i_count)  # relative weight (mass) of shared to uniques for "left"
                else:
                    share_i_mass = 1000.0                               # something big instead of divide by zero
                if unique_j_count > 0:
                    share_j_mass = share_count / float(unique_j_count)  # relative weight (mass) of shared to uniques for "right"
                else:
                    share_j_mass = 1000.0                               # something big instead of divide by zero

                """We compare "left" unique evidence to the shared evidence. If left unique is too small
                OR the shared evidence overwhelms the left unique evidence, we group. We do the same thing
                for the "right" evidence. Either "left" or "right" evidence can trigger grouping.
                """
                # see if we need to combine the pair
                if (((ave_i_unique < self.low) and (ave_share_count > self.share)) or
                    (share_i_mass > self.mass) or
                    ((ave_j_unique < self.low) and (ave_share_count > self.share)) or
                    (share_j_mass > self.mass)):

                    to_combine.append((self.acc_list[i], self.acc_list[j]))     # add the pair to the to_combine list
                    self.log.append('...%s to be combined with %s [%s (%s) %s]' %
                                    (self.acc_list[i], self.acc_list[j], unique_i_count, share_count, unique_j_count))
                    skip = True
                    break   # break out of the inner loop
                else:
                    # also useful to see what is NOT combined
                    line = ('...%s NOT combined with %s [%s (%s) %s]' %
                            (self.acc_list[i], self.acc_list[j], unique_i_count, share_count, unique_j_count))
                    if line not in self.log: # this keep console output under control
                        self.log.append(line)

            if skip:
                break # break out of the outer loop (need to get out of both loops after a grouping)

        return to_combine

    def sort_accessions(self, direction='ascending'):
        """Sorts accession indexes by spectral counts (total then unique)
        then sorts lists in self in the same order.

        Sorts can be either ascending (default) or descending. Sort tuples
        are created with leading counts to achieve correct sorting logic.
        An array of indexes in the sorted order is computed so that several
        other lists can be rearranged to match the sorted order.

        """
        # set flag for sort order
        REV = False
        if direction == 'descending':
            REV = True

        # get new index order
        indexes = []
        for i in range(len(self.acc_list)):
            total = self.total_count_set(self.sets[i], self.tot_count[i]) # get total count
            unique = self.total_count_set(self.sets[i], self.tot_unique[i]) # get unique count
            indexes.append((unique, total, i)) # make a sort tuple
        indexes.sort(reverse=REV) # sort tuples in the desired order
        indexes = [x[2] for x in indexes] # extract the list indexes from the sorted tuples

        # put lists in new order, reassign the results to the originals
        self.acc_list = [self.acc_list[i] for i in indexes]
        self.sets = [self.sets[i] for i in indexes]
        self.tot_count = [self.tot_count[i] for i in indexes]
        self.tot_unique = [self.tot_unique[i] for i in indexes]
        self.tot_share = [self.tot_share[i] for i in indexes]
        self.tot_samples = [self.tot_samples[i] for i in indexes]
        self.uniq_samples = [self.uniq_samples[i] for i in indexes]
        self.share_samples = [self.share_samples[i] for i in indexes]
        return

    def protein_rollup(self, combine_list, why='Share'):
        import copy
        """Combines a list of accessions into a single entry and
        updates data structure.

        "combine_list" is an accession list (they can be simple or
        compound accessions). "why" is a text string to indicate why
        the things were combined (PsRedundant, PsSubset, or Share).

        written by Phil Wilmarth, OHSU, 2011.

        Original protein and peptide data structures are flagged so that
        they do not get written to output files. Combined data strutures
        are created, populated, and added to end of the data structure lists.

        Note that two of the tests are iterative and intermediate pairs may be
        grouped along the way towards even larger groups. Those intermediate groups
        may not be present in the final output if they get subsumed by larger
        groups.
        """
        # make new accession label with accessions sorted by decreasing counts
        new_acc_list = self.sort_combine_list(combine_list)
        new_acc = '&'.join(new_acc_list)

        # make new history strings (for tracebacks, and verbose output)
        if self.history != '':
            self.history = '&'.join([why] + new_acc_list + ['(' + self.history + ')'])
        else:
            self.history = '&'.join([why] + new_acc_list)

        # flag proteins that are being combined (set protein and peptide "keep" flags to False)
        for acc in combine_list:
            if acc == '':
                continue
            if '&' not in acc:
                self.not_to_keep.append(acc)
                self.proteins.tuplets[self.proteins.acc_map[acc]].keep = False # flag to skip for output
                self.peptides.tuplets[self.peptides.acc_map[acc]].keep = False # flag to skip for output

        # make acc:index dictionary for combine_list
        index = {}
        for acc in combine_list:
            index[acc] = self.acc_list.index(acc)

        # get list of index values and sort descending
        index_list = list(index.values())
        index_list.sort(reverse=True)

        # make a union of all peptides associated with combine_list
        # also update total count dictionaries
        pep = set()
        count = {}
        tot_samples = [{} for s in self.samples]
        for acc in combine_list:
            pep = pep.union(self.sets[index[acc]])
            count.update(self.tot_count[index[acc]])
            for j in range(len(self.samples)):
                tot_samples[j].update(self.tot_samples[index[acc]][j])

        # remove (delete) data associated with original proteins
        for i in index_list:
            del self.acc_list[i]
            del self.sets[i]
            del self.tot_count[i]
            del self.tot_samples[i]

        # add new accession, peptide set, and count dictionaries to self
        self.acc_list.append(new_acc)   # NOTE "self.acc_list[-1]" is "new_acc"
        self.sets.append(pep)
        self.tot_count.append(count)
        self.tot_samples.append(tot_samples)

        # for each peptide, replace individual accessions with combined accession
        for pep in self.pep_to_acc:
            affected = False
            for acc in combine_list:
                try:    # if this runs, we have altered accessions
                    self.pep_to_acc[pep].remove(acc)
                    affected = True
                except ValueError: # when above fails, nothing to change
                    pass
            if affected: # add combined accession to altered lists
                self.pep_to_acc[pep].append(new_acc)

        # recompute unique (to protein collection) and shared peptide count dictionaries
        self.tot_unique = []
        self.tot_share = []
        self.uniq_samples = []
        self.share_samples = []
        self.make_unique_share_count_maps()
        #
        return

    def split_counts(self, index, sample):
        """Calculates corrected totals after splitting shared peptide counts.
        This may improve spectral counting total counts, if anyone cares anymore...
        """
        # make accession to acc_list index map and total unique counts for each accession map
        acc_map = {}
        acc_to_uniqtot = {}
        for i, acc in enumerate(self.acc_list):
            acc_map[acc] = i
            """Can drive splitting two ways, by unqiues per sample or by uniques per experiment.
            per sample might have unreliable information if counts are low, per experiment may not
            be the right thing to do if samples are very different from each other
            """
##            acc_to_uniqtot[acc] = self.total_count_set(self.sets[i], self.uniq_samples[i][sample]) # use uniques per sample
            acc_to_uniqtot[acc] = self.total_count_set(self.sets[i], self.tot_unique[i])    # use uniques per experiment

        # loop over all peptides mapped to accession (pointed to by "index")
        prot_total = 0.0
        for pep in self.sets[index]:

            # get total counts for "pep" in sample "sample"
            pep_total = self.total_count_set(set([pep]), self.tot_samples[index][sample])
            if pep_total == 0:
                continue

            # get relative unique totals for all proteins that pep maps to (itself and others)
            other_unique_total = 0
            self_unique_total = 0 # this may not be necessary
            for acc in self.pep_to_acc[pep]:
                i = acc_map[acc]
                if i == index:
                    self_unique_total = acc_to_uniqtot[acc]
                else:
                    other_unique_total += acc_to_uniqtot[acc]
            #
            unique_total = self_unique_total + other_unique_total

            # split peptide's counts based on relative unique counts
            # (or split equally if no unique counts)
            if unique_total > 0:
                prot_total += (pep_total * self_unique_total) / unique_total
            else:
                prot_total += pep_total / float(len(self.pep_to_acc[pep]))
        #
        return ('%0.3f' % (prot_total,)) # note: this returns a string not a number

    def make_peptide_rows(self, acc, newpep):
        """Makes new rows for peptide summary for combined proteins.

        There is not much processing here. It may not actually provide
        a properly merged set of peptides. There might be issues with
        peptides with I and L amino acids. This needs some more testing...
        """
        import re

        seq_seen = {}
        subacc_list = acc.split('&')
        for subacc in subacc_list:
            tuplet_index = self.peptides.acc_map[subacc]
            for row in self.peptides.tuplets[tuplet_index].rows:
                col = row.split('\t')
                col[0] = newpep.group
                col[1] = newpep.acc

                # skip rows if peptide has already been seen in previous protein
                seq = col[2].split('.')[1] + '.' + col[7]   # sequence plus charge

##                if self.IL:
##                    seq = re.sub(r'[IL]', 'j', seq)     # mask I and L with j (make peptide sequence mass spec-like)

                if seq in seq_seen:
                    continue
                else:
                    seq_seen[seq] = True

                # build new OtherLoci string
                new_other = []
                if len(col) == len(self.peptides.col_map): # rows without OtherLoci will be shorter and be skipped
                    new_other = [acc]   # start with self (a compound accession) in accession list
                    for old_other in [x.strip() for x in col[-1].split(',')]:
                        if old_other and old_other not in subacc_list:
                            new_other.append(old_other)     # add any accessions not already in the compound accession

                    # we do not make OtherLoci with just self (use blank string instead)
                    if len(new_other) > 1:
                        col[-1] = ', '.join(new_other)
                    else:
                        col[-1] = ' '
                        new_other = []

                # update unique status (if no OtherLoci, peptide is unique)
                if len(new_other) == 0:
                    col[5] = 'TRUE'
                else:
                    col[5] = 'FALSE'

                # add the new row
                newpep.rows.append('\t'.join(col))
        return

    def sort_combine_list(self, to_combine):
        """Orders individual accessions in combined accessions from high to low unique counts.
        """
        temp = '&'.join(to_combine)
        temp = temp.split('&')  # these two statements effectively "uncombines" all accessions

        to_sort = []
        for acc in temp:
            if acc == '':
                print('...WARNING: empty accession string')
                continue
            
            # make sort tuples so sorting is by counts not aphabetical by accessions
            tuplet = self.proteins.tuplets[self.proteins.acc_map[acc]]
            try:
                to_sort.append((int(tuplet.tot_unique), int(tuplet.tot_count), acc))
            except ValueError:
                to_sort.append((float(tuplet.tot_unique), float(tuplet.tot_count), acc))
        to_sort.sort(reverse=True)
        #
        return [x[2] for x in to_sort] # extract the accessions from the sorted list of tuples

    def get_other_loci(self, index):
        """Gets a list of all other loci sharing any peptides with index's protein.
        """
        other = set()
        for pep in self.sets[index]:
            for acc in self.pep_to_acc[pep]:
                other.add(acc)  # this makes a non-redundant set of accessions
        other = list(other)
        other.sort()

        # we do not include self in OtherLoci lists, use empty list instead
        if len(other) == 1:
            other = []
        return other

    def trap_combined_accessions(self, combine_list):
        import copy
        """Removes any accessions that have already been combined.

        "combine_list" is accessed by reference so when we
        modify combine_list, it is changed in the calling location.
        We are not working on a copy of the list.
        """
        test = {}
        to_delete = []
        for i, acc in enumerate(combine_list):
            try:
                test[acc] = self.acc_list.index(acc)
            except ValueError:
                to_delete.append(i)
        #
        to_delete.sort(reverse=True)
        for i in to_delete:
            del combine_list[i]
        return

    # end class OtherLociCombiner

def main(prot_file=None):
    """
    Processes PAW protein and peptide reports and creates a quantitation
    summary file where similar proteins are grouped.

    written by Phil Wilmarth, OHSU, 2011.
    """
    #=======================================================================
    # Takes identified peptide and protein results files and processes
    # data to make subsequent quantitative analysis easier.
    # Rolls up redundant proteins into first group member row.
    # Creates combined contaminant and decoy rows (originals removed).
    # Groups highly similar proteins and recalculates shared peptide splits.
    # Removes some extraneous columns to simplify the resulting table.
    # written by Phil Wilmarth, OHSU, 2011.
    #=======================================================================
    #
    # print name, version information
    print('===============================================================')
    print(' PAW_protein_grouper.py, %s, written by Phil Wilmarth, OHSU' % VERSION)
    print('===============================================================')

    """Use a dialog box to get the protein summary file. Try to get
    the peptide summary file name based on protein file name. Open
    a log file in the same location as the protein file.
    """
    # browse to protein and peptide results file(s)
    if not prot_file:
        print('Select PROTEIN results file to process')
        default_location = 'F:\PSR_Core_Analysis'
        if not os.path.exists(default_location):
            default_location = os.getcwd()
        ext_list = [('Text files', '*.txt')]
        title = 'Select PROTEIN results file to process'
        prot_file = get_file(default_location, ext_list, title)
        if prot_file == '': sys.exit()

    pep_file = prot_file.replace('protein', 'peptide')
    if not os.path.exists(pep_file):
        default = os.path.dirname(prot_file)
        title = 'Select corresponding PEPTIDE results file'
        pep_file = get_file(default, ext_list, title)
        if pep_file == '': sys.exit()
        if os.path.dirname(prot_file) != os.path.dirname(pep_file):
            print('...WARNING: protein and peptide files should be in the same directory.')

    # open the log file
    prot_path = os.path.dirname(prot_file)
    log_obj = open(os.path.join(prot_path, 'PAW_protein_grouper.log'), 'wt')
    time_stamp_logfile('\n>>> starting: PAW_protein_grouper.py %s' % VERSION, log_obj)
    write = [None, log_obj]

    # create a ProteinSummary object and use its method to read in the protein summary file
    proteins = ProteinSummary()
    proteins.load_results(prot_file)

    # create a PeptideSummary object and use its method to read in the peptide summary file
    peptides = PeptideSummary(IL, diffZ)
    peptides.load_results(pep_file)

    print('files have been read in')

    """Both ProteinSummary and PeptideSummary objects are kind of meta data containers.
    Each of them have list attributes that hold lists of ProteinTuplet objects and
    lists of IdentifiedPeptide objects. The ProteinTuplets and IdentifiedPeptide objects
    are the main data containers. There is one instance of each object for each specific
    protein or specific protein group or specific protein family. Furthermore, they are
    linked by protein accession keys and the two list of objects are kept in sync with
    each other. IdentifiedPeptide objects could have been added as attributes to
    ProteinTuplets since they go together. They were kept separate to facilitate reading in
    two different summary files and creating two new summary files.
    """
    # initialize counters, etc.
    candidates = 0  # counter of all proteins that were compared
    keepers = 0     # candidate count corrected because of combining
    keeper_list = []

    """Loop logic has many options. Proteins that only have unique peptides
    do not need to be tested. The comparison testing is done on small sets
    of proteins that have shared peptides in common. Rather than trying to
    collect all of those small sets together in some structure that could
    be looped over, we will do a linear traversal over all of the proteins
    that were read in. We will use the list of protein tuplets not only to
    keep track of the proteins we started with but we will also add any
    new groups/families to the end of the list. We will use some flag attributes
    in the tuplet objects to keep track of things. When we encounter a tuplet
    that has a list of OtherLoci, we will pass those to the combiner object.
    All of those proteins will get flagged as having been tested so we don't
    test groups multiple times. If proteins are combined into any families,
    a new tuplet is added to the end of the protein list for the new family.
    The members in that familt get flagged to be removed. We don't remove them
    as much as skip writing them to the new output file.
    """
    # loop over protein tuplets created from reading the summary files
    for i in range(proteins.protein_count):
        prot = proteins.tuplets[i]
        if prot.otherloci == '':    # unique proteins are already "done"
            prot.done = True
        elif not prot.done:     # skip proteins that have already been through testing (attribute set during testing)
            acc_list = list(prot.all_loci) # proteins to compare to each other
            acc_list.sort()     # probably don't need to do this except to make console output cleaner
            combiner = OtherLociCombiner(PSEUDO, LOW, SHARE, MASS, IL, diffZ) # create a combiner object
            combiner.load_data(acc_list, proteins, peptides) # load data for the proteins being compared
            candidates += len(acc_list)     # keep count of how many proteins were tested
            combiner.combine_proteins(i+1) # test proteins for similarity and group if needed (pass in index for console and log file)
            keepers += len(combiner.acc_list) # this counts the new, net number of proteins/groups/families
            keeper_list += combiner.acc_list    # list of all final protein/group/family accessions (some simple, some compound)

            # write console and log lines (combiner saves output lines so file lists do not have to be passed in)
            if VERBOSE:     # write more detailed log file if desired
                for obj in write:
                    for line in combiner.log:
                        print(line, file=obj)
            else:
                for obj in write[1:]:
                    print('\nFor protein family %s:' % (i+1,), file=obj)
                    print('In: ', acc_list, file=obj)
                    print('Out:', combiner.acc_list, file=obj)

    # print out new protein and peptide summaries
    new_prot_file = os.path.join(os.path.dirname(prot_file), 'grouped_'+os.path.basename(prot_file))
    new_pep_file = os.path.join(os.path.dirname(pep_file), 'grouped_'+os.path.basename(pep_file))
    new_prot_file_obj = open(new_prot_file, 'w')
    new_pep_file_obj = open(new_pep_file, 'w')
    proteins.print_summary(new_prot_file_obj, write)
    peptides.print_summary(new_pep_file_obj, write)

    # print summary of families
    for obj in write:
        print('\nFamily Member Summaries:', file=obj)
    for new_family in proteins.tuplets[proteins.protein_count:]:
        for obj in write:
            print('\n..Family:', new_family.name, file=obj)
            for acc in new_family.combined:
                member = proteins.tuplets[proteins.acc_map[acc]]
                print('....%s: %s' % (member.name, member.description), file=obj)
                if member.name.endswith(')'):
                    for redundant in member.redundants[1:]:
                        print('......%s: %s' % (redundant, proteins.acc_to_desc[redundant]), file=obj)


    # print some summary stats and parameter settings
    family_count = len(proteins.tuplets[proteins.protein_count:])
    distinguishable_count = (proteins.protein_count - candidates) + (keepers - family_count)
    for obj in write:
        print(file=obj)
        print('PAW_protein_grouper settings:', file=obj)
        print('PSEUDO:\t', PSEUDO, '\t(used in pseudo-redudant and pseudo-subset tests)', file=obj)
        print('LOW:\t', LOW, '\t(used in combination test - min. ave. true unique count)', file=obj)
        print('SHARE:\t', SHARE, '\t(used in combination test - absolute ave. share count)', file=obj)
        print('MASS:\t', MASS, '\t(used in combination test - relative shared-to-unique fraction)', file=obj)       
        print(file=obj)
        print('Number of proteins read in:', proteins.row_count, file=obj)
        print('...%d distinguishable peptide sets (proteins) and %d fully redundant peptide sets (groups)' %
              (proteins.distinguishable_count, proteins.indistinguishable_count), file=obj)
        print('Total number of proteins/groups is:', proteins.protein_count, file=obj)
        print('Total number of proteins/groups containing shared peptides:', candidates, file=obj)
        print('Candidates tested for highly homologous protein families:', candidates, file=obj)
        print('...%d proteins/groups and %d combined family groups' % (keepers-family_count, family_count), file=obj)
        print('Final total number of proteins/groups/families:', distinguishable_count+family_count, file=obj)

    # close files
    new_prot_file_obj.close()
    new_pep_file_obj.close()
    log_obj.close()
    # end

# created main function so grouping can be run from the results script
if __name__ == '__main__':
    try:
        if os.path.exists(sys.argv[1]):
            main(sys.argv[1])
    except IndexError:
        main()
    else:
        main()
