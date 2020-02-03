"""comet_GUI.py
Simple GUI interface for setting common Comet search parameters.

written by Delan Huang, OHSU, 2014
additions by Phil Wilmarth, OHSU, 2014.

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

# need to change license to MIT
# converted for Python 3 -PW 9/15/2017


###########################
# NOTE: there is no input validation checking for numerical fields
###########################

# default Comet params file contents (this is for version 2016.01 rev. 3)
comet_default_params = """# comet_version 2016.01 rev. 3
# Comet MS/MS search engine parameters file.
# Everything following the '#' symbol is treated as a comment.

database_name = /some/path/db.fasta
decoy_search = 0                       # 0=no (default), 1=concatenated search, 2=separate search

num_threads = 0                        # 0=poll CPU to set num threads; else specify num threads directly (max 64)

#
# masses
#
peptide_mass_tolerance = 3.00
peptide_mass_units = 0                 # 0=amu, 1=mmu, 2=ppm
mass_type_parent = 1                   # 0=average masses, 1=monoisotopic masses
mass_type_fragment = 1                 # 0=average masses, 1=monoisotopic masses
precursor_tolerance_type = 0           # 0=MH+ (default), 1=precursor m/z; only valid for amu/mmu tolerances
isotope_error = 0                      # 0=off, 1=on -1/0/1/2/3 (standard C13 error), 2= -8/-4/0/4/8 (for +4/+8 labeling)

#
# search enzyme
#
search_enzyme_number = 1               # choose from list at end of this params file
num_enzyme_termini = 2                 # 1 (semi-digested), 2 (fully digested, default), 8 C-term unspecific , 9 N-term unspecific
allowed_missed_cleavage = 2            # maximum value is 5; for enzyme search

#
# Up to 9 variable modifications are supported
# format:  <mass> <residues> <0=variable/else binary> <max_mods_per_peptide> <term_distance> <n/c-term> <required>
#     e.g. 79.966331 STY 0 3 -1 0 0
#
variable_mod01 = 15.9949 M 0 3 -1 0 0
variable_mod02 = 0.0 X 0 3 -1 0 0
variable_mod03 = 0.0 X 0 3 -1 0 0
variable_mod04 = 0.0 X 0 3 -1 0 0
variable_mod05 = 0.0 X 0 3 -1 0 0
variable_mod06 = 0.0 X 0 3 -1 0 0
variable_mod07 = 0.0 X 0 3 -1 0 0
variable_mod08 = 0.0 X 0 3 -1 0 0
variable_mod09 = 0.0 X 0 3 -1 0 0
max_variable_mods_in_peptide = 5
require_variable_mod = 0

#
# fragment ions
#
# ion trap ms/ms:  1.0005 tolerance, 0.4 offset (mono masses), theoretical_fragment_ions = 1
# high res ms/ms:    0.02 tolerance, 0.0 offset (mono masses), theoretical_fragment_ions = 0
#
fragment_bin_tol = 1.0005              # binning to use on fragment ions
fragment_bin_offset = 0.4              # offset position to start the binning (0.0 to 1.0)
theoretical_fragment_ions = 1          # 0=use flanking peaks, 1=M peak only
use_A_ions = 0
use_B_ions = 1
use_C_ions = 0
use_X_ions = 0
use_Y_ions = 1
use_Z_ions = 0
use_NL_ions = 0                        # 0=no, 1=yes to consider NH3/H2O neutral loss peaks

#
# output
#
output_sqtstream = 0                   # 0=no, 1=yes  write sqt to standard output
output_sqtfile = 0                     # 0=no, 1=yes  write sqt file
output_txtfile = 0                     # 0=no, 1=yes  write tab-delimited txt file
output_pepxmlfile = 1                  # 0=no, 1=yes  write pep.xml file
output_percolatorfile = 0              # 0=no, 1=yes  write Percolator tab-delimited input file
output_outfiles = 0                    # 0=no, 1=yes  write .out files
print_expect_score = 1                 # 0=no, 1=yes to replace Sp with expect in out & sqt
num_output_lines = 5                   # num peptide results to show
show_fragment_ions = 0                 # 0=no, 1=yes for out files only

sample_enzyme_number = 1               # Sample enzyme which is possibly different than the one applied to the search.
                                       # Used to calculate NTT & NMC in pepXML output (default=1 for trypsin).

#
# mzXML parameters
#
scan_range = 0 0                       # start and scan scan range to search; 0 as 1st entry ignores parameter
precursor_charge = 0 0                 # precursor charge range to analyze; does not override any existing charge; 0 as 1st entry ignores parameter
override_charge = 0                    # 0=no, 1=override precursor charge states, 2=ignore precursor charges outside precursor_charge range, 3=see online
ms_level = 2                           # MS level to analyze, valid are levels 2 (default) or 3
activation_method = ALL                # activation method; used if activation method set; allowed ALL, CID, ECD, ETD, PQD, HCD, IRMPD

#
# misc parameters
#
digest_mass_range = 600.0 5000.0       # MH+ peptide mass range to analyze
num_results = 100                      # number of search hits to store internally
skip_researching = 1                   # for '.out' file output only, 0=search everything again (default), 1=don't search if .out exists
max_fragment_charge = 3                # set maximum fragment charge state to analyze (allowed max 5)
max_precursor_charge = 6               # set maximum precursor charge state to analyze (allowed max 9)
nucleotide_reading_frame = 0           # 0=proteinDB, 1-6, 7=forward three, 8=reverse three, 9=all six
clip_nterm_methionine = 0              # 0=leave sequences as-is; 1=also consider sequence w/o N-term methionine
spectrum_batch_size = 0                # max. # of spectra to search at a time; 0 to search the entire scan range in one loop
decoy_prefix = DECOY_                  # decoy entries are denoted by this string which is pre-pended to each protein accession
output_suffix =                        # add a suffix to output base names i.e. suffix "-C" generates base-C.pep.xml from base.mzXML input
mass_offsets =                         # one or more mass offsets to search (values substracted from deconvoluted precursor mass)

#
# spectral processing
#
minimum_peaks = 10                     # required minimum number of peaks in spectrum to search (default 10)
minimum_intensity = 0                  # minimum intensity value to read in
remove_precursor_peak = 0              # 0=no, 1=yes, 2=all charge reduced precursor peaks (for ETD)
remove_precursor_tolerance = 1.5       # +- Da tolerance for precursor removal
clear_mz_range = 0.0 0.0               # for iTRAQ/TMT type data; will clear out all peaks in the specified m/z range

#
# additional modifications
#

add_Cterm_peptide = 0.0
add_Nterm_peptide = 0.0
add_Cterm_protein = 0.0
add_Nterm_protein = 0.0

add_G_glycine = 0.0000                 # added to G - avg.  57.0513, mono.  57.02146
add_A_alanine = 0.0000                 # added to A - avg.  71.0779, mono.  71.03711
add_S_serine = 0.0000                  # added to S - avg.  87.0773, mono.  87.03203
add_P_proline = 0.0000                 # added to P - avg.  97.1152, mono.  97.05276
add_V_valine = 0.0000                  # added to V - avg.  99.1311, mono.  99.06841
add_T_threonine = 0.0000               # added to T - avg. 101.1038, mono. 101.04768
add_C_cysteine = 57.021464             # added to C - avg. 103.1429, mono. 103.00918
add_L_leucine = 0.0000                 # added to L - avg. 113.1576, mono. 113.08406
add_I_isoleucine = 0.0000              # added to I - avg. 113.1576, mono. 113.08406
add_N_asparagine = 0.0000              # added to N - avg. 114.1026, mono. 114.04293
add_D_aspartic_acid = 0.0000           # added to D - avg. 115.0874, mono. 115.02694
add_Q_glutamine = 0.0000               # added to Q - avg. 128.1292, mono. 128.05858
add_K_lysine = 0.0000                  # added to K - avg. 128.1723, mono. 128.09496
add_E_glutamic_acid = 0.0000           # added to E - avg. 129.1140, mono. 129.04259
add_M_methionine = 0.0000              # added to M - avg. 131.1961, mono. 131.04048
add_O_ornithine = 0.0000               # added to O - avg. 132.1610, mono  132.08988
add_H_histidine = 0.0000               # added to H - avg. 137.1393, mono. 137.05891
add_F_phenylalanine = 0.0000           # added to F - avg. 147.1739, mono. 147.06841
add_U_selenocysteine = 0.0000          # added to U - avg. 150.3079, mono. 150.95363
add_R_arginine = 0.0000                # added to R - avg. 156.1857, mono. 156.10111
add_Y_tyrosine = 0.0000                # added to Y - avg. 163.0633, mono. 163.06333
add_W_tryptophan = 0.0000              # added to W - avg. 186.0793, mono. 186.07931
add_B_user_amino_acid = 0.0000         # added to B - avg.   0.0000, mono.   0.00000
add_J_user_amino_acid = 0.0000         # added to J - avg.   0.0000, mono.   0.00000
add_X_user_amino_acid = 0.0000         # added to X - avg.   0.0000, mono.   0.00000
add_Z_user_amino_acid = 0.0000         # added to Z - avg.   0.0000, mono.   0.00000

#
# COMET_ENZYME_INFO _must_ be at the end of this parameters file
#
[COMET_ENZYME_INFO]
0.  No_enzyme              0      -           -
1.  Trypsin                1      KR          P
2.  Trypsin/P              1      KR          -
3.  Lys_C                  1      K           P
4.  Lys_N                  0      K           -
5.  Arg_C                  1      R           P
6.  Asp_N                  0      D           -
7.  CNBr                   1      M           -
8.  Glu_C                  1      DE          P
9.  PepsinA                1      FL          P
10. Chymotrypsin           1      FWYL        P

"""

# global imports
from tkinter import *
import tkinter.ttk as ttk
import PAW_lib
import os
import sys

class staticMods(Toplevel):
    """Creates a static modification top-level widget.
    """
    def __init__(self, parent, default_params, static_list):
        """Constructor
        """
        Toplevel.__init__(self, parent)
        self.transient(parent)
        self.default_params = default_params
        self.static_list = static_list
        self.parent = parent
        self.title('Static modifications')
        self.attributes('-topmost', 1)

        # parse the static modifications from the default params contents
        self.keys = [] # this keeps the order of the static mods
        self.static_dict = {}
        self.parse_static_mods()

        # create a static mods frame with labels and entries
        self.create_static_mods_frame()
        self.initial_focus = self.static_mods_frame
        self.buttonbox()

        # this makes the static mods a modal widget
        self.grab_set()
        self.protocol('WM_DELETE_WINDOW', self.onDone)
        self.geometry('+%d+%d' % (parent.winfo_rootx()+50, parent.winfo_rooty()+50))
        self.initial_focus.focus_set()
        self.wait_window(self)
        return

    def buttonbox(self):
        """Create some action buttons below main GUI elements.
        """
        box = Frame(self)
        w1 = Button(box, text='TMT labels', width=10, command=self.TMT)
        w1.pack(side=LEFT, padx=5, pady=5)
        box.pack()
        w = Button(box, text='Done', width=10, command=self.onDone)
        w.pack(side=LEFT, padx=5, pady=5)
        self.bind('<Return>', self.onDone)
        box.pack()

    def parse_static_mods(self):
        """Parses the static modification information from the default params file.
        """
        for line in self.default_params.splitlines():
            line = line.strip()
            if line.startswith('add_'):
                key = line.split('=')[0].strip()
                self.keys.append(key)
                value = line.split('=')[1].split('#')[0].strip()
                try:
                    comment = line.split('=')[1].split('#')[1].strip()
                except IndexError:
                    comment = ''
                self.static_dict[key] = (value, comment)
            else:
                continue
        return
    
    def create_static_mods_frame(self):
        """Creates a grid layout of the static modifications with entries for deltamass.
        """
        self.static_mods_frame = ttk.Labelframe(self, text='Static modifications:')
        self.static_mods_frame.pack(fill=X, expand=YES, padx=5, pady=5)
        
        #Variables
        self.static_mass = {}
        
        #Creation
        headers = ('Residue/Position', 'Static deltamass', 'Comments')
        locations = ('e', 'we', 'w')
        for i, header in enumerate(headers):
            Label(self.static_mods_frame, text=header).grid(column=i,row=0, sticky=locations[i])
            
        for i, key in enumerate(self.keys):
            # get values from parsed dictionary and set entry fields
            mass, comment = self.static_dict[key]
            self.static_mass[i] = DoubleVar()
            self.static_mass[i].set(float(mass))

            # grid the widgets
            Label(self.static_mods_frame, text=key).grid(column=0, row=i+1, sticky=E)
            Entry(self.static_mods_frame, textvariable=self.static_mass[i]).grid(column=1, row=i+1)
            Label(self.static_mods_frame, text=comment).grid(column=2, row=i+1, sticky=W)

        return

    def TMT(self, event=None):
        """Sets static mods for 10-plex TMT labeling.
        """
        self.static_dict['add_K_lysine'] = ('229.1629320', 'added to K - avg. 128.1723, mono. 128.09496')
        self.static_dict['add_Nterm_peptide'] = ('229.1629320', '')
        for i, key in enumerate(self.keys):
            # get values from parsed dictionary and set entry fields
            mass, comment = self.static_dict[key]
            self.static_mass[i].set(float(mass))
##        print('\nstatic dict:')
##        print(self.static_dict)
        
    def onDone(self, event=None):
        """Reads out the deltamass values and loads them into a passed in list pointer, then exits.
        """
        for i, key in enumerate(self.keys):
            self.static_list.append((key, self.static_mass[i].get(), self.static_dict[key][1]))

        self.withdraw()
        self.update_idletasks()
        self.parent.focus_set()
        self.destroy()

    # end class

class CometGUI:
    """Main GUI for setting a few of the common Comet parameters. A modified
    comet.params file is written and comet searches can be launched.
    """
    def __init__(self):
        """Constructor. Uses a collection of frames for related parameters.
        """
        self.root = Tk()
        self.root.title('Comet Parameters')
        self.root.protocol('WM_DELETE_WINDOW', self.quit_gui)
#        self.root.attributes('-topmost', 1)
#        self.root.attributes('-topmost', 0)

        # set some default attributes
        self.ms2_folder = None
        self.filename = None
        
        #create GUI
        self.create_dir_frame()
        self.create_masses_frame()
        self.create_use_ion_frame()
        self.create_search_enzyme_frame()
        self.create_variable_mods_frame()

        # buttons at bottom of GUI
        Button(self.root, text='Change static modifications', command=self.change_static).pack(pady=2)
        self.save_params = Button(self.root, text='Save Settings and Create Parameters File',
                                  command=self.save_settings, state=NORMAL)
        self.save_params.pack(pady=2)
        self.run_comet = Button(self.root, text='Run Comet', command=self.run_comet)
        self.run_comet.pack(pady=2)
        Button(self.root, text='Quit', command=self.quit_gui).pack(pady=2)


        # place holder for static mod information
        self.static_list = []

        # enter mainloop
        self.root.mainloop()
        return
        
    #Functions to help create widgets             
    def create_entry(self, root, label, variable):
        """Creates a text entry widget.
        """
        frame = Frame(root)
        Label(frame, text=label).pack(side=LEFT)
        Entry(frame, textvariable=variable).pack(side=LEFT)
        return frame

    def create_checkboxes(self, root, label, ions_list):
        """Creates a checkbox widget.
        """
        frame = Frame(root)
        Label(frame, text=label).pack(side=LEFT)
        for button, variable in ions_list:
            Checkbutton(frame, text=button, variable=variable).pack(side=LEFT)
        return frame
                
    def create_radiobuttons(self, root, label, buttons, variable):
        """Creates a radiobutton widget.
        """
        frame = Frame(root)
        Label(frame, text=label).pack(side=LEFT)
        for text,value in (buttons):
            Radiobutton(frame, text=text, variable=variable, value=value).pack(side=LEFT)
        return frame
        
    def create_listbox(self, root, label, items, lines, function):
        """Creates a listbox widget.
        """
        frame = Frame(root)
        Label(frame, text=label).pack(side=LEFT)
        item_var = StringVar(value=items)
        listbox = Listbox(frame, listvariable=item_var, height=lines)
        for item, index in enumerate(items):
            listbox.insert(item, index)
        scrollbar = Scrollbar(frame)
        scrollbar.pack(side=RIGHT, fill=BOTH)
        listbox.config(yscrollcommand=scrollbar.set)
        scrollbar.config(command=listbox.yview)
        listbox.pack(side=LEFT,expand=YES, fill=X)
        listbox.bind('<<ListboxSelect>>', self.update_enzyme)
        return frame
        
    def create_combobox(self, root, label, variable, items):
        """Creates a ttk combobox widget.
        """
        frame = Frame(root)
        Label(frame, text=label).pack(side=LEFT)
        combobox = ttk.Combobox(frame, textvariable=variable)
        combobox['values'] = items
        combobox.pack(side=LEFT,expand=YES, fill=X)
        return frame, combobox
        
    #Helper Methods
    def load_database(self):
        """Maybe just handles a cancel button click?
        """
        database_name = PAW_lib.get_file(self.def_location, self.extension_list, 'Select a FASTA database file')
        if not database_name:
            self.quit_gui()
        self.database.set(database_name)
        self.filename = database_name
        return
          
    def update_enzyme(self, evt):
        """Grabs value of search_enzyme when user click on option in combo box
        """
        w = evt.widget
        self.enzyme_number = int(w.curselection()[0])
        enzyme = w.get(self.enzyme_number)
        print(enzyme)
        return
        
    def callback(self, *args):
        """Callback function creates two lists of all values grabbed from variable mods frame (deltamass and residue)
        """
        self.deltamass_values = []
        self.residue_values = []
        for i in range(self.max_num_mods):
            try:
                self.deltamass_values.append(self.num_vars_deltamass[i].get())
            except ValueError:
                pass
            self.residue_values.append(self.num_vars_residue[i].get())
            
    def check_for_blank_mods(self):
        """Catches if user inputs a variable mod into a slot number higher than the next available one, ie: 
        user puts '16.0 M'  into variable mod 6, when variable mod 1 is still '0.0 X'
        """
        # this does not change the GUI variables but copies them into a companion structure
        self.vars_deltamass = [0.0 for x in range(self.max_num_mods)]
        self.vars_residue = ['' for x in range(self.max_num_mods)]
        mod_counter = 0
        for mod in range(self.max_num_mods):
            if self.num_vars_deltamass[mod].get() > 0.0:
                self.vars_deltamass[mod_counter] = self.num_vars_deltamass[mod].get()
                self.vars_residue[mod_counter] = self.num_vars_residue[mod].get()
                mod_counter += 1
        blank_mod = self.max_num_mods - mod_counter
        for i in range(blank_mod):
            offset = i + mod_counter
            self.vars_deltamass[offset] = 0.0
            self.vars_residue[offset] = 'X'
            
    def update_variable_mods(self, row):
        """Gets the data from the copies of the GUI entry data.
        """
        residue = self.vars_residue[row]
        mass = self.vars_deltamass[row]
        if residue == 'N-term' and mass != 0.0:
            target_variable_mod = '%0.4f %s %s' % (mass, 'n', '0 3 -1 0 0')
        elif residue == 'C-term' and mass != 0.0:
            target_variable_mod = '%0.4f %s %s' % (mass, 'c', '0 3 -1 1 0')
        else:
            target_variable_mod = '%0.4f %s %s' % (mass, residue.upper(), '0 3 -1 0 0')
        return target_variable_mod
            
    def get_ms2_folder(self):
        """get the folder where MS2 files are located"""
        location = os.getcwd()
        self.ms2_folder = PAW_lib.get_folder(location, 'Select folder with MS2 files')
        return

    def save_settings(self):
        """Creates an updated Comet parameters file based on the user input.
        """
        self.check_for_blank_mods()
        user_selected_params = {
                'database_name': (self.database.get(), ''),
                'peptide_mass_tolerance': (self.pep_mass_tol.get(), ''),
                'fragment_bin_tol': (self.frag_bin_tol.get(), '# binning to use on fragment ions'),
                'peptide_mass_units': (self.pep_mass_units.get(), '# 0=amu, 1=mmu, 2=ppm'),
                'mass_type_parent': (self.mass_type_par.get(), '# 0=average masses, 1=monoisotopic masses'),
                'mass_type_fragment': (self.mass_type_frag.get(), '# 0=average masses, 1=monoisotopic masses'),
                'use_A_ions': (self.use_A_ions.get(), ''),
                'use_B_ions': (self.use_B_ions.get(), ''),
                'use_C_ions': (self.use_C_ions.get(), ''),
                'use_X_ions': (self.use_X_ions.get(), ''),
                'use_Y_ions': (self.use_Y_ions.get(), ''),
                'use_Z_ions': (self.use_Z_ions.get(), ''),
                'use_NL_ions': (self.use_NL_ions.get(), '# 0=no, 1=yes to consider NH3/H2O neutral loss peaks'),
                'search_enzyme_number': (self.enzyme.get(), '# choose from list at end of this params file'),
                'variable_mod01': (self.update_variable_mods(0), ''), 
                'variable_mod02': (self.update_variable_mods(1), ''),
                'variable_mod03': (self.update_variable_mods(2), ''),
                'variable_mod04': (self.update_variable_mods(3), ''),
                'variable_mod05': (self.update_variable_mods(4), ''),
                'variable_mod06': (self.update_variable_mods(5), ''),
                'variable_mod07': (self.update_variable_mods(6), ''),
                'variable_mod08': (self.update_variable_mods(7), ''),
                'variable_mod09': (self.update_variable_mods(8), ''),
                'num_threads': ('20', '# 0=poll CPU to set num threads; else specify num threads directly (max 64)'),
                'output_sqtfile': ('1', '# 0=no, 1=yes  write sqt file'),
                'output_pepxmlfile':  ('0', '# 0=no, 1=yes  write pep.xml file'),
                'num_output_lines': ('12', '# num peptide results to show')
            }

        # see if database can be found
        if not os.path.exists(user_selected_params['database_name'][0]):
            self.load_database()
            user_selected_params['database_name'] = (self.database.get(), '')
            
        # see if any static mods were changed
        if self.static_list:
            for key, value, comment in self.static_list:
                if comment:
                    user_selected_params[key] = (value, '# ' + comment)
                else:
                    user_selected_params[key] = (value, '')

        # get the folder where MS2 files are located (where to write params file)
        if not self.ms2_folder:
            self.get_ms2_folder()
            if not self.ms2_folder:
                return
            self.params_filename = os.path.join(self.ms2_folder, 'comet.params')
        
        with open(self.params_filename, 'w') as f:
            for line in comet_default_params.splitlines():              
                key = line.split('=')[0].strip()
                if key in user_selected_params: # come back and do string formatting
                    if key == 'peptide_mass_tolerance':
                        f.write(str(key) + ' = ' + '{0:.2f}'.format(self.pep_mass_tol.get()) + '\n')
                    elif key == 'fragment_bin_tol':
                        f.write('%s = %0.4f%s%s\n' % (key, user_selected_params[key][0], (30-len(key))*' ', user_selected_params[key][1]))
                    elif key.startswith('search_'):
                        idx = self.search_enzyme_list.index(user_selected_params[key][0])
                        user_selected_params[key] = (idx, user_selected_params[key][1])
                        pad = (36 - len(key) - len(str(user_selected_params[key][0]))) * ' '
                        f.write('%s = %s%s%s\n' % (key, user_selected_params[key][0], pad, user_selected_params[key][1]))                        
                    elif key.startswith('use_'):
                        if user_selected_params[key][0] == True:
                            user_selected_params[key] = ('1', user_selected_params[key][1])
                        else:
                            user_selected_params[key] = ('0', user_selected_params[key][1])
                        pad = (36 - len(key) - len(user_selected_params[key][0])) * ' '
                        f.write('%s = %s%s%s\n' % (key, user_selected_params[key][0], pad, user_selected_params[key][1]))
                    elif key.startswith('variable_'):
                        f.write('%s = %s \n' % (key, user_selected_params[key][0]))
                    elif key.startswith('add_'):
                        f.write('%s = %0.4f%s%s\n' % (key, user_selected_params[key][0], (30-len(key))*' ', user_selected_params[key][1]))
                    else:
                        pad = (36 - len(key) - len(str(user_selected_params[key][0]))) * ' '
                        f.write('%s = %s%s%s\n' % (key, user_selected_params[key][0], pad, user_selected_params[key][1]))
                else:
                    f.write(line + '\n')
#        self.run_comet['state'] = NORMAL
        return  
                                                                    
    def run_comet(self):
        """This executes Comet using the user-created params file. Assumes that
        the input file format is MS2 (output is SQT) and that "Comet" is defined
        as an executable command.
        """
        import glob
        import sqt_converter          

        # create a status bar for search progress
        self.progressbar = ttk.Progressbar(self.root)
        self.progressbar.pack(expand=Y, fill=X)
        self.progresstext = Label(self.root, text='Search progress')
        self.progresstext.pack()
        self.progressbar.update()

        # check if ms2 folder is set
        if not self.ms2_folder:
            self.get_ms2_folder()

        # check if comet.params exists
        self.filename = os.path.join(self.ms2_folder, 'comet.params')
        if not os.path.exists(self.filename):
            self.progresstext.configure(text='WARNING: comet.params file not found!!!')
            return

        os.chdir(self.ms2_folder)
        
        step = 100/len(glob.glob('*.ms2'))  # step for progressbar
        for ms2_file in glob.glob('*.ms2'):
            
            # update status bar
            self.progresstext.configure(text='Searching: %s' % (ms2_file,))
            self.progressbar.update()
            self.progressbar.step(step)
            self.progressbar.update()
            
            # run on each MS2 file with a wait for completion
            os.system('START "Comet" /WAIT /MIN /LOW CMD /C COMET2016 ' + ms2_file)

        # this creates top-hit TXT files from the SQT files after Comet has finished
        self.progresstext.configure(text='Searches completed. Starting TXT creation.')
##        self.progressbar.update()
##        self.progressbar.step(step)
        sqt_converter.main(os.path.dirname(self.filename), overwrite=True)
        self.progresstext.configure(text='Conversions completed. Quit when ready...')

    def change_static(self):
        """Creates a window to view/change the static modifications.
        """
        staticMods(self.root, comet_default_params, self.static_list)
        return

    def quit_gui(self):
        """Quits the GUI.
        """
        self.root.withdraw()
        self.root.update_idletasks()
        self.root.destroy()
        sys.exit()
                                                                                                                                                                                                       
    # Methods to create sections of GUI
    def create_dir_frame(self):
        """Lets the user browse to the FASTA database location
        """
        dir_frame = ttk.Labelframe(self.root, text='Database:')
        dir_frame.pack(fill=X, expand=YES, padx=5, pady=5)

        #Variables
        self.database=StringVar()
        self.def_location = r'C:\\'
        self.extension_list = [('FASTA File','*.fasta')]
        
        #Defaults
        self.database.set('Select a database (click the button on the left)')

        #Creation
        ttk.Button(dir_frame, text='Database', command=self.load_database).pack(side=LEFT)
        db_entry = ttk.Entry(dir_frame, textvariable=self.database)
        db_entry.pack(side=LEFT, fill=X, expand=YES)
        return
        
    def create_masses_frame(self):
        """Lets the user change parent ion mass tolerance, type, etc.
        """
        masses_frame = ttk.Labelframe(self.root, text='Mass parameters:')
        masses_frame.pack(fill=X, expand=YES, padx=5, pady=5)
        
        #Variables
        self.pep_mass_tol = DoubleVar()
        self.frag_bin_tol = DoubleVar()
        self.pep_mass_units = IntVar()
        self.mass_type_par = IntVar()
        self.mass_type_frag = IntVar()        
        
        #Creation
        self.create_entry(masses_frame, 'Peptide Mass Tolerance: ', self.pep_mass_tol).pack(fill=X, expand=YES)
        self.create_radiobuttons(masses_frame, 'Peptide Mass Units: ',
                                 [('AMU',0),('MMU',1),('PPM',2)], self.pep_mass_units).pack(fill=X, expand=YES)
        self.create_radiobuttons(masses_frame, 'Parent Ion Mass Type: ',
                                 [('Average',0),('Monoisotopic',1)], self.mass_type_par).pack(fill=X, expand=YES)
        self.create_entry(masses_frame, 'Fragment Bin Tolerance (use 1.0005 for low-res IT) : ', self.frag_bin_tol).pack(fill=X, expand=YES)        
        
        #Set Defaults
        self.pep_mass_tol.set(1.25)
        self.frag_bin_tol.set(1.0005)
        self.pep_mass_units.set(0)
        self.mass_type_par.set(1)
        self.mass_type_frag.set(1)
        return
        
    def create_use_ion_frame(self):
        """Lets the user select the ion series to use in scoring.
        """
        use_ion_frame = ttk.Labelframe(self.root, text='Ion series:')
        use_ion_frame.pack(fill=X, expand=YES, padx=5, pady=5)
        
        #Variables
        self.use_A_ions = BooleanVar()
        self.use_B_ions = BooleanVar()
        self.use_C_ions = BooleanVar()
        self.use_X_ions = BooleanVar()
        self.use_Y_ions = BooleanVar()
        self.use_Z_ions = BooleanVar()
        self.use_NL_ions = BooleanVar()
        ions_list = [
        ('A ions', self.use_A_ions),
        ('B ions', self.use_B_ions),
        ('C ions', self.use_C_ions),
        ('X ions', self.use_X_ions),
        ('Y ions', self.use_Y_ions),
        ('Z ions', self.use_Z_ions),
        ('NL ions', self.use_NL_ions)
        ]
        
        #Creation
        self.create_checkboxes(use_ion_frame, 'Use: ', ions_list).pack(fill=X, expand=YES)
        
        #Set Defaults
        self.use_B_ions.set(True)
        self.use_Y_ions.set(True)
        self.use_NL_ions.set(True)
        return
                
    def create_search_enzyme_frame(self):
        """Lets user select digestion enzyme from pulldown menu.
        """
        enzyme_frame = ttk.Labelframe(self.root, text='Enzyme:')
        enzyme_frame.pack(fill=X, expand=YES, padx=5, pady=5)
        
        #Variables
        self.search_enzyme_list = [
        'No Enzyme', 
        'Trypsin',
        'Trypsin/P', 
        'Lys_C', 
        'Lys_N', 
        'Arg_C', 
        'Asp_N', 
        'CNBr',
        'Glu_C',
        'PepsinA',
        'Chymotrypsin']
        
        #Creation
        self.enzyme = StringVar()
        frame, combobox = self.create_combobox(enzyme_frame, 'Search Enzyme: ', self.enzyme,
                                               self.search_enzyme_list)
        frame.pack(fill=X, expand=YES, padx=5, pady=5)

        # set default
        combobox.set('Trypsin')
        return
        
    def create_variable_mods_frame(self):
        """Lets user specify up to 7 variable mods (last two used for n-term or c-term mods).
        """
        self.variable_mods_frame = ttk.Labelframe(self.root, text='Variable modifications:')
        self.variable_mods_frame.pack(fill=X, expand=YES, padx=5, pady=5)
        
        #Variables
        self.num_variable_mods = IntVar()
        self.variable_mod_binary = IntVar()
        self.variable_mod_maximum = IntVar()
        self.max_num_mods = 9
        self.num_vars_residue = {}
        self.num_vars_deltamass = {}

        
        #Creation
        self.top_label = ('Delta Mass', 'Residues')
        for label in range(1, self.max_num_mods+1):
            if label == 8:
                Label(self.variable_mods_frame, text='N-term Mod:').grid(column=0,row=label+1)
            elif label == 9:
                Label(self.variable_mods_frame, text='C-term Mod:').grid(column=0,row=label+1)
            elif not (4 <= label <= 7):
                Label(self.variable_mods_frame, text='Mod %i:'%label).grid(column=0,row=label+1)

        for label in range(len(self.top_label)):
            Label(self.variable_mods_frame, text=self.top_label[label]).grid(column=label+1,row=0)

        for i in range(self.max_num_mods):
            self.num_vars_residue[i] = StringVar()
            self.num_vars_deltamass[i] = DoubleVar()
            #Sneaky Defaults
            self.num_vars_deltamass[i].set(0.0)
            if i == 0:
                self.num_vars_residue[i].set('M')
                self.num_vars_deltamass[i].set(15.9949)
            elif i == 7:
                self.num_vars_residue[i].set('N-term')
            elif i == 8:
                self.num_vars_residue[i].set('C-term')
            else:
                self.num_vars_residue[i].set('X')
            
            for j in range(len(self.top_label)):
                if 2 < i < 7:
                    continue
                if j == 0:
                    Entry(self.variable_mods_frame, textvariable=self.num_vars_deltamass[i]).grid(row=i+2,column=j+1)
                    self.num_vars_deltamass[i].trace('w', self.callback)
                elif j == 1:
                    if i < 7:
                        Entry(self.variable_mods_frame, textvariable=self.num_vars_residue[i]).grid(row=i+2,column=j+1)
                    else:
                        Entry(self.variable_mods_frame, textvariable=self.num_vars_residue[i], state='readonly').grid(row=i+2,column=j+1)
                    self.num_vars_residue[i].trace('w', self.callback)
                else:
                    pass       
    
if __name__ == '__main__':
    comet = CometGUI()
