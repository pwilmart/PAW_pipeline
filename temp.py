import os, sys
from pprint import pprint
from collections import OrderedDict

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

params = CometParams()
params.load_from_folder("E:/PXD002875_CarbonSources/msn_files")
params._snoop()
