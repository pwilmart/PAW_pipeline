"""MSConvert_GUI.py, version 1.07.
Simple GUI interface for converting RAW files into MSn files.
Uses MSConvert ultility from Protewizard toolkit to convert RAW files.
Creates MS2-format files for searches and extracts reporter ion
intensity peak areas from MS3 scan in multi-notch TMT runs.

written by Phil Wilmarth, OHSU, 2014, 2016, 2018.

The MIT License (MIT)

Copyright (c) 2018 Phillip A. Wilmarth and OHSU

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

Revision history:
version 1.01:
    first version late summer/fall 2014 - PW
version 1.02:
    better support for multi-notch TMT data - PW
version 1.03:
    better output options - 6/2016 PW
version 1.04:
    rewritten for Python 3
    directly integrates MS3 reporter ion peaks - 7/2016 PW
version 1.05:
    now does all 11 TMT tags (131C added)
    renamed tags so they sort in mass order (128_N before 128C) - PW 11/2017
version 1.06:
    separate reporter ion files (one for each RAW) -PW 3/2018
version 1.07:
    support for MS2 or MS3 reporter ions -PW 20180711
    PAW_tmt.txt files are moved to location with MSn files -PW 20180816
version 1.08:
    support for real time search SPS MS3 data -PW 20200520
    support for 16-plex TMT reagents -PW 20200630

"""
# global imports
from tkinter import *
import tkinter.ttk as ttk
# from tkinter import filedialog    # PAW_lib might take care of this

import os
import sys
import gzip
import time
import glob
import shutil

import numpy as np
import pandas as pd
from scipy.integrate import trapz

import PAW_lib

# flag for computing centroids (m/z positions)
COMPUTE_CENTROIDS = False

# flag for computing peak areas (trapezoidal integrations are slow)
COMPUTE_AREAS = False

class Spectra:
    """Container for Spectrum objects."""
    def __init__(self, count=15, intensity=500, msn_level=2, header_block=None, parent=None):
        self.count = count              # minimum peak count
        self.intensity = intensity      # minimum total absolute intensity
        self.parent = parent            # parent object
        self.version = 1.04             # this program version
        self.xcalibur_version = None    # Xcalibur version
        self.msconvert_version = None   # Proteowizard version
        self.instrument = None          # Instrument type
        self.instrument_SN = None       # Instrument serial number
        self.msn = msn_level            # msn level being extracted
        self.spectra = []               # list of Spectrum objects
        self.mz_zero = 0                # how many scans had no m/z value
        self.z_zero = 0                 # how many scans had no charge state
        self.freq = {}                  # QC check on decimal places in m/z values
        self.get_versions(header_block)

    def get_versions(self, header_block):
        """Get Proteowizard version numbers."""
        if header_block:
            for i, line in enumerate(header_block):
                if line.startswith('cvParam: Xcalibur'):
                    self.xcalibur_version = header_block[i-1].split('version:')[1].strip()
                if line.startswith('cvParam: ProteoWizard software'):
                    self.msconvert_version = header_block[i-1].split('version:')[1].strip()
                if line.startswith('cvParam: instrument serial number'):
                    try:
                        self.instrument_SN = line.split('number,')[1].strip()
                    except IndexError:
                        self.instrument_SN = '0'
                    try:
                        self.instrument = header_block[i-1].split('cvParam:')[1].strip()
                    except IndexError:
                        self.instrument = 'unknown'
        return # nothing happens if there is no header_block

    def add(self, spectrum):
        """Add a spectrum if it meets the criteria."""
        if (spectrum.intensity >= self.intensity and
            len(spectrum.mz_array) >= self.count and
            spectrum.msn == self.msn):
            if spectrum.m_over_z == 0.0:
                spectrum.m_over_z = spectrum.cv_moverz
            self.check_spectrum(spectrum)
            if spectrum.m_over_z == 0.0:
                self.substitute_mz(spectrum)
            if spectrum.charge == [] or spectrum.charge[0] == 0:
                self.determine_charge(spectrum)
            self.spectra.append(spectrum)

    def determine_charge(self, spectrum):
        """Tests unknown charge state scans for 1+ or 2+/3+."""
        total_int = sum(spectrum.int_array)
        test_int = 0.0
        for i, mz in enumerate(spectrum.mz_array):
            if mz < spectrum.m_over_z:
                test_int += spectrum.int_array[i]
        if test_int >= 0.90 * total_int:
            spectrum.charge = [1]
        else:
            spectrum.charge = [2, 3]

    def substitute_mz(self, spectrum):
        """Sometimes the m/z value in userParam: line is zero so get value from scan header string."""
        m_over_z = spectrum.filter_string.split('string,')[-1]
        spectrum.m_over_z = float(m_over_z.split('@')[0].split()[-1])

    def check_spectrum(self, spectrum):
        """Checks spectrum for missing charge state, missing m/z, and checks accuracy of m/z."""
        if spectrum.m_over_z == 0.0:
            self.mz_zero += 1
        if spectrum.charge  == [] or spectrum.charge[0] == 0:
            self.z_zero += 1
        places = len(spectrum.mzstr) - (spectrum.mzstr.index('.') + 1)
        self.freq[places] = self.freq.setdefault(places, 0) + 1

    def report(self):
        """Prints some summary stats."""
        print('......%s m/z values were zero and replaced by m/z in scan label' % self.mz_zero)
        print('......%s charge states were zero or missing' % self.z_zero)
        print('......m_over_z decimal places frequency counts:')
        for k in sorted(self.freq.keys()):
            print('.........%d places occurred %d times' % (k, self.freq[k]))

    def write_msn(self, msn_name, low=0.0, high=2000.0):
        """Write spectra in MSn format (MS2 and MS3 are same format)."""
        msn = open(msn_name, 'w')

        # write header lines
        msn.write('H\tCreationDate\t' + time.ctime() + '\n')
        msn.write('H\tInstrument\t' + self.instrument + '\n')
        msn.write('H\tInstrumentSerialNumber\t' + self.instrument_SN + '\n')
        msn.write('H\tExtractor\t' + 'Proteowizard + PAW\n')
        msn.write('H\tExtractorVersion\t1.03\n')
        msn.write('H\tExcaliburVersion\t' + self.xcalibur_version + '\n')
        msn.write('H\tMSConvertVersion\t' + self.msconvert_version + '\n')
        msn.write('H\tSourceFile\t' + os.path.basename(msn_name) + '.RAW\n')

        for spectrum in self.spectra:
            spectrum.write_msn_block(msn, low, high)

        msn.close()

    # end Spectra class

class Spectrum:
    def __init__(self, parent, block):
        """Parses info out of spectrum blocks."""
        # initializations
        self.parent = parent    # pointer to Spectra (parent) object
        get_mz = True           # flag
        self.mz_array = []      # holds m/z values
        self.int_array = []     # holds intensity values
        self.charge = []        # can have multiple charge states for low res data
        self.scan = None        # scan number
        self.cv_moverz = 0.0    # also m/z of precursor?
        self.m_over_z = 0.0     # m/z value of precursor
        self.mzstr = ''         # scan header string (has useful data)

        # parse out the desired items
        for line in block:
            if line.startswith('id:') and not self.scan:
                self.scan = int(line.split('scan=')[-1])
            if line.startswith('defaultArrayLength:'):
                self.count = int(line.split()[1]) # length of data arrays
            if line.startswith('cvParam:'):
                if 'ms level,' in line:
                    self.msn = int(line.split('ms level,')[-1])
                if 'total ion current,' in line:
                    self.intensity = float(line.split('total ion current,')[-1])
                if 'spectrum title,' in line:
                    self.title = line.split('spectrum title,')[-1].split()[0]
                if 'scan start time,' in line:
                    self.rt_min = float(line.split('time,')[-1].split(',')[0])
                if 'filter string,' in line:
                    self.filter_string = line.split('filter string, ')[-1]
                if 'ion injection time,' in line:
                    self.inj_time_ms = float(line.split('ion injection time,')[-1].split(',')[0])
                if 'selected ion m/z,' in line: # this seems to be rounded to 2 decimal places
                    self.cv_moverz = float(line.split('selected ion m/z,')[-1].split(',')[0])
                if 'charge state,' in line:
                    self.charge.append(int(line.split('charge state,')[-1]))
            if line.startswith('userParam:'):
                if 'Monoisotopic M/Z:,' in line:
                    self.m_over_z = float(line.split(',')[1])
                    self.mzstr = line.split(',')[1].strip()
                    if '.' not in self.mzstr:
                        self.mzstr += '.'
            if line.startswith('binary:'):
                if get_mz:
                    self.mz_array = [float(x) for x in line.split(']')[1].split()]
                    get_mz = False
                else:
                    self.int_array = [float(x) for x in line.split(']')[1].split()]

        # check that array lengths match
        if len(self.mz_array) != self.count or len(self.int_array) != self.count \
           or len(self.mz_array) != len(self.int_array):
            print('...WARNING: non-equal ion and intensity counts!')

        return

    def console_dump(self):
        """Prints parsed values to console."""
        print('\nScan: %d, MSn: %d, Int: %0.0f, Title: %s, RT: %0.2f, Inj: %0.2f, M/Z: %0.2f, Charge: %s' %
              (self.scan, self.msn, self.intensity, self.title, self.rt_min,
               self.inj_time_ms, self.m_over_z, self.charge))
        print('number of data points:', len(self.mz_array), len(self.int_array))

    def write_msn_block(self, fout, low, high):
        """Writes self as one spectrum in an MSn file."""
        fout.write('S\t%s\t%s\t%0.5f\n' % (self.scan, self.scan, self.m_over_z))
        fout.write('I\tText\t%s\n' % self.filter_string)
        fout.write('I\tMSn\t%d\n' % self.msn)
        fout.write('I\tRTime\t%0.2f\n' % self.rt_min)
        fout.write('I\tTIC\t%0.1f\n' % self.intensity)
        for z in self.charge:
            MHplus = (self.m_over_z * float(z)) - (z-1)*(1.007825)
            fout.write('Z\t%d\t%0.5f\n' % (z, MHplus))
        for i in range(self.count):
            if low <= self.mz_array[i] <= high:
                fout.write('%0.4f %0.2f\n' % (self.mz_array[i], self.int_array[i]))
        return

    # end Spectrum class

class Reporter_ion:
    """Data container for TMT peak centroids, areas, and heights."""
    def __init__(self, lc_name, ms2_scan, ms3_scan, centroids, areas, heights):
        self.lc_name = lc_name      # LC run name
        self.ms2_scan = ms2_scan    # MS2 scan number
        self.ms3_scan = ms3_scan    # MS3 scan number (if applicable)
        self.centroids = centroids  # list of peak centroid positions
        self.areas = areas          # list of trapezoidal peak integrations
        self.heights = heights      # list of peak heights (intensities)
        return

    # NOTE: channel lables and lengths of centroids, areas, and heights should agree
    # these are set by the caller
    def write_header(self, channels, fout):
        """Writes a tab-delimited header line to 'fout'."""
        headers = (['lc_name', 'ms2_scan', 'ms3_scan'] +
                   [('cent_' + x) for x in channels] +
                   [('area_' + x) for x in channels] +
                   [('height_' + x) for x in channels])
        fout.write('\t'.join(headers) + '\n')
        return

    def write_row(self, fout):
        """Writes a tab-delimited data row to 'fout'.
        """
        row = ([self.lc_name, self.ms2_scan, self.ms3_scan] +
               [round(x, 6) for x in self.centroids] +
               [round(x, 6) for x in self.areas] +
               [round(x, 6) for x in self.heights])
        fout.write('\t'.join([str(x) for x in row]) + '\n')
        return

    # end Reporter_ions class

class MSConvertGUI:
    """Main GUI for running MSConvert in batch mode to create MSn files."""

    def __init__(self):
        """Constructor. Uses a collection of frames for related parameters."""
        self.root = Tk()
        self.root.minsize(500, 250)
        self.root.title('MSConvert MSn creator v1.06')
        self.root.protocol('WM_DELETE_WINDOW', self.quit_gui)

        # some globals
        self.msn_total = 0       # number of MSn scans written to MSn files
        self.reporter_total = 0  # number of reporter ion sets written to text files

        # maybe get the window to the top...
        self.root.attributes('-topmost', 1)
        self.root.attributes('-topmost', 0)

        self.raw_name_list = None
        self.progressbar = None

        # create GUI
        self.create_files_frame()
        self.create_defaults_frame()

        # buttons at bottom of GUI
        Button(self.root, text='Start conversion', command=self.start_processing).pack(pady=2)
        Button(self.root, text='Quit', command=self.quit_gui).pack(pady=2)

        # enter mainloop
        self.root.mainloop()
        return

    def set_tmt_plex(self):
        """Sets the TMt plex when Start conversion button is clicked."""
        # load TMT labels and masses (all 11 - others kits are subsets)
        # assumes we have high resolution data where N- and C-forms are resolved
        tmt_plex = self.tmt_plex.get()
        if (tmt_plex == 0) or (tmt_plex == 1):       # 6-plex
            self.tmt_tuples = [('126C', 126.1278),
                               ('127_N', 127.1248),
                               ('128C', 128.1344),
                               ('129_N', 129.1315),
                               ('130C', 130.1411),
                               ('131_N', 131.1382)]
        elif tmt_plex == 2:     # 10-plex
            self.tmt_tuples = [('126C', 126.1278),
                               ('127_N', 127.1248),
                               ('127C', 127.1311),
                               ('128_N', 128.1281),
                               ('128C', 128.1344),
                               ('129_N', 129.1315),
                               ('129C', 129.1378),
                               ('130_N', 130.1348),
                               ('130C', 130.1411),
                               ('131_N', 131.1382)]
        elif tmt_plex == 3:     # 11-plex
            self.tmt_tuples = [('126C', 126.1278),
                               ('127_N', 127.1248),
                               ('127C', 127.1311),
                               ('128_N', 128.1281),
                               ('128C', 128.1344),
                               ('129_N', 129.1315),
                               ('129C', 129.1378),
                               ('130_N', 130.1348),
                               ('130C', 130.1411),
                               ('131_N', 131.1382),
                               ('131C', 131.1445)]
        elif tmt_plex == 4:     # 16-plex
            self.tmt_tuples = [('126C', 126.1278),
                               ('127_N', 127.1248),
                               ('127C', 127.1311),
                               ('128_N', 128.1281),
                               ('128C', 128.1344),
                               ('129_N', 129.1315),
                               ('129C', 129.1378),
                               ('130_N', 130.1348),
                               ('130C', 130.1411),
                               ('131_N', 131.1382),
                               ('131C', 131.1445),
                               ('132_N', 132.1415),
                               ('132C', 132.1479),
                               ('133_N', 133.1449),
                               ('133C', 133.1512),
                               ('134_N', 134.1482)]
        elif tmt_plex == 5:     # 18-plex
            self.tmt_tuples = [('126C', 126.1278),
                               ('127_N', 127.1248),
                               ('127C', 127.1311),
                               ('128_N', 128.1281),
                               ('128C', 128.1344),
                               ('129_N', 129.1315),
                               ('129C', 129.1378),
                               ('130_N', 130.1348),
                               ('130C', 130.1411),
                               ('131_N', 131.1382),
                               ('131C', 131.1445),
                               ('132_N', 132.1415),
                               ('132C', 132.1479),
                               ('133_N', 133.1449),
                               ('133C', 133.1512),
                               ('134_N', 134.1482),
                               ('134C', 134.1546),
                               ('135_N', 135.1516)]
        self.tmt_labels = [x[0] for x in self.tmt_tuples]
        self.tmt_masses = [x[1] for x in self.tmt_tuples]
        self.zeroes = [0.0 for x in self.tmt_tuples]

        # set up m/z windows for TMT reporter ions
        self.set_windows()
        return

    # functions to help create widgets
    def create_entry(self, parent, label, variable):
        """Creates a text entry widget."""
        frame = Frame(parent)
        Label(frame, text=label).pack(side=LEFT)
        Entry(frame, textvariable=variable).pack(side=LEFT)
        return frame

    def create_checkboxes(self, parent, label, ions_list):
        """Creates a checkbox widget."""
        frame = Frame(parent)
        Label(frame, text=label).pack(side=LEFT)
        for button, variable in ions_list:
            Checkbutton(frame, text=button, variable=variable).pack(side=LEFT)
        return frame

    def create_radiobuttons(self, parent, label, buttons, variable):
        """Creates a radiobutton widget."""
        frame = Frame(parent)
        Label(frame, text=label).pack(side=LEFT)
        for text, value in (buttons):
            Radiobutton(frame, text=text, variable=variable, value=value).pack(side=LEFT)
        return frame

    def create_progressbar(self, parent):
        """create a status bar for conversion progress."""
        self.progressbar = ttk.Progressbar(parent, maximum=101)
        self.progressbar.pack(expand=Y, fill=X)
        self.progresstext = Label(self.root, text='Conversion progress')
        self.progresstext.pack()
        self.progressbar.update()
        return

    # main button methods
    def load_raw_files(self):
        """Allows user to select one or more RAW files. Sets folder path attribute and populates the
        raw_files_list attribute.
        """
        self.raw_name_list = PAW_lib.get_files(self.def_location, self.extension_list, 'Select RAW files(s)')
        if self.raw_name_list:
            self.raw_path = os.path.dirname(self.raw_name_list[0])
            self.raw_path_view.set(' %s   (%i RAW files selected)' % (self.raw_path, len(self.raw_name_list)))
        return

    def start_processing(self):
        """This executes MSConvert using the user-specified options. Then creates the MSn files
        from the gzipped text files.
        """
        starting_time = time.time()
        if not self.raw_name_list:
            self.load_raw_files()

        # set the TMT plex
        self.set_tmt_plex()
        
        # set up log file
        log_name = 'MSConvert_GUI_log.txt'
##        log_name = 'MSConvert_GUI_%s_log.txt' % time.time() # this add a time stamp - will have multiple log files
        log_obj = open(os.path.join(self.raw_path, log_name), mode='at')
        self.log_obj = [None, log_obj]
        for obj in self.log_obj:
            print('...starting conversions at:', time.ctime(), file=obj)

        # build the MSConvert command options
        centroid = [[' --filter "peakPicking true 2-3"', ' --filter "peakPicking true 2"'],
                    [' --filter "peakPicking true 3"', '']]
        level = [' --filter "msLevel 2"',
                 ' --filter "msLevel 3"',
                 ' --filter "msLevel 2"',
                 ' --filter "msLevel 2-3"']
        idx = self.ms2_centroid.get()
        idy = self.ms3_centroid.get()
        if self.msn_level.get() == 0 or self.msn_level.get() == 2: # if ms2 only, can't centroid ms3
            idy = 1
        elif self.msn_level.get() == 1: # if ms3 only, can't centroid ms2
            idx = 1
        centroid_picked = centroid[idx][idy]
        level_picked = level[self.msn_level.get()]

        # create a status bar for conversion progress
        if self.progressbar is None:
            self.create_progressbar(self.root)
        step = 100/(2*len(self.raw_name_list))  # progressbar steps for conversion to text and conversion to MSn

        # set current dir to raw file folder
        dir_loc = self.raw_path
        os.chdir(dir_loc)

        for raw_name in self.raw_name_list:
            # update status bar
            self.progresstext.configure(text='Converting RAW to Text: %s' % (os.path.basename(raw_name),))
            self.progressbar.update()

            # clear any data structures that reset with each RAW file
            self.tmt_data = []
            lc_name = os.path.splitext(os.path.basename(raw_name))[0]
            if raw_name.lower().endswith('.raw'):
                msconvert_name = raw_name[:-4] + '.txt.gz'

            # call MSConvert
            quoted_raw_name = '"%s"' % raw_name  # in case the path has spaces
            command_line = quoted_raw_name + centroid_picked + level_picked + ' --text --gzip'
            if not os.path.exists(msconvert_name):
                for obj in self.log_obj:
                    print('MSConvert ' + command_line, file=obj)
                os.system('START "MSConvert" /WAIT /MIN /LOW CMD /C MSConvert ' + command_line)
            else:
                for obj in self.log_obj:
                    print('...Skipping conversion of:', lc_name, file=obj)
            self.progressbar.step(step)
            self.progressbar.update()

            # compute the expected MSConvert filename from RAW name
            txt_name_short = os.path.splitext(os.path.basename(raw_name))[0] + '.txt.gz'
            self.txt_name = os.path.join(os.path.dirname(raw_name), txt_name_short) # add full path

            # if RAW file is corrupt, MSConvert will not create a ".txt.gz" file
            if os.path.exists(self.txt_name):
                self.progresstext.configure(text='Converting %s to MSn file' % (txt_name_short,))
                self.progressbar.update()

                # create the MSn file first
                if self.msn_level.get() == 1:
                    msn_level = 3
                else:
                    msn_level = 2
                if self.msn_level.get() == 2:   # MS2 TMT experiment
                    reporter_ions = True
                else:
                    reporter_ions = False
                self.process_msn_level(lc_name, msn_level, self.ion_count.get(), self.min_intensity.get(), reporter_ions)
                if self.msn_level.get() == 3:
                    self.process_ms3_reporter_ions(lc_name)
                    self.progresstext.configure(text='Processing reporter ions')
                self.progressbar.step(step)
                self.progressbar.update()

                # save the TMT results (one quant file for each RAW file)
                if self.tmt_data:
                    fout = open(os.path.join(os.path.dirname(self.raw_name_list[0]), lc_name + '.PAW_tmt.txt'), 'w', newline=None)
                    self.tmt_data[0].write_header(self.tmt_labels, fout) # header line
                    for result in self.tmt_data:
                        result.write_row(fout) # data line
                    fout.close()
            else:
                for obj in self.log_obj:
                    print('...WARNING: possible corrupt file:', lc_name, file=obj)

        # move the MSn and PAW_tmt files into separate folder
        ending_time = time.time()
        self.progresstext.configure(text='Moving search and quant files')
        self.progressbar.update()
        self.move_files()

        # write some log data to console
        for obj in self.log_obj:
            print('...Conversions ended at:', time.ctime(), file=obj)
            print('...Conversions took', ending_time - starting_time, 'seconds', file=obj)
            print("\nThere were %s scans written to MSn files and %s reporter ion regions processed." %
                  (self.msn_total, self.reporter_total), file=obj)

        # update status line when done
        self.progresstext.configure(text='Conversions completed. Quit when ready...')
        self.progressbar.update()
        self.progressbar.step(step)
        self.progressbar.update()

        # try and close log files
        for obj in self.log_obj:
            try:
                obj.close()
            except:
                pass

    def process_ms3_reporter_ions(self, lc_name):
        """Parses ms2 and ms3 spectrum blocks; gets scan numbers, and reporter ion data.
        """
        ms2_count = 0   # MS2 scan counter
        ms3_count = 0   # MS3 scan counter

        # parse files to get ms2 scan numbers, ms3 scan numbers, and data
        ms_dict = {}            # used to link MS3 scan numbers to ms2 scan numbers
        spectrum_flag = False   # limits parsing to spectrum blocks
        msn_level = None        # the MSn level of the scan
#        mz_arr_flag = False     # True if data array is m/z values
        ms1_prev = 0            # previous (maybe current?) MS1 scan number
#        moverz_key = ''         # relevant string from scan header line to link MS2 and MS3
        dict_list = []          # keeps track of two scan cycles to support RTS data
        buff = []

        # read lines in MSConvert file
        for k, line in enumerate(gzip.open(self.txt_name, 'rt')):
            line = line.strip()    # remove leading, trailing white space

           # look for each spectrum block
            if line.startswith('spectrum:') or line.startswith('chromatogramList'):
                spectrum_flag = True
                if buff: # parse the previous spectrum block
                    for buff_line in buff:
                        if buff_line.startswith('cvParam: ms level,'):   # get MSn level (2 or 3)
                            msn_level = int(buff_line.split()[-1])

                            # MS2 scan parsing
                            if msn_level == 2:
                                ms2_count += 1
                                if (ms2_count % 1000) == 0:
                                    for obj in self.log_obj:
                                        print('......%d MS2 scans processed...' % ms2_count, file=obj)
                                ms1_scan = self.parse_ms2_scan(buff, ms_dict)

                                # look for next instrument cycle (a new MS1 scan)
                                if ms1_scan != ms1_prev:
                                    dict_list.append(ms_dict)
                                    ms_dict = {}
                                    ms1_prev = ms1_scan

                            # MS3 scan parsing
                            elif msn_level == 3:
                                ms3_count += 1
                                self.parse_ms3_scan(buff, ms_dict, dict_list, lc_name)

                            break

                # reset buffer
                buff = []

            # spectrum_flag skips header lines until first spectrum block. Stops at chromato info.
            if line.startswith('chromatogramList'):
                spectrum_flag = False

            # save lines in buffer if inside a spectrum block
            if spectrum_flag:
                buff.append(line)

        # update total MS3 counter
        self.reporter_total += ms3_count

    def parse_ms2_scan(self, buffer, ms_dict):
        """Parses an MS2 scan block.
        """
        # read lines
        for line in buffer:

            # get the scan number
            if line.startswith('id:') and 'scan=' in line:
                scan_num = int(line.split()[-1].split('=')[-1])

            # get MS1 scan number
            if line.startswith('spectrumRef: '):
                ms1_scan = int(line.split()[-1].split('=')[-1])

            # get MS2 dissociation key (m/z value) and link to MS2 scan number
            if line.startswith('cvParam: filter string'):
                if '@cid3' in line:
                    moverz_key = line.split('@cid3')[0].split()[-1]
                elif '@hcd3' in line:
                    moverz_key = line.split('@hcd3')[0].split()[-1]
                else:
                    for obj in self.log_obj:
                        print('WARNING: dissociation key (@cid or @hcd) not found', file=obj)
                    return
                ms_dict[moverz_key] = scan_num

        return ms1_scan

    def parse_ms3_scan(self, buffer, ms_dict, dict_list, lc_name):
        """Parses an MS2 scan block.
        """
        # read lines
        for line in buffer:

            # get the scan number
            if line.startswith('id:') and 'scan=' in line:  # get scan number
                scan_num = int(line.split()[-1].split('=')[-1])

            # get MS2 dissociation key (m/z value) accosicate with the MS3 scan
            if line.startswith('cvParam: filter string'):
                if '@cid3' in line:
                    moverz_key = line.split('@cid3')[0].split()[-1]
                elif '@hcd3' in line:
                    moverz_key = line.split('@hcd3')[0].split()[-1]
                else:
                    for obj in self.log_obj:
                        print('WARNING: dissociation key (@cid or @hcd) not found', file=obj)
                    return

                # try and lookup the MS2 scan number
                ms2_scan = None
                if moverz_key in ms_dict:
                    ms2_scan = ms_dict[moverz_key]
                else:
                    for i in range(10):
                        if moverz_key in dict_list[-(i+1)]:
                            ms2_scan = dict_list[-(i+1)][moverz_key]
                            break

                # might have small precursor mass error (mass correction?)
                if not ms2_scan:
                    print('starting fuzzy lookup :', moverz_key, scan_num)
                    fuzzy_list = []
                    lookup = float(moverz_key)
                    for key in ms_dict:                 # try current dictionary
                        if abs(lookup - float(key)) <= 0.02:
                            fuzzy_list.append(key)
                    if len(fuzzy_list) == 1:
                        ms2_scan = ms_dict[fuzzy_list[0]]
                        print('fuzzy matching key was:', fuzzy_list[0])
                    if not fuzzy_list:
                        for key in dict_list[-1]:       # try last dictionary in saved list
                            if abs(lookup - float(key)) <= 0.015:
                                fuzzy_list.append(key)
                        if len(fuzzy_list) == 1:
                            ms2_scan = dict_list[-1][fuzzy_list[0]]
                            print('fuzzy matching key was:', fuzzy_list[0])
                    if not fuzzy_list:
                        for key in dict_list[-2]:       # also try one farther back
                            if abs(lookup - float(key)) <= 0.015:
                                fuzzy_list.append(key)
                        if len(fuzzy_list) == 1:
                            ms2_scan = dict_list[-2][fuzzy_list[0]]
                            print('fuzzy matching key was:', fuzzy_list[0])

                # now it is time to give up
                if not ms2_scan:
                    print('*** MS2 lookup failed! scan:', scan_num, 'key:', moverz_key)
                    print(moverz_key in dict_list[-1])
                    print(ms_dict)
                    print(dict_list[-1])
                    return

            # get the scan data
            if line == 'cvParam: m/z array, m/z':
                mz_arr_flag = True
            if line.startswith('binary: ') and mz_arr_flag:
                mz_list = [float(x) for x in line.split()[2:]]
                mz_arr_flag = False
            elif line.startswith('binary: ') and not mz_arr_flag:
                int_list = [float(x) for x in line.split()[2:]]
                centroids, areas, heights = self.process_tmt_data(mz_list, int_list, scan_num)
                self.tmt_data.append(Reporter_ion(lc_name, ms2_scan, scan_num, centroids, areas, heights))
        return

    def process_tmt_data(self, mz_list, int_list, scan_num):
        """Computes peak centroid and integral within each TMT window.
        This will work for either MS2 reporter ions or MS3 reporter ions."""
        # initialize variables
        df = pd.DataFrame({'m/z': mz_list, 'intensity': int_list})
        centroids = []
        areas = []
        heights = []

        # integrate peak inside each window
        for window in self.windows:
            left, center, right = [x for x in window]
            mask = (df['m/z'] >= left) & (df['m/z'] <= right)
            df_window = df[mask].copy() # copy needed if doing the centroid calculation
            df_window['weighted intensity'] = df_window['m/z'] * df_window['intensity']
            """
            =========================================================
            should skip centroid calc if using instrument centroiding
            =========================================================
            """
            if COMPUTE_CENTROIDS:
                try:
                    centroid = df_window['weighted intensity'].sum() / df_window['intensity'].sum()
                except ZeroDivisionError:
                    centroid = center
                    print('...WARNING: No Intensities for scan', scan_num)
                centroids.append(centroid)
            else:
                centroids = self.zeroes
            if COMPUTE_AREAS:
                areas.append(trapz(df_window['intensity'], df_window['m/z']))
            else:
                areas = self.zeroes
            h = df_window['intensity'].max()
            h = h if h is not np.nan else 0.0 # replace NAN with zero
            heights.append(h)
        return centroids, areas, heights

    def set_windows(self):
        """Sets narrow integration windows around each reporter ion's predicted m/z position."""
        self.windows = np.array([[x-0.003, x, x+0.0025] for x in self.tmt_masses])
        return

    def process_msn_level(self, lc_name, msn_level=2, ion_count=15, min_intensity=100.0, reporter_ions = False):
        """Converts one Proteowizard TEXT formatted file to MSn format.
        Includes RT, etc. in MSn Information (I) lines. -PW June 2014
        Added summary statistics: missing charge states and m/z values -PW 6/15/2014
        Added extraction of reporter ions (if applicable) -PW 20180711
        """
        # get stuff for filenames
        if msn_level == 2:
            msn_extension = '.ms2'
        else:
            msn_extension = '.ms3'

        # get the header lines for Spectra setup
        header_block = []
        block = []
        in_spec = False     # this skips lines until first spectrum block
        spectra = None      # initial value
        for obj in self.log_obj:
            print('...Starting %s.txt.gz file scan' % (lc_name,), file=obj)
        for line in gzip.open(self.txt_name, 'rt'):
            line = line.strip()
            if line.startswith('chromatogramList'):
                break   # towards end of file after spectrum info
            if spectra is None:
                header_block.append(line)
            if line.startswith('spectrum:'):
                # create container for all spectra when at first "spectrum" line
                if spectra is None:
                    spectra = Spectra(ion_count, min_intensity, msn_level, header_block, self)
                # regular spectrum block processing
                in_spec = True
                if block:   # process previous spectrum block
                    spectrum = Spectrum(spectra, block)
                    spectra.add(spectrum)
                    if reporter_ions:
                        centroids, areas, heights = self.process_tmt_data(spectrum.mz_array, spectrum.int_array, spectrum.scan)
                        self.tmt_data.append(Reporter_ion(lc_name, spectrum.scan, spectrum.scan, centroids, areas, heights))
                    block = []  # reset block
            if in_spec:
                block.append(line)

        # need to parse last block
        if block:
            spectrum = Spectrum(spectra, block)
            spectra.add(spectrum)
            if reporter_ions:
                centroids, areas, heights = self.process_tmt_data(spectrum.mz_array, spectrum.int_array, spectrum.scan)
                self.tmt_data.append(Reporter_ion(lc_name, spectrum.scan, spectrum.scan, centroids, areas, heights))

        # write diagnostic stats from conversion
        for obj in self.log_obj:
            print('...Diagnostics for:', lc_name, file=obj)
        if spectra:
            spectra.report()

            # write data in desired formats
            total_scans = len(spectra.spectra)
            self.msn_total += total_scans
            for obj in self.log_obj:
                print('...writing MS%s file: %d scans passed cutoffs' % (msn_level, total_scans), file=obj)
            msn_name = os.path.join(self.raw_path, lc_name + msn_extension)
            spectra.write_msn(msn_name)
        spectra = None
        return

    def move_files(self):
        """Moves .ms2 and PAW_tmt.txt file into msn_files folder"""
        # get paths for source and destination folders
        source_folder = self.raw_path
        destination_folder = os.path.join(os.path.dirname(source_folder), 'msn_files')
        if not os.path.exists(destination_folder):
            os.mkdir(destination_folder)

        # get a list of all files to move
        files_to_move = []
        for pattern in ['*.ms2', '*.ms3', '*.PAW_tmt.txt']:
            files_to_move += glob.glob(pattern)

        # move the files (checks and deletes any files with the same names)
        for file in files_to_move:
            source = os.path.join(source_folder, file)
            dest = os.path.join(destination_folder, file)
            if os.path.exists(dest):
                os.remove(dest)
            shutil.move(source, dest)
        return

    def quit_gui(self):
        """Quits the GUI."""
        self.root.withdraw()
        self.root.update_idletasks()
        self.root.quit()
        return

    # methods to create sections of GUI
    def create_files_frame(self):
        """Lets the user select the RAW files."""
        dir_frame = ttk.Labelframe(self.root, text='Select files:')
        dir_frame.pack(fill=X, expand=YES, padx=5, pady=5)

        # variables
        self.raw_path_view = StringVar()
        self.def_location = r'C:\\'
        self.extension_list = [('RAW File(s)','*.RAW')]

        # defaults
        self.raw_path_view.set(r' Load RAW files (click the button on the left)')

        # creation
        Button(dir_frame, text=' Load files ', command=self.load_raw_files).pack(side=LEFT)
        Entry(dir_frame, textvariable=self.raw_path_view).pack(side=LEFT, fill=X, expand=YES)
        return

    def create_defaults_frame(self):
        """Lets the user change minimum on count, intensity, etc."""
        defaults_frame = ttk.Labelframe(self.root, text='MSConvert parameters:')
        defaults_frame.pack(fill=X, expand=YES, padx=5, pady=5)

        # variables
        self.ion_count = IntVar()
        self.min_intensity = DoubleVar()
        self.msn_level = IntVar()
        self.tmt_plex = IntVar()
        self.ms2_centroid = IntVar()
        self.ms3_centroid = IntVar()

        # creation
        self.create_entry(defaults_frame, 'Minimum ion count: ', self.ion_count).pack(fill=X, expand=YES)
        self.create_entry(defaults_frame, 'Minimum Intensity: ', self.min_intensity).pack(fill=X, expand=YES)
        self.create_radiobuttons(defaults_frame, 'Data to extract: ',
                                 [('MS2', 0), ('MS3', 1), ('MS2 TMT', 2), ('MS3 TMT', 3)],
                                 self.msn_level).pack(fill=X, expand=YES)
        self.create_radiobuttons(defaults_frame, 'TMT Reagents: ',
                                 [('None', 0), ('6-plex', 1), ('10-plex', 2),
                                  ('11-plex', 3), ('16-plex', 4), ('18-plex', 5)],
                                 self.tmt_plex).pack(fill=X, expand=YES)
        self.create_radiobuttons(defaults_frame, 'Centroid MS2 data: ',
                                 [('Yes', 0), ('No', 1)], self.ms2_centroid).pack(fill=X, expand=YES)
        self.create_radiobuttons(defaults_frame, 'Centroid MS3 data: ',
                                 [('Yes', 0), ('No', 1)], self.ms3_centroid).pack(fill=X, expand=YES)

        # set defaults
        self.ion_count.set(15)
        self.min_intensity.set(100)
        self.msn_level.set(0)
        self.tmt_plex.set(0)
        self.ms2_centroid.set(0)
        self.ms3_centroid.set(0)
        return

if __name__ == '__main__':
    convert = MSConvertGUI()
