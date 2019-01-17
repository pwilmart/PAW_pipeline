"""histo_GUI.py: Written by Billy Rathje, OHSU, 2014.
Also Phil Wilmarth, OHSU.
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
"""###############################
###############################
# converting to Python 3 -PW 9/16/2017


from tkinter import *
import tkinter.ttk as ttk
from tkinter import filedialog
from tkinter import messagebox

import os
import numpy as np

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from PAW_lib import FigureGenerator
from PAW_lib import DataInfoAndFilter
from PAW_lib import Threshold

import PAW_lib

import pickle    

class Histogram:
    ''' A histogram object tracks information needed to plot histogram data, including
    the textstring to display, the matplotlib, plot itself, the plot widget, x ranges, the placement of
    the vertical bar, the placement of curve fitting information, and a pointer to the original object
    storing numpy data (histogram) '''
    def __init__(self, h, notebook):
        self.ntt = h.ntt
        self.z = h.z
        self.histogram = h                                  # Pointer to data object
        self.plot = None                                    # The figure to plot
        self.xmin = -6                                      # position of x axis
        self.xmax = 12
        self.vertLine = None                                # green tracking bar
        self.vertLineThresholdSet = None                      # location of black dotted line
        self.textString = ''                                # The string to plot
        self.headerString = ''                              # first header line about text view
        self.FDRThreshold = 1.0                             # Flag for initial FDR threshold 'guess' value
        self.sparseCutoff = 50                              # Cutoff for number of Forward - Reverse spectra qualifying as spare
        self.histogram.sparseData = False                   # Does the current histogram contain spare data? Used in filtering step
        self.ln = self.initLn()                             # Placement of the line in the text view
        self.notebook = notebook                            # Pointer to parent notebook
        
        # Gui pointers
        self.can = None                                 # Pointer to canvas widget for this plot
        
        # Builds text string and plot
        self.makeTextString()
        self.makePlot()
        self.checkSparse()     # Remove lines from sparse figures
        
    def checkSparse(self):
        if(self.histogram.sparseData):
            self.vertLine.set_visible(False)
            self.vertLineThresholdSet.set_visible(False)
        
    def initLn(self):
        ''' Sets the initial value of the line in the table view to the first case of 1% FDR,
            or determines if the data is sparse and then sets the sparseData flag to True.
        '''
#        print('Z, ntt, mod:', self.z, self.ntt, self.histogram.mod)
        x = self.histogram.histo[self.histogram.histo.FDR.astype(float) <= self.FDRThreshold]
#        print('length of x:', len(x))
#        test = x[x['DiscScore'] > 2]
        test = self.histogram.histo[self.histogram.histo['DiscScore'] > 2]
#        print('length of test:', len(test))
        if len(test.index) > 0:
            # Sparse data - set discScore to 100.0
            if(test[test.index == test.index[0]].RRForward.item() - test[test.index == test.index[0]].RRReverse.item() <= self.sparseCutoff):
                self.histogram.sparseData = True
                # not sure about the next line
                self.zeroData = True
                print('Warning, sparse data frame...', ' ntt: ', self.ntt, ' z:', self.z, ' mod:', self.histogram.mod)
        if len(x.index) > 0:
            # Normal data, find first FDR <= self.FDRThreshold
            self.histogram.threshold = x.index[0]
            return x.index[0] + 1     # Off by one, I think because of removing header line
        # Frame is empty (no points where FDR is below cutoff)
        else:
            print('Warning, zero data frame...', ' ntt: ', self.ntt, ' z:', self.z, ' mod:', self.histogram.mod)
            self.histogram.zeroData = True
            self.histogram.sparseData = True
            return 0
            
    def setLn(self, val):
        '''
            Sets the position of the line. This is called by the gui widget's callbacks.
        '''
        if(self.histogram.sparseData):
            self.vertLine.set_visible(True)      # Make green line visible 
            
        self.ln = val
        #self.histogram.threshold = self.ln - 1   
        
    def updateVerticalLine(self):
        '''
            Updates the position of the green vertical line. Has to look for the discscore value
            corresponding to the current line in the text file. Ln in the text file corresponds to an
            index in the discscore column. This is called by the gui widget's callbacks.
        '''
        xval = self.histogram.histo['DiscScore'][self.ln]
        self.vertLine.set_xdata((xval, xval))
    
    def setThresh(self, event):
        '''
            This updates the black dotted line and sets the plot's score threshold to the value at the
            current location in the table/plot. This is called from one of the gui widget's callbacks.
        '''
        if(self.histogram.sparseData):        # Manually override sparse data cutoff
            self.histogram.sparseData = False
            self.vertLineThresholdSet.set_visible(True)   # Make black dotted line visible
            
        self.histogram.threshold = self.ln - 1
        xval = self.histogram.histo['DiscScore'][self.ln] 
        self.vertLineThresholdSet.set_xdata((xval, xval))
        self.can.draw()
        
    def makeTextString(self):
        '''
            Construct the text string
        '''
        
        # I use a cStringIO buffer in order to speed up the construction of the string. This buffer is passsed
        # into panda's toString method, which takes a 'buf' param
        from io import StringIO
        output = StringIO()

        # Make the text string, but don't display the index values.
        self.textString = self.histogram.histo.to_string(buf = output, index=False, col_space=12)
        self.textString = output.getvalue()
        output.close() # Close the buffer
        
        splitString = self.textString.split('\n')
        self.textString = '\n'.join(splitString[1:])   
        self.headerString = splitString[0]              # Make a seperate header line for the titles
        
    def makePlot(self):
        '''
            Construct a matplotlib plot for this histogram.
        '''
        self.plot = Figure(figsize=(5,3), dpi=100) 
        target = self.plot.add_subplot(111)
        decoy = self.plot.add_subplot(111)
        
        # center is used for the x axis, and is calculated from the bins computed with linspace.
        mybins = np.linspace(-8, 12, 301)
        center = (mybins[:-1] + mybins[1:]) / 2
        if(self.notebook.plotType.get() == 'Smoothed Plot'):
            target.fill(center, self.histogram.Smforward, color='b')
            decoy.fill(center, self.histogram.Smreverse, color='r')
        else:
            target.fill(center, self.histogram.forward, color='b')
            decoy.fill(center, self.histogram.reverse, color='r')
        
        # Adds the green vertical line. L is a tuple returned by the plot call that can be used to update the position later
        # with set_xdata.    
        greenVerticalLine = self.plot.add_subplot(111)
        l = greenVerticalLine.plot((0.0, 0.0), (target.axis()[2], target.axis()[3]), '-g', linewidth = 2.0)
        self.vertLine = l[0]
        
        # The black dotted line. Set the placement of it to the 1% FDR position.
        vertLineThresholdSet = self.plot.add_subplot(111)
        xval = self.histogram.histo['DiscScore'][self.ln]
        l = vertLineThresholdSet.plot((xval, xval), (target.axis()[2], target.axis()[3]), '--k', linewidth = 2.0)
        self.vertLineThresholdSet = l[0]

        # If there are xmin or max values supplied, set them.
        if self.xmin and self.xmax:
            target.set_xlim(self.xmin, self.xmax)
            decoy.set_xlim(self.xmin, self.xmax)
                    
        # set up labels
        target.set_title('Z=%s, NTT=%s, %s' % (str(self.z), str(self.ntt), self.histogram.mod))
        target.set_xlabel("Disc Score")
        target.set_ylabel("Counts")
        
        # makes sure the bottom axis is fully visible
        self.plot.tight_layout()
        
    def makePlotWidget(self, frcanvas):
        ''' Generate a tkinter canvas widget from the matplotlib plot. This is called
        from setup() in the BRNotebook class, and frcanvas is the gui frame into which the widget
        should be inserted. It should have already been constructed by the BRNotebook instance.
        '''
        canvas = BRCanvas(self.plot, master=frcanvas, gui=self.notebook.gui)
        canvas.show()
        canvas.get_tk_widget().pack(side=LEFT, fill=BOTH, expand=YES)
        canvas.get_tk_widget().bind('<Button-1>', canvas.focusCanvas)   # These callbacks need to be set here because
        canvas.get_tk_widget().bind('<Return>', self.setThresh)         # they are handled differently for deltamass plots.
        frcanvas.pack(side=TOP, expand=YES, fill=BOTH)                  # Other callbacks are set in the BRCanvas class that don't vary.
        canvas.histogram = self    # Retain reference to self           # Allow access to this class from the canvas widget.
        self.can = canvas                                               # Keep track of the canvas widget in an instance variable.
        
class DeltaMassHistogram(Histogram):
    ''' A DeltaMassHistogram object tracks information needed to plot deltamasshistogram data, including
    the textstring to display, the matplotlib, plot itself, the plot widget, x ranges, the placement of
    the vertical bar, the placement of curve fitting information, and a pointer to the original object
    storing numpy data (histogram).
    
    It is a subclass of histogram, so any functions not present here will default to histogram's functions. Most
    functions required custom implementations. Subclassing allows for function ovverriding, however, which makes for
    easier calls to things like setting the placement of the vertical selection bar from outside of the class (the function
    implementation will vary based on whether the histogram is a score or mass histogram, but outside the interface the call is the same)'''
    def __init__(self, h, notebook, FULL_RANGE = False):
        self.histogram = h        # Original histogram data
        self.notebook = notebook  # The notebook the plot resides in
        self.ln = self.initLn() # The current line in the text file
        self.z = h.z
        self.dm = h.dm
        self.plot = None        # The figure to plot
        self.target = None
        self.decoy = None  # Currently tracking target and decoy plots to get data for zoom
        self.xmin = None   # For plotting x axis
        self.xmax = None   # For plotting x axis
        self.vertLine = None    # The moving green vertical line
        self.textString = ''    # The string to plot
        self.headerString = ''  # The first title line above the large string
        self.FULL_RANGE = FULL_RANGE # Full range plot
        
        self.zoomFactor = 2     # Factor by which to zoom plots
        
        # 0 = low or left DM, 1 = high or right DM
        self.currentDM = 0     # Tracks the current threshold being modified (left or right)
        
        self.histogram.sparseData = False # for histogram/dmhistogram issues, remove later
        
        self.makeTextString()
        self.makePlot()
        
    
    def makeTextString(self):
        '''
        Builds the text string to display as data
        '''
        from io import StringIO
        output = StringIO()
        self.textString = self.histogram.histo.to_string(buf=output, index=False, col_space=12)
        self.textString = output.getvalue()
        output.close()
        
        splitString = self.textString.split('\n')
        self.textString = '\n'.join(splitString[1:])
        self.headerString = splitString[0]
        
    def makePlot(self):
        '''
        Builds the dm plots
        '''
        self.plot = Figure(figsize=(5,3), dpi=100)
        target = self.plot.add_subplot(111)
        decoy = self.plot.add_subplot(111)
        
        if self.dm == 0:
            target.set_xlim(-0.05, 0.05)
            decoy.set_xlim(-0.05, 0.05)
            self.xmin = -0.05
            self.xmax = 0.05
            Forward = self.histogram.forwardDeltaMassZero
            Reverse = self.histogram.reverseDeltaMassZero
            smForward = self.histogram.smForwardDeltaMassZero
            smReverse = self.histogram.smReverseDeltaMassZero
            mybins = np.linspace(-0.05, 0.05, 201)
        elif self.dm == 1:
            target.set_xlim(0.9, 1.1)
            decoy.set_xlim(0.9, 1.1)
            self.xmin = 0.9
            self.xmax = 1.1
            Forward = self.histogram.forwardDeltaMassOne
            Reverse = self.histogram.reverseDeltaMassOne
            smForward = self.histogram.smForwardDeltaMassOne
            smReverse = self.histogram.smReverseDeltaMassOne
            mybins = np.linspace(0.90, 1.10, 401)
        elif self.dm == 'ALL':
            self.xmin = -self.histogram.dmRange
            self.xmax = self.histogram.dmRange
            Forward = self.histogram.forwardDeltaMass
            Reverse = self.histogram.reverseDeltaMass
            smForward = self.histogram.smForwardDeltaMass
            smReverse = self.histogram.smReverseDeltaMass
            mybins = self.histogram.mybins
        elif self.dm == 2:
            self.xmin = -self.histogram.dmRange
            self.xmax = self.histogram.dmRange
            Forward = self.histogram.forwardDeltaMass
            Reverse = self.histogram.reverseDeltaMass
            smForward = self.histogram.smForwardDeltaMass
            smReverse = self.histogram.smReverseDeltaMass        
            mybins = self.histogram.mybins
                
        center = (mybins[:-1] + mybins[1:]) / 2
        #print('dm:', self.dm)
        #print('len y:', len(smForward))
        #print('len bins:', len(mybins))
        #print('len x:', len(center))
        if(self.notebook.plotType.get() == 'Smoothed Plot'):
            #target.fill(center, smForward, color='b')
            #decoy.fill(center, smReverse, color='r')
            target.plot(center, smForward, color='b')
            decoy.plot(center, smReverse, color='r')
        else:
            #target.fill(center, Forward, color='b')
            #decoy.fill(center, Reverse, color='r')
            target.plot(center, Forward, color='b')
            decoy.plot(center, Reverse, color='r')
        
        # Sets up the vertical line to display. Initializes it at the low threshold position.
        # Target.axis()[2] and [3] correspond to the current placements of the y-axis, which the line should mirror
        greenLine = self.plot.add_subplot(111)
        l = greenLine.plot((self.histogram.thresholdLow, self.histogram.thresholdLow), (target.axis()[2], target.axis()[3]), '-g', linewidth = 2.0)
        self.vertLine = l[0]
            
        # Sets up the low and high threshold black dotted lines
        # Ignores if no low/high thresholds present, like the out region that's not displayed.
        if(self.histogram.thresholdLow and self.histogram.thresholdHigh):
            low = self.plot.add_subplot(111)
            l= low.plot((self.histogram.thresholdLow, self.histogram.thresholdLow), (target.axis()[2], target.axis()[3]), '--k', linewidth = 2.0)
            high = self.plot.add_subplot(111)
            h = high.plot((self.histogram.thresholdHigh, self.histogram.thresholdHigh), (target.axis()[2], target.axis()[3]), '--k', linewidth = 2.0)
            self.low = l[0]
            self.high = h[0]
    
        # Set title and axis labels
        if(self.dm == 0 or self.dm == 1 or self.dm == 2):
            if self.dm != 2:
                target.set_title(str(self.dm) + ' Da Delta Mass')
            else:
                target.set_title('Full Range Delta Mass')
        else:
            target.set_title('Full range Delta Mass')
        target.set_xlabel("Deltamass (Da)")
        target.set_ylabel("Counts")
        
        target.set_xlim((self.xmin, self.xmax))
        
        self.plot.tight_layout()
        
        #track these for zooming
        self.target = target
        self.decoy = decoy
        
    def makePlotWidget(self, frcanvas):
        '''
        Build the actual tkinter widget and assign callbacks
        '''
        canvas = BRCanvas(self.plot, master=frcanvas, gui=self.notebook.gui)
        canvas.show()
        canvas.get_tk_widget().pack(side=LEFT, fill=BOTH, expand=YES)
        canvas.get_tk_widget().bind('<Button-1>', canvas.focusCanvas)
        canvas.histogram = self    # Retain reference to self
        self.can = canvas
        canvas.get_tk_widget().bind('<Double-Button-1>', self.zoom)
        canvas.get_tk_widget().bind('<Double-Button-3>', self.unZoom)
        canvas.get_tk_widget().bind('<Return>', self.setThresh)
        canvas.get_tk_widget().bind('<Left>', self.goToLeft)
        canvas.get_tk_widget().bind('<Right>', self.goToRight)
        
    
    def setThresh(self, event):
        '''
        Sets the threshold when the enter key is pressed depending on whether the left
        or right threshold is selected.
        '''
        if self.currentDM == 0:
            self.setLeftDM()
        else:
            self.setRightDM()
            
    def setLeftDM(self):
        '''
        Helper method for setThresh. Sets the threshold attribute for left (low) threshold. 
        '''
        self.histogram.thresholdLow = self.histogram.histo['deltaMass'][self.ln-1]
        self.low.set_xdata((self.histogram.thresholdLow, self.histogram.thresholdLow))
        self.can.draw()
    
    def setRightDM(self):
        '''
        Helper method for set Thresh. Sets the threshold attribute for right (high) threshold. 
        '''
        self.histogram.thresholdHigh = self.histogram.histo['deltaMass'][self.ln-1]
        self.high.set_xdata((self.histogram.thresholdHigh, self.histogram.thresholdHigh))
        self.can.draw()
    
    def goToLeft(self, event):
        '''
        Sets left threshold as active when left arrow key is pressed, and jumps green display line to
        that threshold.
        '''
        self.currentDM = 0
        self.vertLine.set_xdata((self.histogram.thresholdLow, self.histogram.thresholdLow))
        self.can.draw()
        xval = self.histogram.histo[abs(self.histogram.histo['deltaMass'] - self.histogram.thresholdLow) < .0005]
        self.ln = xval.index[0]
        
    def goToRight(self, event):
        '''
        Sets right threshold as active when right arrow key is pressed, and jumps green display line to
        that threshold.
        '''
        self.currentDM = 1
        self.vertLine.set_xdata((self.histogram.thresholdHigh, self.histogram.thresholdHigh))
        self.can.draw()
        xval = self.histogram.histo[abs(self.histogram.histo['deltaMass'] - self.histogram.thresholdHigh) < .0005]
        self.ln = xval.index[0]
        
    def zoom(self, event):
        ymin, ymax = self.target.get_ylim()
        self.target.set_ylim(0, ymax - ymax/(self.zoomFactor))
        self.can.draw()
        
    def unZoom(self, event):
        ymin, ymax = self.target.get_ylim()
        self.target.set_ylim(0, ymax + ymax*(self.zoomFactor))
        self.can.draw()
        
    def updateVerticalLine(self):
        #print('vertical line method')
        #print('len frame:', len(self.histogram.histo))
        #print('line:', self.ln)
        if self.ln > 1:
            xval = self.histogram.histo['deltaMass'][self.ln-1]
        else:
            xval = self.histogram.histo['deltaMass'][1]
        #print('xval:', xval)
        self.vertLine.set_xdata((xval, xval))
        
    def setLn(self, val):
        '''
        Setter for current line position in table
        '''
        self.ln = val
        
    def initLn(self):
        if(self.histogram.thresholdLow and self.histogram.thresholdHigh):
            # Setting the deltamass location in the table to the low threshold calculated will fail because the calculated
            # number is to a different level of float accuracy, so make sure they're within less than five thousandths of each other.
            xval = self.histogram.histo[abs(self.histogram.histo['deltaMass'] - self.histogram.thresholdLow) < .0005]
            if xval.empty:    # If there's no location for low threshold, just set to 1. Probably will happen for all
                return 1   # low mass data
            else:
                return xval.index[0]
        else:
            return 1
               
class BRNotebook(ttk.Notebook):
    ''' BRNotebook manages a ttk notebook object'''
    def __init__(self, gui=None, container = [], plotType = '', **args):
        ttk.Notebook.__init__(self, **args)
        self.histograms = {}      # Dictionary of histograms
        self.deltaMassHistograms = {}      # Dictionary of histograms
        self.fr = None         # main frame
        self.gui = gui          # reference to main GUI object
        self.container = container.container    # container of histogram data to process
        self.containerDeltaMass = container.dmContainer   # container of delatamass histograms
        self.containerStats = container.globalStats
        self.rawContainer = container                     # This is just a pointer back to the main containe - it has list and flag information
        #self.txtStats = container.txtStats
        self.plotType = plotType                          # Smoothed or not
        
        self.setup_deltaMass()    # First setup deltamass plots by default.
        #self.setup()            # initialization function to set up view...
        #self.setup_stats()
    
    def saveDMFigures(self):
        '''
        Save dm pdf figures
        '''
        sqt_container = os.path.dirname(os.getcwd())   # assumes we have set location to the folder with SQT files
        filter_container = os.path.join(sqt_container, 'filtered_files') # put the threshold figures in the filtered_files folder
        self.fig_container = os.path.join(filter_container, 'ThresholdFigures') # put the threshold figures in the filtered_files folder
        if not os.path.exists(filter_container):
            os.mkdir(filter_container)
        if not os.path.exists(self.fig_container):
            os.mkdir(self.fig_container)
        for figs in self.deltaMassHistograms.values():
            for fig in figs:    # using relative folder paths
                fig_file_name = os.path.join(self.fig_container, 'Mass_Figure_dm=' + str(fig.dm) + '_z=' + str(fig.z) + '.pdf')
                fig.plot.savefig(fig_file_name)
                
    def saveScoreFigures(self):
        '''
        Save score pdf figures (we have already set up folder for DM figures)
        '''
        for figs in self.histograms.values():
            for fig in figs:
                fig_file_name = os.path.join(self.fig_container, 'Score_Figure_dm=' + str(fig.histogram.dm) + '_z=' + str(fig.z) +
                                             '_ntt=' + str(fig.ntt) + '_mod=' + str(fig.histogram.mod) + '.pdf')
                fig.plot.savefig(fig_file_name)
        
            
    def setup(self):
        '''
        Setup score view
        '''
        
        # remove deltamass tabs from window
        for tab in self.tabs():
            self.forget(tab)
        
        # Create names for dm windows to reference when building tabs
        if self.rawContainer.accurateMass:
            daWindow = {0: '0 Da', 1: '1 Da', 2: 'out'}
        else:
            daWindow = {0: 'All'}
        
        # Loop through score data and set up the canvases
        for dm, dmFig in enumerate(self.container): 
            for f, fig in enumerate(dmFig):  
                theZ = self.rawContainer.zList[f]     # get z value for place in list
                # set up canvases
                self.fr = Frame(self.gui.root)        # self.fr is a large canvas for the whole tab
                frcanvas = Frame(self.fr)             # frcanvas is a frame for the entire mod notebook on top of the textview
               
                nb = ttk.Notebook(frcanvas)           # notebook to hold each mod
                
                frames = {}
                for mod in self.rawContainer.modList:
                    frames[mod] = Frame()               # Make a frame for each mod
                
                for fig in fig:     # Loops over ntt
                    for fig in fig:  # loops over mod
                        h = Histogram(fig, self)
                        h.makePlotWidget(frames[fig.mod])     # add plot to frame for specific mod
                        if f not in self.histograms:
                            self.histograms[f] = [h]
                        else:
                            self.histograms[f].append(h)

                for frcan in frames:
                    string = frcan
                    if frcan == ' ':      #
                        string = 'Unmod'  # I think this if is no longer needed - unmod is set in loading_TXT_files...
                    nb.add(frames[frcan], text = string)     # add each mod frame to notebook as seperate tab
                nb.pack(side=TOP, expand=YES, fill=BOTH)    
                frcanvas.pack(side=TOP, expand=YES, fill=BOTH)
                self.add(self.fr, text = (theZ, '+', '_', str(daWindow[dm])))    # Add the whole frame with textview and canvases to main notebook
                
                histPointer = h
                
                # Now set up the text view

                # Get just header line
                headerString = histPointer.headerString
                # set up header text view
                frtext = Frame(self.fr)              # Add text view as a frame in the frame holding frcanvas (the canvases notebook) and the text view
                text = Text(frtext, relief=SUNKEN)
                text.insert('1.0', headerString)
                text.pack(side=TOP, expand=NO, fill=X)
                text.config(width=1, height=1)
                text['wrap'] = 'none' 
                text['state'] = 'disabled'
                
                # Set up main text view
                textString = histPointer.textString
                text = BRText(frtext, relief=SUNKEN, notebook=self, gui=self.gui)
                text.insert('1.0', textString)
                
                text.focus()
                text.bind('<Button-1>', text.select)
                text.bind('<Up>', text.upKey)
                text.bind('<Down>', text.downKey)
                text.bind('<Return>', text.setThresh)
                text.pack(side=LEFT, expand=YES, fill=BOTH)
                text['state'] = 'disabled'  # prevents text editing
                text['wrap'] = 'none' 
                sbar = Scrollbar(frtext)
                sbar.config(command=text.yview)
                text.config(yscrollcommand=sbar.set)
                sbar.pack(side=RIGHT, fill=Y)
                hbar = Scrollbar(frtext, orient='horizontal')
                hbar.config(command=text.xview)
                text.config(xscrollcommand=hbar.set)
                hbar.pack(side=BOTTOM, fill=X)
            
                frtext.pack(side=BOTTOM, expand=YES, fill=BOTH)
                
                for histo in self.histograms[f]:
                    if histo.histogram.dm == daWindow[dm]:
                        histo.can.text = text              # Add pointer to current text view in canvas
                        text.canvas = histo.can            # Keep pointer to canvas in text view
                        text.see("%d.0" % text.canvas.histogram.ln)   # Go to current set line
                        text.refreshView()
            
    def setup_deltaMass(self):
        import pprint
        for f in range(len(self.containerDeltaMass[0])):
            self.fr = Frame(self.gui.root)
            plotsContainer = Frame(self.fr)
            plotsContainer.pack(side=TOP, expand=YES, fill=BOTH)
            frcanvas = Frame(plotsContainer)
            bottom_fr = Frame(plotsContainer)
            frcanvas.pack(side=TOP, expand=YES, fill=X)
            bottom_fr.pack(side=BOTTOM, expand=YES, fill=X)
           
            for x, fig in enumerate(self.containerDeltaMass):
                # set up canvases
                 
                # Since there's nothing in the container for full mass range, at the end of the list,
                # make an extra plot for full mass range ONLY if data has accurate mass. The low mass
                # container will only have 1 plot, the full range plot, so no need in that case.
                if self.gui.ACCURATE_MASS and x == (len(self.containerDeltaMass) - 1):
                    h = DeltaMassHistogram(fig[f], self, FULL_RANGE=True)
                    h.makePlotWidget(frcanvas)
                else:
                    h = DeltaMassHistogram(fig[f], self)
                    h.makePlotWidget(bottom_fr)
                    
                if f not in self.deltaMassHistograms:
                    self.deltaMassHistograms[f] = [h]
                else:
                    self.deltaMassHistograms[f].append(h)

            # Get just header line
            headerString = self.deltaMassHistograms[f][len(self.deltaMassHistograms[f])-1].headerString
            # set up related text view
            frtext = Frame(self.fr)
            text = Text(frtext, relief=SUNKEN)
            text.insert('1.0', headerString)
            text.pack(side=TOP, expand=NO, fill=X)
            text.config(width=1, height=1)
            text['wrap'] = 'none' 
            text['state'] = 'disabled'
            
            textString = self.deltaMassHistograms[f][len(self.deltaMassHistograms[f])-1].textString
            text = BRText(frtext, relief=SUNKEN, notebook=self, gui=self.gui)
            text.insert('1.0', textString)
            
            #text.tag_add(SEL, '1.0', '1.200')
            text.focus()
            text.bind('<Button-1>', text.select)
            text.bind('<Up>', text.upKey)
            text.bind('<Down>', text.downKey)
            text.bind('<Return>', text.setThresh)
            text.bind('<Left>', text.goToLeft)
            text.bind('<Right>', text.goToRight)
            text.pack(side=LEFT, expand=YES, fill=BOTH)
            text['wrap'] = 'none' 
            text['state'] = 'disabled'  # prevents text editing
            sbar = Scrollbar(frtext)
            sbar.config(command=text.yview)
            text.config(yscrollcommand=sbar.set)
            sbar.pack(side=RIGHT, fill=Y)
            hbar = Scrollbar(frtext, orient='horizontal')
            hbar.config(command=text.xview)
            text.config(xscrollcommand=hbar.set)
            hbar.pack(side=BOTTOM, fill=X)
        
            frtext.pack(side=BOTTOM, expand=YES, fill=BOTH)
            
            for histo in self.deltaMassHistograms[f]:
                histo.can.text = text
                text.canvas = histo.can
                
                #if not text.canvas.histogram.ln:
                #    continue
                    
                text.see("%d.0" % text.canvas.histogram.ln)
                text.refreshView()
            
            theZ = self.rawContainer.zList[f]     
            self.add(self.fr, text = (theZ, '+_DM'))
            f += 1
     
    def setup_stats(self):
        self.fr = Frame(self.gui.root)
        statsContainer = Frame(self.fr)
        Label(statsContainer, text = '1+ Target\t').grid(column = 1, row = 0)
        Label(statsContainer, text = '1+ Decoy\t').grid(column = 2, row = 0)
        Label(statsContainer, text = '2+\t').grid(column = 3, row = 0)
        Label(statsContainer, text ='2+\t').grid(column = 4, row = 0)
        Label(statsContainer, text ='3+\t').grid(column = 5, row = 0)
        Label(statsContainer, text ='3+\t').grid(column = 6, row = 0)
        Label(statsContainer, text ='4+\t').grid(column = 7, row = 0)
        Label(statsContainer, text ='4+\t').grid(column = 8, row = 0)
        Label(statsContainer, text ='Unmod').grid(column = 0, row = 1)
        
        Label(statsContainer, text ='Full').grid(column = 0, row = 2)
        
        DM = len(self.containerDeltaMass) - 1
        NTT = len(self.container[0][0])
        if(NTT == 3):
            for z in range(len(self.containerDeltaMass[0])):
                Label(statsContainer, text = str(self.containerStats.target_subclass[DM][z][0][0])).grid(column = z+1, row = 2)
                Label(statsContainer, text = str(self.containerStats.decoy_subclass[DM][z][0][0])).grid(column = z+2, row = 2)  
            Label(statsContainer, text ='Semi').grid(column = 0, row = 3)
            for z in range(len(self.containerDeltaMass[0])):
                Label(statsContainer, text = str(self.containerStats.target_subclass[DM][z][1][0])).grid(column = z+1, row = 3)
                Label(statsContainer, text = str(self.containerStats.decoy_subclass[DM][z][1][0])).grid(column = z+2, row = 3)  
            Label(statsContainer, text ='Non').grid(column = 0, row = 4)
            for z in range(len(self.containerDeltaMass[0])):
                Label(statsContainer, text = str(self.containerStats.target_subclass[DM][z][2][0])).grid(column = z+1, row = 4)
                Label(statsContainer, text = str(self.containerStats.decoy_subclass[DM][z][2][0])).grid(column = z+2, row = 4)
        if(NTT < 3):
            for z in range(len(self.containerDeltaMass[0])):
                Label(statsContainer, text = str(self.containerStats.target_subclass[DM][z][0][0])).grid(column = z+1, row = 2)
                Label(statsContainer, text = str(self.containerStats.decoy_subclass[DM][z][0][0])).grid(column = z+2, row = 2)  
            Label(statsContainer, text ='Semi').grid(column = 0, row = 3)
            for z in range(len(self.containerDeltaMass[0])):
                Label(statsContainer, text = str(self.containerStats.target_subclass[DM][z][1][0])).grid(column = z+1, row = 3)
                Label(statsContainer, text = str(self.containerStats.decoy_subclass[DM][z][1][0])).grid(column = z+2, row = 3)  
            Label(statsContainer, text ='Non').grid(column = 0, row = 4)
            for z in range(len(self.containerDeltaMass[0])):
                Label(statsContainer, text = "----------").grid(column = z+1, row = 4)
                Label(statsContainer, text = "----------").grid(column = z+2, row = 4)    
        
        Label(statsContainer, text ='Totals').grid(column = 0, row = 5)
        Label(statsContainer, text = self.containerStats.target_filtered).grid(column = 1, row = 5)
        Label(statsContainer, text = self.containerStats.decoy_filtered).grid(column = 3, row = 5)
        statsContainer.pack()
        self.add(self.fr, text = 'Stats')   
    

class BRCanvas(FigureCanvasTkAgg):
    ''' Manages a FigureCanvasTkAgg object '''
    def __init__(self, parent=None, gui=None, **args):
        FigureCanvasTkAgg.__init__(self, parent, **args)
        self.text = None    # Associated text view
        self.ntt = 0        # Ntt for canvas
        self.charge = 0     # Charge for canvas
        self.gui = gui      # Reference to main gui object
        self.histogram = None
        
    def focusCanvas(self, event):
        ''' Sets focus to current canvas and updates text view '''
        self.text.canvas = self
        self.text['state'] = 'normal'
        self.text.focus()
        self.text.delete('0.0', END)
        self.text.insert('1.0', self.histogram.textString)
        self.text.see("%d.0" % self.text.canvas.histogram.ln)
        #self.text.tag_remove(SEL, '0.0', END)
        
        self.text.tag_remove('highlight', '0.0', END)
        self.text.tag_add('highlight', "%d.0" % self.text.canvas.histogram.ln, "%d.0" % (self.text.canvas.histogram.ln + 1))
        self.text.tag_configure('highlight', background = 'sky blue')
        #self.text.tag_add(SEL, "%d.0" % self.text.canvas.histogram.ln, "%d.0" % (self.text.canvas.histogram.ln + 1))
        self.text['state'] = 'disabled'
        self.get_tk_widget().focus_set()
        return 'break'
    
class BRText(Text):
    ''' Manages a text view object '''
    def __init__(self, parent=None, notebook=None, gui=None, **args):
        Text.__init__(self, parent, takefocus=0, **args)
        self.canvas = None            # Associated canvas
        self.notebook = notebook      # Reference to notebook textview is in
        self.gui = gui                # Reference to main GUI object
        
    def refreshView(self):
        ''' Helper method, updates figure for text view and redraws canvas '''
        #self.tag_remove(SEL, '0.0', END)
        #self.tag_add(SEL, "%d.0" % self.canvas.histogram.ln, "%d.0" % (self.canvas.histogram.ln + 1))
        self.tag_remove('highlight', '0.0', END)
        self.tag_add('highlight', "%d.0" % self.canvas.histogram.ln, "%d.0" % (self.canvas.histogram.ln + 1))
        self.tag_configure('highlight', background = 'sky blue')
        self.focus()
        
        # Updates differently depending on whether class is deltamass or score histogram
        self.canvas.histogram.updateVerticalLine()
        self.canvas.draw()

    def select(self, event):
        ''' Callback for selection of textview line with mouse'''
        self.canvas.histogram.setLn(int(self.index(CURRENT).split('.')[0]))  # Gets just the line number
        self.refreshView()
        return 'break'      # 'break' overrides widget's default behavior
    
    def upKey(self, event):
        ''' Callback for selection of textview line with upKey'''
        self.canvas.histogram.setLn(int(self.canvas.histogram.ln) - 1)
        self.see("%d.0" % self.canvas.histogram.ln)
        self.refreshView()
        return 'break'      # 'break' overrides widget's default behavior
    
    def downKey(self, event):
        ''' Callback for selection of textview line with downKey'''
        self.canvas.histogram.setLn(int(self.canvas.histogram.ln) + 1)
        self.see("%d.0" % self.canvas.histogram.ln)
        self.refreshView()
        return 'break'      # 'break' overrides widget's default behavior

    def setThresh(self, event):
        self.canvas.histogram.setThresh(event)
        return 'break'
    
    def goToRight(self, event):
        self.canvas.histogram.goToRight(event)
        return 'break'
    
    def goToLeft(self, event):
        self.canvas.histogram.goToLeft(event)
        return 'break'
    
class GUI:
    ''' Manages the main GUI window. 
    Also starts parsing in files, loading pandas structures'''
    
    def __init__(self, folder=None):

        self.root = Tk()
        self.root.title('PAW Histogram GUI')
        if not folder:
            folder = os.getcwd()
        self.folder = folder

        # this is the starter window         
        self.modal = Toplevel(self.root)
        self.modal.geometry("%dx%d%+d%+d" % (300, 200, 250, 125))
        self.modal.title('PAW Set Up Dialog')
        ttk.Button(self.modal, text="Select Top Hit Summary Files", command=self.select_files).pack(pady=5)
        variable = StringVar(self.modal)
        variable.set("Plot") # default value
        ttk.OptionMenu(self.modal, variable, "Standard Plots", "Smoothed Plots").pack(pady=5)
        self.massAccuracy = StringVar(self.modal)
        self.massAccuracy.set("High") # default value (gets over-written during file loading)
        ttk.OptionMenu(self.modal, self.massAccuracy, "High Resolution", "Low Resolution").pack(pady=5)
        ttk.Button(self.modal, text="Load and Plot Histograms", command=self.exit_modal).pack(pady=5)
                
        self.root.protocol('WM_DELETE_WINDOW', self.onExit) # handle exit
        self.modal.protocol('WM_DELETE_WINDOW', self.exit_modal) # cannot get setup window to delete on mac
        self.root.withdraw()
        self.modal.attributes('-topmost', 1)
#        self.modal.attributes('-topmost', 0)   
        self.root.wait_window(self.modal)   # this waits for the user to set the files, resolution, etc.
        self.root.deiconify()

        # when we get here, we are starting the histogramming        
        # Setup flags
        self.sparseDiscScore = 100.0        
        self.ACCURATE_MASS = True
        self.SMOOTHED = True
        if self.massAccuracy.get() == 'Low':
            self.ACCURATE_MASS = False
        if variable.get() == 'Plot':
            self.SMOOTHED = False
            
        # Make histograms
        self.container = FigureGenerator(files=self.txt_files, accurateMass=self.ACCURATE_MASS, smoothed=self.SMOOTHED)
        
        print('Generating GUI...')
        
        # Main frames
        self.buttonFrame = Frame(self.root)
        self.buttonFrame.pack(side=TOP, fill=X)
        self.computeScoreHistograms = ttk.Button(self.buttonFrame, text = 'Compute Score Histograms', takefocus=0, command=self.compute_score_histograms)
        self.computeScoreHistograms.pack(side=LEFT, padx=5, pady=2)
        ttk.Button(self.buttonFrame, text = 'Show mass windows', takefocus=0, command=self.get_masses).pack(side=LEFT, padx=5, pady=2)
        
        self.notebook = BRNotebook(gui=self, container=self.container, plotType = variable, takefocus=0)
        self.notebook.pack(side=BOTTOM)
        
        #f = Frame(self.root)
        #Label(f, text = 'Sidebar: Additional widgets here...').pack()
        #f.pack(side=BOTTOM)
        
        mainloop()

    def onExit(self):
        import sys
        """Properly closes down the GUI program."""
        self.root.withdraw()
        self.root.update_idletasks()
        self.root.destroy()
        sys.exit()
        
    def pickleDeltaMass(self):
        with open("output.pkl", "wb") as fout:
            pickle.dump(self.container, fout)
        
    def pickleScore(self):
        with open("output_scores.pkl", "wb") as fout:
            pickle.dump(self.container, fout)
    
    def exit_modal(self):
        self.modal.withdraw()
        self.modal.update_idletasks()
        self.modal.destroy()
    
    def select_files(self):
        self.txt_files = PAW_lib.get_files(self.folder, [('Text files', '*.txt'), ('PAW Text files', '*.PAW.txt')],
                                           'Select the Top-hit TXT files')   # returns full paths
        if not self.txt_files: sys.exit() # cancel button response
        self.folder = os.path.dirname(self.txt_files[0])
        os.chdir(self.folder)
    
    def get_scores(self):
        s = ''
        for score in self.container.container:
            for score in score:
                for score in score:
                    for score in score:
                        s += (str(score.dm) + " , " 
                            + str(score.z) + " +, " 
                            + str(score.ntt) + " tt, " +
                            str(score.mod) + " mod: ")
                        if(score.sparseData):
                            s += str(100) + '\n'
                        else:
                            s += ('%0.4f' % score.histo.DiscScore[score.threshold]) + '\n'
        messagebox.showinfo(title = "Scores", message = s)
        
        #self.get_stats()
    
    def get_masses(self):
        s = ''
        for dm in range(len(self.container.dmList)):
            for z in range(len(self.container.zList)):
                if dm == 2:
                    s += '\n' + ('%d' % self.container.dmContainer[dm][z].z) + "+ , " + " Outside"
                    s += "\n\tLow: " + ('%0.4f' % -self.container.dmContainer[dm][z].dmRange) + " \n\tHigh: " + ('%0.4f' % self.container.dmContainer[dm][z].dmRange)
                elif dm == 0 or dm == 1:
                    s += '\n' + ('%d' % self.container.dmContainer[dm][z].z) + "+ , " + ('%d' % self.container.dmContainer[dm][z].dm) + " Da"
                    s += "\n\tLow: " + ('%0.4f' % self.container.dmContainer[dm][z].thresholdLow) + " \n\tHigh: " + ('%0.4f' % self.container.dmContainer[dm][z].thresholdHigh)
                      
        messagebox.showinfo(title = "Mass Thresholds", message = s)
    
    def compute_score_histograms(self):
        self.container.regenerateScorePlots()
        self.notebook.setup()
        ttk.Button(self.buttonFrame, text = 'Show score thresholds', takefocus=0, command=self.get_scores).pack(side=LEFT, padx=5, pady=2)
        ttk.Button(self.buttonFrame, text = 'Filter Data', takefocus=0, command=self.exportToFilterer).pack(side=LEFT, padx=5, pady=2)
        self.computeScoreHistograms['state'] = 'disabled'
        
        # Save figures
        self.notebook.saveDMFigures()
        
    def get_stats(self):
        self.container.get_stats_helper()
        
    def exportToFilterer(self):
        import time
        
        filterer = DataInfoAndFilter(self.folder, self.container.f.getFrame(), self.container.txtObjects,
            self.container.dmList, self.container.zList, self.container.nttList, self.container.specialCharsList,
            self.container.minLength, self.container.maxMods, self.container.peptideMassTol)
       
        filterer.get_pre_stats()
        
        masses = [[Threshold() for dm in self.container.dmList] for z in self.container.zList]
        for dm in range(len(self.container.dmList)):
            for z in range(len(self.container.zList)):
                if(dm == 2):
                    masses[z][dm].low = -1 * self.container.dmContainer[dm][z].dmRange
                    masses[z][dm].high = self.container.dmContainer[dm][z].dmRange
                else:
                    masses[z][dm].low = self.container.dmContainer[dm][z].thresholdLow
                    masses[z][dm].high = self.container.dmContainer[dm][z].thresholdHigh
        
        
        scores = [[[[100.0 for mod in self.container.specialCharsList] for ntt in self.container.nttList] for z in self.container.zList] for dm in self.container.dmList]
        for dm in range(len(self.container.dmList)):
            for z in range(len(self.container.zList)):
                for ntt in range(len(self.container.nttList)):
                    for mod in range(len(self.container.specialCharsList)):
                        ref = self.container.container[dm][z][ntt][mod]
                        if ref.sparseData:
                            scores[dm][z][ntt][mod] = self.sparseDiscScore
                        else:
                            scores[dm][z][ntt][mod] = ref.histo.DiscScore[ref.threshold]
    
        filter_frame = filterer.filter_with_stats(masses, scores)   # probably do not need to have a returned dataframe
        filterer.copy_params_files()
        
        for obj in filterer.write:
            print('\n...finished.', time.asctime(), file=obj)
        filterer.log_file.close()
        
        # Save figures
        self.notebook.saveScoreFigures()
         
    def _quit(self):
        root.quit()     # stops mainloop
        root.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate
#########################
# default folder location
folder = os.getcwd()   # this is a safe default
folder = "E:"   # or set to something useful for your system
#########################
gui = GUI(folder)
