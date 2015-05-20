#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     rmPipeViewer.py                                                   #
#                                                                             #
# PURPOSE:  A graphical interface designed to view the results of the Level 5 #
#           RM-pipeline prototype.                                            #
#                                                                             #
# MODIFIED: 19-May-2015 by cpurcell                                           #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
# App               ... class to create the main window and sub-windows       #
# SessChooseFrame   ... frame to browse to and load a session / results       #
# TabPanelFrame     ... frame holding the tabbed notebook widget              #
# PipeInputsFrame   ... frame to display the pipeline inputs                  #
# ResultsManager    ... class to hold the database tables in memory           #
# SpectraPlotFrame  ... frame defining the spectra plots                      #
# RMsynthPlotFrame  ... frame defining the RM-synthesis plots                 #
#                                                                             #
#=============================================================================#

# Default session
defaultSessionDir = "testSession"

# Window geometry
geometryBrowseWin = "1024x700"
geometryPlotWin = "900x800"

#-----------------------------------------------------------------------------#

import os
import sqlite3
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import Tkinter as tk
import ttk
import tkFileDialog
import tkMessageBox
import tkFont

from Imports.util_PPC import PipelineInputs
from Imports.util_DB import *
from Imports.util_tk import *
from Imports.util_plotTk import *

# Turn off print statements buffering
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)


#-----------------------------------------------------------------------------#
class App:
    """Class defining the RM Pipeline Viewer."""

    def __init__(self, root):
        self.root = root
        self.root.title("RM Pipeline Viewer - Browser Window")
        self.root.geometry(geometryBrowseWin)
        self.root.resizable(True, True)
        self.root.protocol("WM_DELETE_WINDOW", self.applicationExit)
        
        # Create the session chooser panel
        self.sessFrm = SessChooseFrame(self.root)
        self.sessFrm.grid(row=0, column=0, padx=5, pady=5, sticky="NSEW")
        self.sessFrm.sesDir.set(defaultSessionDir)

        # Create the tabbed "notebook" to hold the main pages
        self.tabPanel = TabPanelFrame(self.root)
        self.tabPanel.grid(row=1, column=0, padx=5, pady=5, sticky="NSEW")
        
        # Allow the action frame to expand
        self.root.rowconfigure(1, weight=1)
        self.root.columnconfigure(0, weight=1)
        
        # Bind virtual events generated by sub-widgets
        self.root.bind("<<load_session>>", self.on_load_session)
        self.root.bind("<<cat_row_selected>>", self.on_select_catalogue)
        self.root.bind("<<specparm_row_selected>>", self.on_select_specparm)
        self.root.bind("<<rmsynth_row_selected>>", self.on_select_rmsynth)

        # Create the visualisation window and set the focus back to root
        self.visWin = tk.Toplevel(self.root)
        self.visWin.title("RM Pipeline Viewer - Plotting Window")
        self.visWin.geometry(geometryPlotWin)
        self.visWin.resizable(True, True)
        self.visWin.protocol("WM_DELETE_WINDOW", self.applicationExit)        
        self.visWin.columnconfigure(0, weight=1)
        self.visWin.rowconfigure(0, weight=1)

        # Create a master frame for the visualisation window
        self.visWinFrame = tk.Frame(self.visWin)
        self.visWinFrame.grid(row=0, column=0, padx=5, pady=5, sticky="NSEW")
        self.visWinFrame.columnconfigure(0, weight=1)
        self.visWinFrame.rowconfigure(0, weight=1)
        self.tmpVisLab = tk.Label(self.visWinFrame, font=("Helvatica", 30),
                                  text="Plotting Window")
        self.tmpVisLab.grid(row=0, column=0, padx=5, pady=5, sticky="NSEW")

        # Force a minimum size on the windows
        self.root.update()
        self.root.minsize(self.root.winfo_width(), 
                          self.root.winfo_height())
        self.visWin.update()
        self.visWin.minsize(self.visWin.winfo_width(),
                            self.visWin.winfo_height())
        
        # Raise the main window and start the event loop
        self.root.after(1, lambda: self.root.focus_force())
        self.root.after(1, lambda: self.root.lift())
        self.root.mainloop()

    def reset_plotting_window(self):
        """Reset the plotting window."""
        self.visWinFrame.destroy()
        self.visWinFrame = tk.Frame(self.visWin)        

    def applicationExit(self):
        """Exit the application cleanly."""
        self.root.destroy()

    # Event handlers ---------------------------------------------------------#
    
    def on_load_session(self, event=None):
        """Load in the chosen pipeline session into the GUI.
        Reset the GUI, read the pipeline inputs from the ASCII input file
        and calculate basic parameters. Populate the inputs tab and load 
        the tables from the database."""
        
        # Check that the session exists and has an input file
        sessionPath = self.sessFrm.sesDir.get()
        if not os.path.exists(sessionPath):
            errStr = "The directory does not exist:\n\n%s" % sessionPath
            tkMessageBox.showinfo("Error", errStr)
            return
        inputFile = sessionPath + "/inputs.config"
        if not os.path.exists(inputFile):
            errStr = "Missing input file.\n\nThe file 'inputs.config' " + \
                "is missing from the session directory."
            tkMessageBox.showinfo("Error", errStr)
            return
        
        # Create a PipelineInput class, in the process reading in the 
        # pipeline input file and calculating the derived parameters
        pipeInpObj = PipelineInputs(inputFile)
        pipeInpObj.calculate_derived_parms(resetPhiSamp=False)        
        pDict = pipeInpObj.get_flat_dict(includeDerived=True)
        
        # Set the pipeline input parameters in the GUI
        self.tabPanel.inputsNBTab.clear_entries()
        self.tabPanel.inputsNBTab.set_entries(pDict)

        # Create a ResultsManager class, in the process loading the tables
        # from the database into numpy record arrays
        print "Loading results into memory ...",
        self.resultsMan = ResultsManager(sessionPath)
        print "done."
        
        # Load the tables into the GUI
        self.tabPanel.catNBTab.clear_entries()
        print "Loading catalogue into GUI ...",
        self.tabPanel.catNBTab.insert_table(self.resultsMan.catRec)
        print "done."
        self.tabPanel.specParmNBTab.clear_entries()
        print "Loading spectral table into GUI ...",
        self.tabPanel.specParmNBTab.insert_table(self.resultsMan.specParmRec)
        print "done."
        self.tabPanel.rmsynthNBTab.clear_entries()
        print "Loading RM-synthesis table into GUI ...",
        self.tabPanel.rmsynthNBTab.insert_table(self.resultsMan.dirtyFDFRec)
        print "done."

        # Reset the plotting window
        self.reset_plotting_window()
        
    def on_select_catalogue(self, event=None):
        """Plot the source details when a catalogue entry is selected."""
        pass
        
    def on_select_specparm(self, event=None):
        """Plot the spectra when a catalogue entry is selected."""
        
        self.reset_plotting_window()

        # Create a new SpectraPlotFrame panel
        indx = event.widget.get_indx_selected()
        self.specPlotFrame = SpectraPlotFrame(self.visWinFrame, 
                                              self.resultsMan, indx)
        self.specPlotFrame.grid(row=0, column=0, padx=5, pady=5,
                                sticky="NSEW")

        # Grid the container frame last
        self.visWinFrame.columnconfigure(0, weight=1)
        self.visWinFrame.rowconfigure(0, weight=1)
        self.visWinFrame.grid(row=0, column=0, padx=5, pady=5, sticky="NSEW")
        
    def on_select_rmsynth(self, event=None):
        """Plot the dirty FDF when a catalogue entry is selected."""
        
        self.reset_plotting_window()
        
        # Create a new RMsynthPlotFrame panel
        indx = event.widget.get_indx_selected()
        self.rmsynthPlotFrame = RMsynthPlotFrame(self.visWinFrame,
                                                 self.resultsMan, indx)
        self.rmsynthPlotFrame.grid(row=0, column=0, padx=5, pady=5,
                                   sticky="NSEW")

        # Grid the container frame last
        self.visWinFrame.columnconfigure(0, weight=1)
        self.visWinFrame.rowconfigure(0, weight=1)
        self.visWinFrame.grid(row=0, column=0, padx=5, pady=5, sticky="NSEW")


#-----------------------------------------------------------------------------#
class SessChooseFrame(tk.Frame):
    """Frame containing a session chooser widget. Allows the user to enter or
    browse to a session directory and load the parameters and results into 
    the GUI."""
    
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        
        # Create the session chooser frame
        self.frmChoose = tk.LabelFrame(self, text=" Current Session ")
        self.frmChoose.grid(column=0, row=0, padx=5, pady=2, sticky="NSEW")
        self.labIns = tk.Label(self.frmChoose, justify="left", anchor="e",
            text="Enter or browse to a session directory and click 'Load':")
        self.labIns.grid(row=0, column=0, padx=5, pady=5, sticky="NW")
        self.sesDir = tk.StringVar()
        self.entSes = ttk.Entry(self.frmChoose, textvariable=self.sesDir,
                                width=50)
        self.entSes.grid(row=1, column=0, padx=5, pady=5, sticky="W")
        self.btnBrowse = ttk.Button(self.frmChoose, text="Browse",
                                    command=self._handlerBrowseButton)
        self.btnBrowse.grid(row=0, column=1, padx=5, pady=2,sticky="NW" )
        self.btnLoad = ttk.Button(self.frmChoose, text="Load",
                                  command=self._handlerLoadButton)
        self.btnLoad.grid(row=1, column=1, padx=5, pady=5, sticky="NW" )
        
        # Create the session status frame
        self.frmSess = tk.LabelFrame(self, text=" Session Status ")
        self.frmSess.grid(column=1, row=0, padx=5, pady=2, sticky="NSEW")
        #self.labName = tk.Label(self.frmSess, justify="left", anchor="e",
        #                        text="Session Name:")
        #self.labName.grid(row=0, column=0, padx=5, pady=5, sticky="NW")
        #self.entSes = ttk.Entry(self.frmSess, width=20, state="disabled")
        #self.entSes.grid(row=0, column=1, padx=5, pady=5, sticky="NW")
        #self.labDate = tk.Label(self.frmSess, justify="left", anchor="e",
        #                        text="Creation Date:")
        #self.labDate.grid(row=1, column=0, padx=5, pady=5, sticky="NW")
        #self.entDate = ttk.Entry(self.frmSess, width=20, state="disabled")
        #self.entDate.grid(row=1, column=1, padx=5, pady=5, sticky="NW")
        
        # Set the expansion behaviour
        #self.frmChoose.columnconfigure(0, weight=1)
        #self.frmChoose.rowconfigure(0, weight=1)
        #self.frmChoose.rowconfigure(1, weight=1)
        self.columnconfigure(1, weight=1)

    def _handlerBrowseButton(self):
        """Open the file selection dialog."""
        sesDir = tkFileDialog.askdirectory(parent=self, initialdir=".",
                                    title="Please select a session directory")
        if sesDir!="":        
            self.sesDir.set(sesDir)
            
    def _handlerLoadButton(self):      
        """Raise a <<load_session>> event in the parent window."""
        self.event_generate("<<load_session>>")


#-----------------------------------------------------------------------------#
class TabPanelFrame(tk.Frame):
    """Frame containing the a tabbed 'notebook' used to layout the inputs and
    result tables of the RM pipeline in a logical order."""

    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        
        # Create the notebook and grid
        self.nb = ttk.Notebook(self, name="notebook")
        self.nb.enable_traversal()
        self.nb.grid(row=0, column=0, padx=5, pady=5, sticky="NSEW")
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        # Add the panels in order
        self._add_pipeinputs_tab()
        self._add_catalogue_tab()
        self._add_specparm_tab()
        self._add_rmsynth_tab()
        
    def _add_pipeinputs_tab(self):
        """Add the 'Pipeline Inputs' table to the notebook."""
        tabFrame = tk.Frame(self.nb)
        tabFrame.rowconfigure(0, weight=1)
        tabFrame.columnconfigure(0, weight=1)
        self.inputsNBTab = PipeInputsFrame(tabFrame)
        self.inputsNBTab.grid(column=0, row=0, sticky="NSEW")
        self.nb.add(tabFrame, text="Pipeline Inputs", padding=5)
        
    def _add_catalogue_tab(self):
        """Add the 'Source Catalogue' table to the notebook."""
        tabFrame = tk.Frame(self.nb)
        tabFrame.rowconfigure(0, weight=1)
        tabFrame.columnconfigure(0, weight=1)
        title = "Input Catalogue:"
        foot = "Click on a column header to sort up or down."
        foot += "\nClick on a row to plot the data for that entry."
        self.catNBTab = TableFrame(tabFrame, title=title,  foot=foot,
                                   virtEvent="<<cat_row_selected>>" )
        self.catNBTab.grid(column=0, row=0, sticky="NSEW")
        self.nb.add(tabFrame, text="Source Catalogue", padding=5)
        
    def _add_specparm_tab(self):
        """Add the 'Spectral Parameters' table to the notebook."""
        tabFrame = tk.Frame(self.nb) 
        tabFrame.rowconfigure(0, weight=1)
        tabFrame.columnconfigure(0, weight=1)
        title = "Parameters of the extracted spectra:"
        foot = "Click on a column header to sort up or down."
        foot += "\nClick on a row to plot the spectra for that entry."
        self.specParmNBTab = TableFrame(tabFrame, title=title,  foot=foot,
                                        virtEvent="<<specparm_row_selected>>" )
        self.specParmNBTab.grid(column=0, row=0, sticky="NSEW")
        self.nb.add(tabFrame, text="Spectral Parameters", padding=5)
        
    def _add_rmsynth_tab(self):
        """Add the 'RM-Synthesis Results' table to the notebook."""
        tabFrame = tk.Frame(self.nb) 
        tabFrame.rowconfigure(0, weight=1)
        tabFrame.columnconfigure(0, weight=1)
        title = "RM-synthesis results:"
        foot = "Click on a column header to sort up or down."
        foot += "\nClick on a row to plot the FDF for that entry."
        self.rmsynthNBTab = TableFrame(tabFrame, title=title, foot=foot,
                                       virtEvent="<<rmsynth_row_selected>>" )
        self.rmsynthNBTab.grid(column=0, row=0, sticky="NSEW")
        self.nb.add(tabFrame, text="RM-Synthesis Results", padding=5), 


#-----------------------------------------------------------------------------#
class PipeInputsFrame(tk.Frame):
    """Frame presenting the pipeline inputs to the user."""
    
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        
        # Layout the Dataset & Extraction column
        self.labTitle1 = tk.Label(self, justify="center", anchor="n",
                                  text="Dataset & Extraction Parameters",
                                  font=("Helvatica", 10))
        self.labTitle1.grid(row=0, column=0, columnspan=2, padx=5, pady=3)
        self.labDataset = tk.Label(self, justify="left", text="Dataset:")
        self.labDataset.grid(row=1, column=0, padx=5, pady=3, sticky="W")
        self.entDataset = ttk.Entry(self, width=25, state="disabled")
        self.entDataset.grid(row=1, column=1,  padx=5, pady=3)
        self.labDataType = tk.Label(self, justify="left", text="Data Type:")
        self.labDataType.grid(row=2, column=0,  padx=5, pady=3, sticky="W")
        self.entDataType = ttk.Entry(self, width=25, state="disabled")
        self.entDataType.grid(row=2, column=1,  padx=5, pady=3)
        self.labFreqRng = tk.Label(self, justify="left", text=u"\u03bd Range:")
        self.labFreqRng.grid(row=3, column=0,  padx=5, pady=3, sticky="W")
        self.entFreqRng = ttk.Entry(self, width=25, state="disabled")
        self.entFreqRng.grid(row=3, column=1, padx=5, pady=3)
        self.labFreqChanWidth = tk.Label(self, justify="left", 
                                         text=u"\u03bd Chan. Width:")
        self.labFreqChanWidth.grid(row=4, column=0,  padx=5, pady=3, 
                                   sticky="W")
        self.entFreqChanWidth = ttk.Entry(self, width=25, state="disabled")
        self.entFreqChanWidth.grid(row=4, column=1,  padx=5, pady=3)
        self.labLamSqRng = tk.Label(self, justify="left",
                                    text=u"\u03bb\u00b2 Range:")
        self.labLamSqRng.grid(row=5, column=0,  padx=5, pady=3, sticky="W")
        self.entLamSqRng = ttk.Entry(self, width=25, state="disabled")
        self.entLamSqRng.grid(row=5, column=1,  padx=5, pady=3)
        self.labNFreqChan = tk.Label(self, justify="left",
                                     text=u"Num. \u03bd Channels:")
        self.labNFreqChan.grid(row=6, column=0,  padx=5, pady=3, sticky="W")
        self.entNFreqChan = ttk.Entry(self, width=25, state="disabled")
        self.entNFreqChan.grid(row=6, column=1,  padx=5, pady=3)
        self.labBoxScale = tk.Label(self, justify="left", text="Extract Box:")
        self.labBoxScale.grid(row=7, column=0,  padx=5, pady=3, sticky="W")
        self.entBoxScale = ttk.Entry(self, width=25, state="disabled")
        self.entBoxScale.grid(row=7, column=1,  padx=5, pady=3)
        sep = ttk.Separator(self, orient="vertical")
        sep.grid(row=0, column=2, rowspan=8, padx=5, pady=5, sticky="NS")
        self.columnconfigure(2, weight=1)
        
        # Layout the RM-Synthesis & RM-Clean column
        self.labTitle2 = tk.Label(self, justify="center", anchor="n",
                                  text="RM-Synthesis & RM-Clean Parameters",
                                  font=("Helvatica", 10))
        self.labTitle2.grid(row=0, column=3, columnspan=2, padx=5, pady=3),
        self.labPhiRng = tk.Label(self, justify="left", text=u"\u0278 Range:")
        self.labPhiRng.grid(row=1, column=3,  padx=5, pady=3, sticky="W")
        self.entPhiRng = ttk.Entry(self, width=25, state="disabled")
        self.entPhiRng.grid(row=1, column=4,  padx=5, pady=3)
        self.labPhiChanWidth = tk.Label(self, justify="left", 
                                        text=u"\u0278 Chan. Width:")
        self.labPhiChanWidth.grid(row=2, column=3,  padx=5, pady=3, sticky="W")
        self.entPhiChanWidth = ttk.Entry(self, width=25, state="disabled")
        self.entPhiChanWidth.grid(row=2, column=4,  padx=5, pady=3)
        self.labNPhiChan = tk.Label(self, justify="left",
                                    text=u"Num. \u0278 Channels:")
        self.labNPhiChan.grid(row=3, column=3,  padx=5, pady=3, sticky="W")
        self.entNPhiChan = ttk.Entry(self, width=25, state="disabled")
        self.entNPhiChan.grid(row=3, column=4,  padx=5, pady=3)
        self.labWtType = tk.Label(self, justify="left", text=u"Weighting:")
        self.labWtType.grid(row=4, column=3,  padx=5, pady=3, sticky="W")
        self.entWtType = ttk.Entry(self, width=25, state="disabled")
        self.entWtType.grid(row=4, column=4,  padx=5, pady=3)
        self.labClnCut = tk.Label(self, justify="left", text=u"Clean Cutoff:")
        self.labClnCut.grid(row=5, column=3,  padx=5, pady=3, sticky="W")
        self.entClnCut = ttk.Entry(self, width=25, state="disabled")
        self.entClnCut.grid(row=5, column=4,  padx=5, pady=3)
        self.labClnGain = tk.Label(self, justify="left", text=u"Clean Gain:")
        self.labClnGain.grid(row=6, column=3,  padx=5, pady=3, sticky="W")
        self.entClnGain = ttk.Entry(self, width=25, state="disabled")
        self.entClnGain.grid(row=6, column=4,  padx=5, pady=3)
        self.labMaxIter = tk.Label(self, justify="left", 
                                   text=u"Max. Iterations:")
        self.labMaxIter.grid(row=7, column=3,  padx=5, pady=2, sticky="W")
        self.entMaxIter = ttk.Entry(self, width=25, state="disabled")
        self.entMaxIter.grid(row=7, column=4,  padx=5, pady=3)
        sep = ttk.Separator(self, orient="vertical")
        sep.grid(row=0, column=5, rowspan=8, padx=5, pady=5, sticky="NS")
        self.columnconfigure(5, weight=1)
        
        # Layout the Processing & Thresholds column
        self.labTitle3 = tk.Label(self, justify="center", anchor="n",
                                  text="Processing & Flagging Thresholds",
                                  font=("Helvatica", 10))
        self.labTitle3.grid(row=0, column=6, columnspan=2, padx=5, pady=3),
        self.labDetectT = tk.Label(self, justify="left", text=u"Detection:")
        self.labDetectT.grid(row=1, column=6,  padx=5, pady=3, sticky="W")
        self.entDetectT = ttk.Entry(self, width=25, state="disabled")
        self.entDetectT.grid(row=1, column=7,  padx=5, pady=3)
        self.labDebiasT = tk.Label(self, justify="left", text=u"De-bias:")
        self.labDebiasT.grid(row=2, column=6,  padx=5, pady=3, sticky="W")
        self.entDebiasT = ttk.Entry(self, width=25, state="disabled")
        self.entDebiasT.grid(row=2, column=7,  padx=5, pady=3)
        
        sep = ttk.Separator(self, orient="horizontal")
        sep.grid(row=8, column=0, columnspan=9, padx=5, pady=15, sticky="EW")
        
    def clear_entries(self, disable=False):
        """Clear all the entry boxes in the frame."""

        self.entDataset.delete(0,tk.END)
        self.entDataType.delete(0,tk.END)
        self.entFreqRng.delete(0,tk.END)
        self.entFreqChanWidth.delete(0,tk.END)
        self.entLamSqRng.delete(0,tk.END)
        self.entNFreqChan.delete(0,tk.END)
        self.entBoxScale.delete(0,tk.END)
        self.entPhiRng.delete(0,tk.END)
        self.entPhiChanWidth.delete(0,tk.END)
        self.entNPhiChan.delete(0,tk.END)
        self.entWtType.delete(0,tk.END)
        self.entClnCut.delete(0,tk.END)
        self.entClnGain.delete(0,tk.END)
        self.entMaxIter.delete(0,tk.END)
        self.entDetectT.delete(0,tk.END)
        self.entDebiasT.delete(0,tk.END)
        if disable:
            self.entDataset.configure(state="disabled")
            self.entDataType.configure(state="disabled")
            self.entFreqRng.configure(state="disabled")
            self.entFreqChanWidth.configure(state="disabled")
            self.entLamSqRng.configure(state="disabled")
            self.entNFreqChan.configure(state="disabled")
            self.entBoxScale.configure(state="disabled")
            self.entPhiRng.configure(state="disabled")
            self.entPhiChanWidth.configure(state="disabled")
            self.entNPhiChan.configure(state="disabled")
            self.entWtType.configure(state="disabled")
            self.entClnCut.configure(state="disabled")
            self.entClnGain.configure(state="disabled")
            self.entMaxIter.configure(state="disabled")
            self.entDetectT.configure(state="disabled")
            self.entDebiasT.configure(state="disabled")

    def set_entries(self, pDict):
        """Insert values from a dictionary into each entry box."""

        self.entDataset.configure(state="enabled")
        self.entDataset.insert(0, str(pDict["dataPath"]))
        self.entDataType.configure(state="enabled")
        self.entDataType.insert(0, str(pDict["dataType"]))        
        tmpStr = u"%.3f \u2192 %.3f GHz" % (pDict["freqArr_Hz"][0]/1e9, 
                                            pDict["freqArr_Hz"][-1]/1e9)
        self.entFreqRng.configure(state="enabled")
        self.entFreqRng.insert(0, tmpStr)
        tmpStr = "%.2f MHz" % (pDict["dFreq_Hz"]/1e6)
        self.entFreqChanWidth.configure(state="enabled")
        self.entFreqChanWidth.insert(0, tmpStr)
        tmpStr = u"%.3f \u2192 %.3f m\u00b2" % (pDict["lambdaSqArr_m2"][0], 
                                                pDict["lambdaSqArr_m2"][-1])
        self.entLamSqRng.configure(state="enabled")
        self.entLamSqRng.insert(0, tmpStr)
        self.entNFreqChan.configure(state="enabled")
        self.entNFreqChan.insert(0, str(pDict["nChanFreq"]))
        tmpStr = u"%.1f \u2192 %.1f rad/m\u00b2" % (pDict["phiArr_radm2"][0], 
                                                    pDict["phiArr_radm2"][-1])
        self.entBoxScale.configure(state="enabled")
        self.entBoxScale.insert(0, str(pDict["sumBoxPix"]) + " pixels")
        self.entPhiRng.configure(state="enabled")
        self.entPhiRng.insert(0, tmpStr)
        tmpStr = u"%.1f rad/m\u00b2" % float(pDict["dPhi_radm2"])
        self.entPhiChanWidth.configure(state="enabled")
        self.entPhiChanWidth.insert(0, tmpStr)
        self.entNPhiChan.configure(state="enabled")
        self.entNPhiChan.insert(0, str(pDict["nChanRM"]))
        self.entWtType.configure(state="enabled")
        self.entWtType.insert(0, str(pDict["weightType"]))
        self.entClnCut.configure(state="enabled")
        self.entClnCut.insert(0, str(pDict["cleanCutoff_sigma"]) 
                              + u" \u03C3")
        self.entClnGain.configure(state="enabled")
        self.entClnGain.insert(0, str(pDict["gain"]))
        self.entMaxIter.configure(state="enabled")
        self.entMaxIter.insert(0, str(pDict["maxCleanIter"]))
        self.entDetectT.configure(state="enabled")
        self.entDetectT.insert(0, str(pDict["thresholdSignalPI_sigma"]) 
                               + u" \u03C3")
        self.entDebiasT.configure(state="enabled")
        self.entDebiasT.insert(0, str(pDict["thresholdPolBias_sigma"]) 
                               + u" \u03C3")

        
#-----------------------------------------------------------------------------#
class ResultsManager:
    """Class load the results from the database into memory."""

    def __init__(self, sessionPath):
        self.sessionPath = sessionPath
        self.catRec = None
        self.specParmRec = None
        self.dirtyFDFRec = None
        
        # Connect to the database
        dbFile = sessionPath + "/session.sqlite"
        conn = sqlite3.connect(dbFile)
        cursor = conn.cursor()

        # Load the tables into memory
        self._load_sourceCat(cursor)
        self._load_spectraParms(cursor)
        self._load_dirtyFDFparms(cursor)

        # Clean up
        cursor.close()
        conn.close()
        
    def _load_sourceCat(self, cursor):
        sql = "SELECT * FROM sourceCat"
        self.catRec = select_into_arr(cursor, sql)
        
    def _load_spectraParms(self, cursor):
        sql = "SELECT * FROM spectraParms"
        self.specParmRec = select_into_arr(cursor, sql)
        
    def _load_dirtyFDFparms(self, cursor):
        sql = "SELECT * FROM dirtyFDFparms"
        self.dirtyFDFRec = select_into_arr(cursor, sql)
        

#-----------------------------------------------------------------------------#
class SpectraPlotFrame(tk.Frame):
    """Class to show a summary of the input spectra for one source."""
    
    def __init__(self, parent, resultsManager, indx):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        self.frame = tk.Frame(self.parent)
        self.frame.grid(row=0, column=0, padx=5, pady=5, sticky="NSEW")
        self.resultsMan = resultsManager
        self.indx = indx
        self.uniqueName = self.resultsMan.catRec["uniqueName"][indx]
        
        # Title of table
        self.titleLab = tk.Label(self.frame, justify="left", anchor="nw",
                                 text="Stokes I, Q & U Spectra",
                                 font=("Helvatica", 10))
        self.titleLab.grid(row=0, column=0, padx=5, pady=3, sticky="NW")
        
        # Create the figure of the Stokes I spectrum and embed in a canvas
        #self.fig1 = plotSpecIerrs(self.resultsMan.sessionPath, self.uniqueName)
        self.fig1 = plotSpecParms(self.resultsMan.sessionPath, self.uniqueName)
        figCanvas = FigureCanvasTkAgg(self.fig1, master=self.frame)
        figCanvas.show()
        self.myCan = figCanvas.get_tk_widget()
        self.myCan.grid(row=1, column=0, padx=5, pady=5, sticky="NSEW")
        
        # Create the figure of the Stokes QU&V spectra and embed in a canvas
        #self.fig2 = plotSpecPQUerrs(self.resultsMan.sessionPath, 
        #                            self.uniqueName)
        #figCanvas = FigureCanvasTkAgg(self.fig2, master=self.frame)
        #figCanvas.show()
        #self.myCan = figCanvas.get_tk_widget()
        #self.myCan.grid(row=2, column=0, padx=5, pady=5, sticky="NSEW")
        
        # Create the figure of the RMS noise and embed in a canvas
        #self.fig3 = plotSpecRMS(self.resultsMan.sessionPath, self.uniqueName)
        #figCanvas = FigureCanvasTkAgg(self.fig3, master=self.frame)
        #figCanvas.show()
        #self.myCan = figCanvas.get_tk_widget()
        #self.myCan.grid(row=3, column=0, padx=5, pady=5, sticky="NSEW")

        # Configure the layout behaviour
        self.frame.columnconfigure(0, weight=1)
        self.frame.rowconfigure(1, weight=1)
        #self.frame.rowconfigure(2, weight=1)
        #self.frame.rowconfigure(3, weight=1)
        
#-----------------------------------------------------------------------------#
class RMsynthPlotFrame(tk.Frame):
    
    def __init__(self, parent, resultsManager, indx):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        self.frame = tk.Frame(self.parent)
        self.frame.grid(row=0, column=0, padx=5, pady=5, sticky="NSEW")
        self.resultsMan = resultsManager
        self.indx = indx
        self.uniqueName = self.resultsMan.catRec["uniqueName"][indx]

        # Title of table
        self.titleLab = tk.Label(self.frame, justify="left", anchor="nw",
                                 text="RMSF and Dirty FDF",
                                 font=("Helvatica", 10))
        self.titleLab.grid(row=0, column=0, padx=5, pady=3, sticky="NW")

        # Create the figure of the RMSF and dirty FDF and embed in a canvas
        self.fig1 = plotDirtyFDF(self.resultsMan.sessionPath, self.uniqueName)
        figCanvas = FigureCanvasTkAgg(self.fig1, master=self.frame)
        figCanvas.show()
        self.myCan = figCanvas.get_tk_widget()
        self.myCan.grid(row=1, column=0, padx=5, pady=5, sticky="NSEW")

        # Configure the layout behaviour
        self.frame.columnconfigure(0, weight=1)
        self.frame.rowconfigure(1, weight=1)
        

#-----------------------------------------------------------------------------#
if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    

