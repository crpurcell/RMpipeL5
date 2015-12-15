#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     rmPipeViewer.py                                                   #
#                                                                             #
# PURPOSE:  A graphical interface designed to view the results of the Level 5 #
#           RM-pipeline prototype.                                            #
#                                                                             #
# REQUIRED: Requires numpy, astropy, matplotlib                               #
#                                                                             #
# MODIFIED: 07-December-2015 by cpurcell                                      #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
# App               ... class to create the root and plotting window          #
# SessChooseFrame   ... frame to choose and load a session (results)          #
# NotebookFrame     ... frame holding the tabbed notebook interface           #
# PipeInputsFrame   ... frame to display the pipeline inputs                  #
# ResultsFrame      ... frame presenting the summary of results to the user   #
# DatabaseFrame     ... frame presenting the database tables to the user      #
# PlotSctHstFrame   ... frame presenting a plotting interface to the user     #
# SingleFigFrame    ... frame presenting pre-defined figures for a source     #
#                                                                             #
#=============================================================================#
#                                                                             #
# The MIT License (MIT)                                                       #
#                                                                             #
# Copyright (c) 2015 Cormac R. Purcell                                        #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the "Software"),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
#=============================================================================#

# Default session
defaultSessionDir = "testSessionImage"

# Window geometry
geometryBrowseWin = "1024x700"
geometryPlotWin = "900x800"

#-----------------------------------------------------------------------------#

import os
import sqlite3
import traceback
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
import Tkinter as tk
import ttk
import tkFileDialog as filedialog
import tkSimpleDialog
import tkMessageBox
import tkFont
from Imports.ScrolledTextTTK import ScrolledText

from Imports.util_PPC import DataManager
from Imports.util_PPC import PlotParms
from Imports.util_PPC import read_dictfile
from Imports.util_PPC import cleanup_str_input
from Imports.util_PPC import xfloat
from Imports.util_DB import *
from Imports.util_tk import *
from Imports.util_plotTk import *

# Turn off print statements buffering
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)


#-----------------------------------------------------------------------------#
class App:
    """
    Class defining the RM Pipeline Viewer application.
 
    This class creates a 'control' root window and a secondary TopLevel window,
    used for displaying plots and tables. The layout of individual sections in
    the GUI is done using other classes (e.g., DatabaseFrame) and instances of
    these classes are created and gridded into the two windows.

    All actions (e.g., selection, button-clicks) in sub-sections of the GUI
    generate virtual events with specific names. These virtual events are
    bound to event handler methods in this class. Virtual events cannot
    pass information, just an event object and name, so if additional
    information is needed (e.g., row number selected in a listbox) the event
    handler uses the 'event.widget.<method>' syntax to run <method> code
    belonging to the class which generated the event. This code can be
    configured to supply the necessary information (e.g., return a the
    currently selected row number).

    For example, when the button to show the Stokes I cutout image is clicked
    a handler method in the ResultsFrame class generates a virtual event called
    <<plot_stampI>>. This event is acted on by this App class (having
    been bound below using 'self.root.bind') and the 'self.on_show_preset_plot'
    method called with arguments 'event=<event_object>' and
    'plotType=plot_stampI'. This method calls the 'get_indx_selected()'
    method of the ResultsFrame class (via 'event.widget.get_indx_selected()')
    to determine which row is selected in the Source Summary Table. The method
    then clears the plotting window and creates an instance of the class
    'SingleFigFrame' to generate the Stokes I postage stamp image.

    Data access is performed through a DataManager class, a single instance of
    which is created when this App class is initiated. Methods in the
    DataManager object are used to access the FITS files on disk and perform
    queries on the SQLite3 database. Abstracting the data access through a
    single class simplifies modifying the data format in the future by
    requiring code changes only in one place. However, since the data
    storage model is based on one file per source, there is an I/O penalty
    when looping over a catalogue of sources (an open & close operation per
    iteration). This should not be a problem for the current usage model.
    """

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
        self.NBFrm = NotebookFrame(self.root)
        self.NBFrm.grid(row=1, column=0, padx=5, pady=5, sticky="NSEW")
        
        # Create a message panel
        self.msgFrm = MessageFrame(self.root)
        self.msgFrm.grid(row=2, column=0, padx=5, pady=5, sticky="NSEW")
        
        # Allow the notebook row and column to expand
        self.root.rowconfigure(1, weight=1)
        self.root.columnconfigure(0, weight=1)
        
        # Bind virtual events generated by sub-widgets
        self.root.bind("<<load_session>>", self.on_load_session)
        self.root.bind("<<plot_stampI>>", lambda event, 
                       plotType="plot_stampI" : 
                       self.on_show_preset_plot(event, plotType))
        self.root.bind("<<plot_stampP>>", lambda event, 
                       plotType="plot_stampP" :
                       self.on_show_preset_plot(event, plotType))
        self.root.bind("<<plot_IPQU>>", lambda event, 
                       plotType="plot_IPQU" : 
                       self.on_show_preset_plot(event, plotType))
        self.root.bind("<<plot_IQUrms>>", lambda event, 
                       plotType="plot_IQUrms" : 
                       self.on_show_preset_plot(event, plotType))
        self.root.bind("<<plot_polang>>", lambda event, 
                       plotType="plot_polang" : 
                       self.on_show_preset_plot(event, plotType))
        self.root.bind("<<plot_qup>>", lambda event, 
                       plotType="plot_qup" : 
                       self.on_show_preset_plot(event, plotType))
        self.root.bind("<<plot_q_vs_u>>", lambda event, 
                       plotType="plot_q_vs_u" : 
                       self.on_show_preset_plot(event, plotType))
        self.root.bind("<<plot_rmsf>>", lambda event,
                       plotType="plot_rmsf" : 
                       self.on_show_preset_plot(event, plotType))
        self.root.bind("<<plot_polsum>>", lambda event, 
                       plotType="plot_polsum" : 
                       self.on_show_preset_plot(event, plotType))
        self.root.bind("<<plot_dirty_fdf>>", lambda event, 
                       plotType="plot_dirty_fdf" : 
                       self.on_show_preset_plot(event, plotType))
        self.root.bind("<<plot_clean_fdf>>", lambda event, 
                       plotType="plot_clean_fdf" : 
                       self.on_show_preset_plot(event, plotType))        
        self.root.bind("<<view_sql_table>>", lambda event :
                       self.on_view_sql_table(event))  
        self.root.bind("<<export_sql_table>>", lambda event :
                       self.on_export_sql_table(event))
        self.root.bind("<<run_custom_sql>>", lambda event :
                       self.on_run_custom_sql(event))
        self.root.bind("<<plot_query_scthst>>", lambda event :
                       self.on_plot_query(event))
        
        # Create the visualisation window and set the focus back to root
        self.visWin = tk.Toplevel(self.root)
        self.visWin.title("RM Pipeline Viewer - Plotting Window1")
        self.visWin.geometry(geometryPlotWin)
        self.visWin.resizable(True, True)
        self.visWin.protocol("WM_DELETE_WINDOW", self.applicationExit)        
        self.visWin.columnconfigure(0, weight=1)
        self.visWin.rowconfigure(0, weight=1)

        # Create a master frame for the visualisation window
        self.visWinFrm = tk.Frame(self.visWin)
        self.visWinFrm.grid(row=0, column=0, padx=5, pady=5, sticky="NSEW")
        self.visWinFrm.columnconfigure(0, weight=1)
        self.visWinFrm.rowconfigure(0, weight=1)
        self.tmpVisLab = tk.Label(self.visWinFrm, font=("Helvatica", 30),
                                  text="Plotting Window")
        self.tmpVisLab.grid(row=0, column=0, padx=0, pady=0, sticky="NSEW")

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
        self.visWinFrm.destroy()
        self.visWinFrm = tk.Frame(self.visWin)        

    def applicationExit(self):
        """Exit the application cleanly."""
        self.root.destroy()

    # Event handlers ---------------------------------------------------------#
    #
    # on_load_session
    # on_show_preset_plot
    # on_view_sql_table
    # on_export_sql_table
    # on_run_custom_sql
    # on_plot_query
    #
    #
    #
    
    def on_load_session(self, event=None):
        """Load in the chosen pipeline session into the GUI."""
        
        # Check that the session exists and has an input file
        sessionPath = self.sessFrm.sesDir.get()
        if not os.path.exists(sessionPath):
            errStr = "Missing session directory.\n\nCheck that the " + \
                     "path to the session directory is correct.\n\n" + \
                     "[...%s]" % sessionPath[-27:]
            tkMessageBox.showerror("Error", errStr)
            return
        inputFile = sessionPath + "/inputs.config"
        if not os.path.exists(inputFile):
            errStr = "Missing input file.\n\nThe file 'inputs.config' " + \
                "is missing from the session directory.\n\nCheck that " + \
                "the chosen directory is a valid pipeline session."
            tkMessageBox.showerror("Error", errStr)
            return
        
        # Read the pipeline status file and set the status lights
        statusFile = sessionPath + "/status.json"
        statusDict = read_dictfile(statusFile)
        if int(statusDict["session"])<1:
            errStr = "Invalid session.\n\nSession status file reports " + \
                "that the session  was not created successfully.\n\n" + \
                "Check for errors when running 'create_session.py'."
            tkMessageBox.showerror("Error", errStr)
            return
        for key, val in statusDict.iteritems():
            self.sessFrm.set_status(process=key, status=int(val))

        
        # Create a DataManager instance to interface with DB and data
        # The DataManager loads the pipeline inputs, calculates secondary
        # parameters and loads a summary table of results to memory
        self.dataMan = DataManager(sessionPath)
                
        # Set the pipeline input parameters in the GUI
        self.NBFrm.inputsNB.clear_entries()
        self.NBFrm.inputsNB.set_entries(self.dataMan.pDict)
        
        # Load the summary table into the GUI
        self.NBFrm.resNB.load(self.dataMan.summaryRec)
        
        # Load the database structure into the GUI
        self.NBFrm.dbNB.load(self.dataMan)

        # Reset the plotting window
        self.reset_plotting_window()
        self.msgFrm.insert("Session loaded [%s]" % sessionPath)

    def on_show_preset_plot(self, event=None, plotType=""):
        """Create the plot requested by the user clicking a button."""

        # Get the index of the row selected
        indx = event.widget.get_indx_selected()
        if indx is None:
            errStr = "Please select a row in the table before plotting."
            tkMessageBox.showinfo("Error", errStr)
            return
        
        # Create the chosen plot
        if plotType=="plot_IPQU":
            titleStr = "Extracted Stokes I, Q & U Spectra:"
            fig =  plotSpecIPQU(self.dataMan, indx)
        elif plotType=="plot_IQUrms":
            titleStr = "Measured RMS Noise Spectra:"
            fig = plotSpecRMS(self.dataMan, indx)
        elif plotType=="plot_polang":
            titleStr = "Polarisation Angle versus Wavelength-Squared:"
            fig = plotPolang(self.dataMan, indx)
        elif plotType=="plot_qup":
            titleStr = "Polarisation Fraction versus Wavelength-Squared:"
            fig = plotFracPol(self.dataMan, indx)
        elif plotType=="plot_q_vs_u":
            titleStr = "Fractional Stokes Q versus U:"
            fig = plotFracQvsU(self.dataMan, indx)
        elif plotType=="plot_polsum":
            titleStr = "Polarisation Summary Plots:"
            fig = plotPolsummary(self.dataMan, indx)
        elif plotType=="plot_rmsf":
            titleStr = "Rotation Measure Spread Function:"
            fig = plotRMSF(self.dataMan, indx)
        elif plotType=="plot_dirty_fdf":
            titleStr = "Dirty FDF"
            fig = plotDirtyFDF(self.dataMan, indx)
        elif plotType=="plot_clean_fdf":
            titleStr = "Clean FDF and CC spectrum"
            fig = plotCleanFDF(self.dataMan, indx)
        elif plotType=="plot_stampI":
            titleStr = "Stokes I postage stamp image (channel 1)"
            fig = plotStampI(self.dataMan, indx)
        elif plotType=="plot_stampP":
            titleStr = "Polarised intensity P postage stamp image (channel 1)"
            fig = plotStampP(self.dataMan, indx)
        else:
            return
        
        # Clear the plotting window and create the requested plot
        self.reset_plotting_window()
        resultFrm = SingleFigFrame(self.visWinFrm, fig, titleStr)
        
        # Grid the plot in the visualisation window
        resultFrm.grid(row=0, column=0, padx=5, pady=5, sticky="NSEW")
        self.visWinFrm.columnconfigure(0, weight=1)
        self.visWinFrm.rowconfigure(0, weight=1)
        self.visWinFrm.grid(row=0, column=0, padx=0, pady=0, sticky="NSEW")

    def on_view_sql_table(self, event=None):
        """View the results of an SQL query in the plotting window"""

        # Run the query
        try:
            sql = event.widget.sql
            tableName = event.widget.tableName
            rowLimit = event.widget.rowLimit
            if sql is None or sql=="":
                return
            self.dataMan.query_database(sql)
        except Exception:
            return

        # Push results into GUI
        self.reset_plotting_window()
        titleStr = "Table '%s' (%s rows):" % (tableName, rowLimit)
        footerStr = "Click on a column header to sort up or down."
        resultFrm = SingleTabFrame(self.visWinFrm, titleStr, footerStr)
        resultFrm.table.insert_recarray(self.dataMan.tempRec)
        
        # Grid the table in the visualisation window
        resultFrm.grid(row=0, column=0, padx=5, pady=5, sticky="NSEW")
        self.visWinFrm.columnconfigure(0, weight=1)
        self.visWinFrm.rowconfigure(0, weight=1)
        self.visWinFrm.grid(row=0, column=0, padx=0, pady=0, sticky="NSEW")

    def on_export_sql_table(self, event=None):
        """Save the result of a SQL query"""

        # Get the new table name and save directory from the GUI
        tableName = event.widget.tableName
        outDir = event.widget.exp1Dir.get()
        if not os.path.exists(outDir):
            errStr = "Directory does not existL '%s'." % outDir
            tkMessageBox.showinfo("Error", errStr)
            return
        formatStr = event.widget.fmt1Comb.get()
        self.dataMan.export_table(tableName, formatStr, outDir)

    def on_run_custom_sql(self, event=None):
        """Run a custom SQL query"""
        
        # Fetch the query and output choice from the widget
        sql = event.widget.sql
        if sql is None or sql=="":
            errStr = "Please enter a query before clicking 'Go'."
            tkMessageBox.showinfo("Error", errStr)
            return
        outType = event.widget.sqlOutChoice.get()
        formatStr = event.widget.fmt2Comb.get()
        outDir = event.widget.exp2Dir.get()
        if not os.path.exists(outDir):
            errStr = "Directory does not existL '%s'." % outDir
            tkMessageBox.showinfo("Error", errStr)
            return
        newTableName = event.widget.table1Ent.get()
        if newTableName is None or newTableName=="":                
            errStr = "Please enter a new table name in the box."
            tkMessageBox.showinfo("Error", errStr)
            return
        # Prepend a 'CREATE TABLE' statement if creating a new table
        if outType=="db":
            prep = "CREATE TABLE IF NOT EXISTS %s AS " % newTableName
            sql = prep + sql
        
        # Run the query
        try:
            if sql is None or sql=="":
                return
            self.dataMan.query_database(sql, buffer=True)
        except Exception:
            return

        # Process the results
        if outType=="db":
            event.widget.load(self.dataMan) # Reload the schema
        elif outType=="exp":
            self.dataMan.export_table(newTableName, formatStr, outDir, True)
        else:
            
            # Display the results in the plotting window
            self.reset_plotting_window()
            titleStr = "Results of custom SQL Query"
            footerStr = "Click on a column header to sort up or down."
            resultFrm = SingleTabFrame(self.visWinFrm, titleStr, footerStr)
            resultFrm.table.insert_recarray(self.dataMan.tempRec)
            resultFrm.grid(row=0, column=0, padx=5, pady=5, sticky="NSEW")
            self.visWinFrm.columnconfigure(0, weight=1)
            self.visWinFrm.rowconfigure(0, weight=1)
            self.visWinFrm.grid(row=0, column=0, padx=0, pady=0, sticky="NSEW")

    def on_plot_query(self, event=None):
        """Plot the results of a user defined query."""

        # Refresh the plot parameter object from the GUI & create the plot
        event.widget.update_plotparm_obj()
        fig = plotSctHstQuery(self.dataMan, event.widget.plotParm)

        # Clear the plotting window and create the requested plot
        self.reset_plotting_window()
        titleStr = "Custom Query Plot:"
        resultFrm = SingleFigFrame(self.visWinFrm, fig, titleStr)
        
        # Grid the plot in the visualisation window
        resultFrm.grid(row=0, column=0, padx=5, pady=5, sticky="NSEW")
        self.visWinFrm.columnconfigure(0, weight=1)
        self.visWinFrm.rowconfigure(0, weight=1)
        self.visWinFrm.grid(row=0, column=0, padx=0, pady=0, sticky="NSEW")


#-----------------------------------------------------------------------------#
class SessChooseFrame(tk.Frame):
    """Frame allowing the user to choose which session to load and to display
    'lights' showing the status of the modules in the pipeline."""
    
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        
        # Create the session chooser frame
        self.chooseFrm = ttk.LabelFrame(self, text="  Session Chooser  ")
        self.chooseFrm.grid(column=0, row=0, padx=5, pady=2, sticky="NSEW")
        self.instLab = tk.Label(self.chooseFrm, anchor="e",
               text="Enter or browse to a session directory and click 'Load':")
        self.instLab.grid(row=0, column=0, padx=5, pady=5, sticky="NW")
        self.sesDir = tk.StringVar()
        self.entSes = ttk.Entry(self.chooseFrm, textvariable=self.sesDir)
        self.entSes.grid(row=1, column=0, padx=5, pady=5, sticky="EW")
        self.browseBtn = ttk.Button(self.chooseFrm, text="Browse",
                                    command=self._handlerBrowseButton)
        self.browseBtn.grid(row=0, column=1, padx=5, pady=2,sticky="NW" )
        self.btnLoad = ttk.Button(self.chooseFrm, text="Load",
                                  command=self._handlerLoadButton)
        self.btnLoad.grid(row=1, column=1, padx=5, pady=5, sticky="NW" )
        
        # Create the session status frame
        self.statFrm = ttk.LabelFrame(self, text="  Session Status  ")
        self.statFrm.grid(column=1, row=0, padx=5, pady=2, sticky="NSEW")
        self.extractStatLab = ttk.Label(self.statFrm, justify="center",
                                        foreground="black",
                                        background="orange",
                                        text="Waiting", relief="solid",
                                        anchor="center", padding=5, width=10)
        self.extractStat = ttk.Label(self.statFrm, justify="center",
                                         text="Extract Spectra",
                                         anchor="center", padding=0)
        self.extractStatLab.grid(row=0, column=0, padx=15, pady=5,
                                 sticky="S")
        self.extractStat.grid(row=1, column=0, padx=15, pady=5,
                              sticky="N")        
        self.rmsynthStatLab = ttk.Label(self.statFrm, justify="center",
                                        foreground="black",
                                        background="orange",
                                        text="Waiting", relief="solid",
                                        anchor="center", padding=5, width=10)
        self.rmsynthStat = ttk.Label(self.statFrm, justify="center",
                                     text="RM-Synthesis",
                                     anchor="center", padding=0)
        self.rmsynthStatLab.grid(row=0, column=1, padx=15, pady=5,
                                  sticky="S")
        self.rmsynthStat.grid(row=1, column=1, padx=15, pady=5,
                              sticky="N")
        
        self.rmcleanStatLab = ttk.Label(self.statFrm, justify="center",
                                        foreground="black",
                                        background="orange",
                                        text="Waiting", relief="solid",
                                        anchor="center", padding=5, width=10)
        self.rmcleanStat = ttk.Label(self.statFrm, justify="center",
                                     text="RM-CLEAN",
                                     anchor="center", padding=0)
        self.rmcleanStatLab.grid(row=0, column=2, padx=5, pady=5,
                                  sticky="S")
        self.rmcleanStat.grid(row=1, column=2, padx=15, pady=5,
                              sticky="N")
        
        self.cmplxStatLab = ttk.Label(self.statFrm, justify="center",
                                      foreground="black",
                                      background="orange",
                                      text="Waiting", relief="solid",
                                      anchor="center", padding=5, width=10)
        self.cmplxStat = ttk.Label(self.statFrm, justify="center",
                                   text="Complexity",
                                   anchor="center", padding=0)
        self.cmplxStatLab.grid(row=0, column=3, padx=5, pady=5,
                               sticky="S")
        self.cmplxStat.grid(row=1, column=3, padx=15, pady=5,
                            sticky="N")
        
        # Set the expansion behaviour
        self.chooseFrm.columnconfigure(0, weight=1)
        self.statFrm.columnconfigure(0, weight=1)
        self.statFrm.columnconfigure(1, weight=1)
        self.statFrm.columnconfigure(2, weight=1)
        self.statFrm.columnconfigure(3, weight=1)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
    
    def _handlerBrowseButton(self):
        """Open the directory selection dialog and set the session dir."""
        sesDir = filedialog.askdirectory(parent=self, initialdir=".",
                                    title="Please select a session directory")
        if not len(sesDir)==0:
            self.sesDir.set(sesDir)
            
    def _handlerLoadButton(self):      
        """Raise a <<load_session>> virtual event in the parent window."""
        self.event_generate("<<load_session>>")

    def set_status(self, process="extract", status=0):

        """Set the colour and text in the pipeline status 'lights'"""        
        if  process=="extract":
            targetLab = self.extractStatLab
        elif  process=="rmsynth":
            targetLab = self.rmsynthStatLab
        elif  process=="rmclean":
            targetLab = self.rmcleanStatLab
        elif  process=="complexity":
            targetLab = self.cmplxStatLab
        else:
            return
        if status==2:            
            targetLab.config(background="orange")
            targetLab.config(foreground="black")
            targetLab.config(text="NA")
        if status==1:            
            targetLab.config(background="green")
            targetLab.config(foreground="black")
            targetLab.config(text="OK")
        if status==0:            
            targetLab.config(background="red")
            targetLab.config(foreground="white")
            targetLab.config(text="Not OK")


#-----------------------------------------------------------------------------#
class NotebookFrame(tk.Frame):
    """
    Frame containing the a tabbed 'notebook' used to layout the sections of
    the graphical user interface: [inputs, results, database].
    """

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
        self._add_results_tab()
        self._add_database_tab()
        self._add_plotting_tab()
        
    def _add_pipeinputs_tab(self):
        tabFrm = tk.Frame(self.nb)
        tabFrm.rowconfigure(0, weight=1)
        tabFrm.columnconfigure(0, weight=1)
        self.inputsNB = PipeInputsFrame(tabFrm)
        self.inputsNB.grid(column=0, row=0, sticky="NSEW")
        self.nb.add(tabFrm, text="  Pipeline Inputs  ", padding=5)
        
    def _add_results_tab(self):
        tabFrm = tk.Frame(self.nb)
        tabFrm.rowconfigure(0, weight=1)
        tabFrm.columnconfigure(0, weight=1)
        self.resNB = ResultsFrame(tabFrm)
        self.resNB.grid(column=0, row=0, sticky="NSEW")
        self.nb.add(tabFrm, text="  Results by Source  ", padding=5)

    def _add_database_tab(self):
        tabFrm = tk.Frame(self.nb)
        tabFrm.rowconfigure(0, weight=1)
        tabFrm.columnconfigure(0, weight=1)
        self.dbNB = DatabaseFrame(tabFrm)
        self.dbNB.grid(column=0, row=0, sticky="NSEW")
        self.nb.add(tabFrm, text="  Database Query & Export  ", padding=5)

    def _add_plotting_tab(self):
        tabFrm = tk.Frame(self.nb)
        tabFrm.rowconfigure(0, weight=1)
        tabFrm.columnconfigure(0, weight=1)
        self.plotNB = PlotSctHstFrame(tabFrm)
        self.plotNB.grid(column=0, row=0, sticky="NSEW")
        self.nb.add(tabFrm, text="  Scatter & Histogram Plotting  ", padding=5)

        
#-----------------------------------------------------------------------------#
class MessageFrame(tk.Frame):
    """Frame presenting a message box to the user."""

    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        
        self.msgTextBox = ScrolledText(self, height=3, foreground="white",
                                       background="black")
        self.msgTextBox.grid(row=0, column=0, padx=5, pady=5, sticky="EW")
        self.msgTextBox.insert("1.0", "Pipeline Viewer Messages:\n")
        self.msgTextBox.configure(state="disabled")
        
        # Set the expansion behaviour
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        # Set focus so that the user can copy text
        self.msgTextBox.bind("<1>", lambda event: self.msgTextBox.focus_set())
    def insert(self, msg):        
        self.msgTextBox.configure(state="normal")
        self.msgTextBox.insert("end", msg + "\n")
        self.msgTextBox.see(tk.END)
        self.msgTextBox.configure(state="disabled")
    

#-----------------------------------------------------------------------------#
class PipeInputsFrame(tk.Frame):
    """Frame presenting the pipeline inputs to the user.
    
    TODO: Greatly simplify this interface by presenting the pipeline input
    file on the left and derived parameters on the right. The pipeline input
    file is already well formatted and so could be displayed using syntax
    highlighting.    
    """
    
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        
        # Layout the Dataset & Extraction column
        self.title1Lab = tk.Label(self, justify="center", anchor="n",
                                  text="Dataset & Extraction Parameters")
        self.title1Lab.grid(row=0, column=0, columnspan=2, padx=5, pady=3)
        self.datasetLab = tk.Label(self, justify="left", text="Dataset:")
        self.datasetLab.grid(row=1, column=0, padx=5, pady=3, sticky="W")
        self.datasetEnt = ttk.Entry(self, width=20, state="disabled")
        self.datasetEnt.grid(row=1, column=1,  padx=5, pady=3)
        self.datatypeLab = tk.Label(self, justify="left", text="Data Type:")
        self.datatypeLab.grid(row=2, column=0,  padx=5, pady=3, sticky="W")
        self.datatypeEnt = ttk.Entry(self, width=20, state="disabled")
        self.datatypeEnt.grid(row=2, column=1,  padx=5, pady=3)
        self.freqrngLab = tk.Label(self, justify="left", text=u"\u03bd Range:")
        self.freqrngLab.grid(row=3, column=0,  padx=5, pady=3, sticky="W")
        self.freqrngEnt = ttk.Entry(self, width=20, state="disabled")
        self.freqrngEnt.grid(row=3, column=1, padx=5, pady=3)
        self.freqChanWidthLab = tk.Label(self, justify="left", 
                                         text=u"\u03bd Chan. Width:")
        self.freqChanWidthLab.grid(row=4, column=0,  padx=5, pady=3, 
                                   sticky="W")
        self.freqChanWidthEnt = ttk.Entry(self, width=20, state="disabled")
        self.freqChanWidthEnt.grid(row=4, column=1,  padx=5, pady=3)
        self.lamSqRngLab = tk.Label(self, justify="left",
                                    text=u"\u03bb\u00b2 Range:")
        self.lamSqRngLab.grid(row=5, column=0,  padx=5, pady=3, sticky="W")
        self.lamSqRngEnt = ttk.Entry(self, width=20, state="disabled")
        self.lamSqRngEnt.grid(row=5, column=1,  padx=5, pady=3)
        self.nFreqChanLab = tk.Label(self, justify="left",
                                     text=u"Num. \u03bd Channels:")
        self.nFreqChanLab.grid(row=6, column=0,  padx=5, pady=3, sticky="W")
        self.nFreqChanEnt = ttk.Entry(self, width=20, state="disabled")
        self.nFreqChanEnt.grid(row=6, column=1,  padx=5, pady=3)
        self.boxScaleLab = tk.Label(self, justify="left", text="Extract Box:")
        self.boxScaleLab.grid(row=7, column=0,  padx=5, pady=3, sticky="W")
        self.boxScaleEnt = ttk.Entry(self, width=20, state="disabled")
        self.boxScaleEnt.grid(row=7, column=1,  padx=5, pady=3)
        
        sep = ttk.Separator(self, orient="vertical")
        sep.grid(row=0, column=2, rowspan=8, padx=5, pady=5, sticky="NS")
        self.columnconfigure(2, weight=1)
        
        # Layout the RM-Synthesis & RM-Clean column
        self.title2Lab = tk.Label(self, justify="center", anchor="n",
                                  text="RM-Synthesis & RM-Clean Parameters")
        self.title2Lab.grid(row=0, column=3, columnspan=2, padx=5, pady=3),
        self.phiRngLab = tk.Label(self, justify="left", text=u"\u0278 Range:")
        self.phiRngLab.grid(row=1, column=3,  padx=5, pady=3, sticky="W")
        self.phiRngEnt = ttk.Entry(self, width=20, state="disabled")
        self.phiRngEnt.grid(row=1, column=4,  padx=5, pady=3)
        self.phiChanWidthLab = tk.Label(self, justify="left", 
                                        text=u"\u0278 Chan. Width:")
        self.phiChanWidthLab.grid(row=2, column=3,  padx=5, pady=3, sticky="W")
        self.phiChanWidthEnt = ttk.Entry(self, width=20, state="disabled")
        self.phiChanWidthEnt.grid(row=2, column=4,  padx=5, pady=3)
        self.nPhiChanLab = tk.Label(self, justify="left",
                                    text=u"Num. \u0278 Channels:")
        self.nPhiChanLab.grid(row=3, column=3,  padx=5, pady=3, sticky="W")
        self.nPhiChanEnt = ttk.Entry(self, width=20, state="disabled")
        self.nPhiChanEnt.grid(row=3, column=4,  padx=5, pady=3)
        self.wtTypeLab = tk.Label(self, justify="left", text=u"Weighting:")
        self.wtTypeLab.grid(row=4, column=3,  padx=5, pady=3, sticky="W")
        self.wtTypeEnt = ttk.Entry(self, width=20, state="disabled")
        self.wtTypeEnt.grid(row=4, column=4,  padx=5, pady=3)
        self.clnCutLab = tk.Label(self, justify="left", text=u"Clean Cutoff:")
        self.clnCutLab.grid(row=5, column=3,  padx=5, pady=3, sticky="W")
        self.clnCutEnt = ttk.Entry(self, width=20, state="disabled")
        self.clnCutEnt.grid(row=5, column=4,  padx=5, pady=3)
        self.clnGainLab = tk.Label(self, justify="left", text=u"Clean Gain:")
        self.clnGainLab.grid(row=6, column=3,  padx=5, pady=3, sticky="W")
        self.clnGainEnt = ttk.Entry(self, width=20, state="disabled")
        self.clnGainEnt.grid(row=6, column=4,  padx=5, pady=3)
        self.maxIterLab = tk.Label(self, justify="left", 
                                   text=u"Max. Iterations:")
        self.maxIterLab.grid(row=7, column=3,  padx=5, pady=2, sticky="W")
        self.maxIterEnt = ttk.Entry(self, width=20, state="disabled")
        self.maxIterEnt.grid(row=7, column=4,  padx=5, pady=3)
        
        sep = ttk.Separator(self, orient="vertical")
        sep.grid(row=0, column=5, rowspan=8, padx=5, pady=5, sticky="NS")
        self.columnconfigure(5, weight=1)
        
        # Layout the Processing & Thresholds column
        self.title3Lab = tk.Label(self, justify="center", anchor="n",
                                  text="Processing & Flagging Thresholds")
        self.title3Lab.grid(row=0, column=6, columnspan=2, padx=5, pady=3),
        self.detectTLab = tk.Label(self, justify="left", text=u"Detection:")
        self.detectTLab.grid(row=1, column=6,  padx=5, pady=3, sticky="W")
        self.detectTEnt = ttk.Entry(self, width=20, state="disabled")
        self.detectTEnt.grid(row=1, column=7,  padx=5, pady=3)
        self.debiasTLab = tk.Label(self, justify="left", text=u"De-bias:")
        self.debiasTLab.grid(row=2, column=6,  padx=5, pady=3, sticky="W")
        self.debiasTEnt = ttk.Entry(self, width=20, state="disabled")
        self.debiasTEnt.grid(row=2, column=7,  padx=5, pady=3)
        
        sep = ttk.Separator(self, orient="horizontal")
        sep.grid(row=8, column=0, columnspan=9, padx=5, pady=15, sticky="EW")
        
    def clear_entries(self, disable=False):
        """Clear all the entry boxes in the frame."""
        self.datasetEnt.delete(0,tk.END)
        self.datatypeEnt.delete(0,tk.END)
        self.freqrngEnt.delete(0,tk.END)
        self.freqChanWidthEnt.delete(0,tk.END)
        self.lamSqRngEnt.delete(0,tk.END)
        self.nFreqChanEnt.delete(0,tk.END)
        self.boxScaleEnt.delete(0,tk.END)
        self.phiRngEnt.delete(0,tk.END)
        self.phiChanWidthEnt.delete(0,tk.END)
        self.nPhiChanEnt.delete(0,tk.END)
        self.wtTypeEnt.delete(0,tk.END)
        self.clnCutEnt.delete(0,tk.END)
        self.clnGainEnt.delete(0,tk.END)
        self.maxIterEnt.delete(0,tk.END)
        self.detectTEnt.delete(0,tk.END)
        self.debiasTEnt.delete(0,tk.END)
        if disable:
            self.datasetEnt.configure(state="disabled")
            self.datatypeEnt.configure(state="disabled")
            self.freqrngEnt.configure(state="disabled")
            self.freqChanWidthEnt.configure(state="disabled")
            self.lamSqRngEnt.configure(state="disabled")
            self.nFreqChanEnt.configure(state="disabled")
            self.boxScaleEnt.configure(state="disabled")
            self.phiRngEnt.configure(state="disabled")
            self.phiChanWidthEnt.configure(state="disabled")
            self.nPhiChanEnt.configure(state="disabled")
            self.wtTypeEnt.configure(state="disabled")
            self.clnCutEnt.configure(state="disabled")
            self.clnGainEnt.configure(state="disabled")
            self.maxIterEnt.configure(state="disabled")
            self.detectTEnt.configure(state="disabled")
            self.debiasTEnt.configure(state="disabled")

    def set_entries(self, pDict):
        """Insert values from a dictionary into each entry box."""

        self.datasetEnt.configure(state="enabled")
        self.datasetEnt.insert(0, str(pDict["dataPath"]))
        self.datatypeEnt.configure(state="enabled")
        self.datatypeEnt.insert(0, str(pDict["dataType"]))        
        tmpStr = u"%.3f \u2192 %.3f GHz" % (pDict["freqArr_Hz"][0]/1e9, 
                                            pDict["freqArr_Hz"][-1]/1e9)
        self.freqrngEnt.configure(state="enabled")
        self.freqrngEnt.insert(0, tmpStr)
        tmpStr = "%.2f MHz" % (pDict["dFreq_Hz"]/1e6)
        self.freqChanWidthEnt.configure(state="enabled")
        self.freqChanWidthEnt.insert(0, tmpStr)
        tmpStr = u"%.3f \u2192 %.3f m\u00b2" % (pDict["lambdaSqArr_m2"][0], 
                                                pDict["lambdaSqArr_m2"][-1])
        self.lamSqRngEnt.configure(state="enabled")
        self.lamSqRngEnt.insert(0, tmpStr)
        self.nFreqChanEnt.configure(state="enabled")
        self.nFreqChanEnt.insert(0, str(pDict["nChanFreq"]))
        tmpStr = u"%.1f \u2192 %.1f rad/m\u00b2" % (pDict["phiArr_radm2"][0], 
                                                    pDict["phiArr_radm2"][-1])
        self.boxScaleEnt.configure(state="enabled")
        self.boxScaleEnt.insert(0, str(pDict["sumBoxPix"]) + " pixels")
        self.phiRngEnt.configure(state="enabled")
        self.phiRngEnt.insert(0, tmpStr)
        tmpStr = u"%.1f rad/m\u00b2" % float(pDict["dPhi_radm2"])
        self.phiChanWidthEnt.configure(state="enabled")
        self.phiChanWidthEnt.insert(0, tmpStr)
        self.nPhiChanEnt.configure(state="enabled")
        self.nPhiChanEnt.insert(0, str(pDict["nChanRM"]))
        self.wtTypeEnt.configure(state="enabled")
        self.wtTypeEnt.insert(0, str(pDict["weightType"]))
        self.clnCutEnt.configure(state="enabled")
        self.clnCutEnt.insert(0, str(pDict["cleanCutoff_sigma"]) 
                              + u" \u03C3")
        self.clnGainEnt.configure(state="enabled")
        self.clnGainEnt.insert(0, str(pDict["gain"]))
        self.maxIterEnt.configure(state="enabled")
        self.maxIterEnt.insert(0, str(pDict["maxCleanIter"]))
        self.detectTEnt.configure(state="enabled")
        self.detectTEnt.insert(0, str(pDict["thresholdSignalPI_sigma"]) 
                               + u" \u03C3")
        self.debiasTEnt.configure(state="enabled")
        self.debiasTEnt.insert(0, str(pDict["thresholdPolBias_sigma"]) 
                               + u" \u03C3")
        
        
#-----------------------------------------------------------------------------#
class ResultsFrame(tk.Frame):
    """
    Frame presenting a summary table of results to the user and a list of
    further actions, e.g., plotting the spectra, FDF etc. for individual
    sources.""" 
    
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        self.rowSelected = None

        # Table Frame
        self.tableFrm = tk.Frame(self)
        self.tableFrm.grid(column=0, row=0, rowspan=13, padx=5, pady=2,
                           sticky="NSEW")
        self.title1Lab = tk.Label(self.tableFrm, justify="center", anchor="nw",
                                  text="Source Summary Table:")
        self.title1Lab.grid(row=0, column=0, columnspan=2, padx=5, pady=3,
                            sticky="EW")
        self.resultsTab = ScrolledTreeTab(self.tableFrm,
                                   virtEvent="<<cat_row_selected>>")
        self.resultsTab.grid(column=0, row=1, padx=5, pady=0, sticky="NSEW")
        self.footer1Lab = tk.Label(self.tableFrm, justify="center",
                                   anchor="nw",
                           text="Click on a column header to sort up or down.")
        self.footer1Lab.grid(row=2, column=0, columnspan=2, padx=5, pady=3,
                            sticky="EW")
        self.tableFrm.columnconfigure(0, weight=1)
        self.tableFrm.rowconfigure(1, weight=1)

        # Action frame -------------------------------------------------------#
        self.title2Lab = tk.Label(self, justify="center", anchor="nw",
                                  text="Plotting Actions:")
        self.title2Lab.grid(row=0, column=1, columnspan=2, padx=5, pady=3,
                            sticky="EW")
        # Sumary page
        self.act1Lab = tk.Label(self, justify="left",
                                text="Show polarised intensity cutout image")
        self.act1Lab.grid(row=1, column=1, padx=5, pady=3, sticky="W")
        self.act1Btn = ttk.Button(self, text="Go",
                                  command=lambda action="plot_stampP" :
                                  self._handlerButtonPress(action))
        self.act1Btn.grid(row=1, column=2, padx=5, pady=2,sticky="NW")
        # Postage stamp images
        self.act2Lab = tk.Label(self, justify="left",
                                text="Show Stokes I cutout image")
        self.act2Lab.grid(row=2, column=1, padx=5, pady=3, sticky="W")
        self.act2Btn = ttk.Button(self, text="Go",
                                  command=lambda action="plot_stampI" :
                                  self._handlerButtonPress(action))
        self.act2Btn.grid(row=2, column=2, padx=5, pady=2,sticky="NW")        
        # Stokes IQU spectra
        self.act3Lab = tk.Label(self, justify="left",
                                text=u"Plot Stokes I, P, Q & U versus \u03bd")
        self.act3Lab.grid(row=3, column=1, padx=5, pady=3, sticky="W")
        self.act3Btn = ttk.Button(self, text="Go",
                                  command=lambda action="plot_IPQU" :
                                  self._handlerButtonPress(action))
        self.act3Btn.grid(row=3, column=2, padx=5, pady=2,sticky="NW")
        # Stokes IQU RMS
        self.act4Lab = tk.Label(self, justify="left",
                                text=u"Plot I, Q & U RMS noise versus \u03bd")
        self.act4Lab.grid(row=4, column=1, padx=5, pady=3, sticky="W")
        self.act4Btn = ttk.Button(self, text="Go",
                                  command=lambda action="plot_IQUrms" :
                                  self._handlerButtonPress(action))
        self.act4Btn.grid(row=4, column=2, padx=5, pady=2,sticky="NW")
        # Polangle versus lambda^2
        self.act5Lab = tk.Label(self, justify="left",
                                text=u"Plot \u03a8 versus \u03bb\u00b2")
        self.act5Lab.grid(row=5, column=1, padx=5, pady=3, sticky="W")
        self.act5Btn = ttk.Button(self, text="Go",
                                  command=lambda action="plot_polang" :
                                  self._handlerButtonPress(action))
        self.act5Btn.grid(row=5, column=2, padx=5, pady=2,sticky="NW")
        # Fractional spectra ves lambda^2
        self.act6Lab = tk.Label(self, justify="left",
                         text=u"Plot fractional q, u & p versus \u03bb\u00b2")
        self.act6Lab.grid(row=6, column=1, padx=5, pady=3, sticky="W")
        self.act6Btn = ttk.Button(self, text="Go",
                                  command=lambda action="plot_qup" :
                                  self._handlerButtonPress(action))
        self.act6Btn.grid(row=6, column=2, padx=5, pady=2,sticky="NW")
        # Fractional q vs u
        self.act7Lab = tk.Label(self, justify="left",
                          text="Plot fractional q versus u")
        self.act7Lab.grid(row=7, column=1, padx=5, pady=3, sticky="W")
        self.act7Btn = ttk.Button(self, text="Go",
                                  command=lambda action="plot_q_vs_u" :
                                  self._handlerButtonPress(action))
        self.act7Btn.grid(row=7, column=2, padx=5, pady=2,sticky="NW")        
        # Polarisation summary
        self.act8Lab = tk.Label(self, justify="left",
             text=u"Plot 4-panel polarisation summary \n( I vs \u03bd,   " + \
             u"\u03a8 vs \u03bb\u00b2,   p|q|u vs \u03bb\u00b2  &  q vs u )")
        self.act8Lab.grid(row=8, column=1, padx=5, pady=3, sticky="W")
        self.act8Btn = ttk.Button(self, text="Go",
                                  command=lambda action="plot_polsum" :
                                  self._handlerButtonPress(action))
        self.act8Btn.grid(row=8, column=2, padx=5, pady=2,sticky="NW")
        # RMSF
        self.act9Lab = tk.Label(self, justify="left",
                                text="Plot RMSF & Gaussian fit to main lobe")
        self.act9Lab.grid(row=9, column=1, padx=5, pady=3, sticky="W")
        self.act9Btn = ttk.Button(self, text="Go",
                                  command=lambda action="plot_rmsf" :
                                  self._handlerButtonPress(action))
        self.act9Btn.grid(row=9, column=2, padx=5, pady=2,sticky="NW")
        # Dirty FDF
        self.act10Lab = tk.Label(self, justify="left",
                          text="Plot dirty Faraday Dispersion Function (FDF)")
        self.act10Lab.grid(row=10, column=1, padx=5, pady=3, sticky="W")
        self.act10Btn = ttk.Button(self, text="Go",
                                   command=lambda action="plot_dirty_fdf" :
                                   self._handlerButtonPress(action))
        self.act10Btn.grid(row=10, column=2, padx=5, pady=2,sticky="NW")
        # CLEAN FFD
        self.act11Lab = tk.Label(self, justify="left",
                           text="Plot CLEANed FDF & CLEAN-component spectrum")
        self.act11Lab.grid(row=11, column=1, padx=5, pady=3, sticky="W")
        self.act11Btn = ttk.Button(self, text="Go",
                                  command=lambda action="plot_clean_fdf" :
                                  self._handlerButtonPress(action))
        self.act11Btn.grid(row=11, column=2, padx=5, pady=2,sticky="NW")

        # Set the expansion behaviour
        self.columnconfigure(0, weight=1)
        self.rowconfigure(12, weight=1)

    def _handlerButtonPress(self, action):
        """Store the selected row index and generate a virtual event. The
        event is named for the action assoicated with each button above and
        will be processed by the parent object."""
        
        self.rowSelected = self.resultsTab.get_indx_selected()
        virtualEvent = "<<" + action + ">>"
        self.event_generate(virtualEvent)
       
    def get_indx_selected(self):
        """Return the index of the last row selected."""        
        if self.rowSelected is None:
            return None
        else:
            return int(self.rowSelected)

    def load(self, recArr):
        """Reset the summary table and load a recarray"""
        self.resultsTab.clear_entries()
        self.resultsTab.insert_recarray(recArr)
        self.rowSelected = None

   
#-----------------------------------------------------------------------------#
class DatabaseFrame(tk.Frame):
    """
    Frame presenting the schema of the SQLite database and the structure
    the individual tables. The interface provides the ability to export
    individual tables in CSV format or view the tables in the results window.

    The key feature here is the ability to execute custom SQL queries and
    save the results as a new table, export to disk or view in the plotting
    window.
    """
    
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        self.dataMan = None
        self.sql = None
        self.tableName = ""
        self.rowLimit = ""

        # TreeView showing table structures
        self.titleSchemaLab = tk.Label(self, anchor="nw",
                                  text="Database Structure:")
        self.titleSchemaLab.grid(row=0, column=0, padx=5, pady=3, sticky="EW")
        self.schemaFrm = tk.Frame(self)
        self.schemaFrm.grid(row=1, column=0, rowspan=8, padx=5, pady=2,
                            sticky="NSEW")
        self.schemaTree = ScrolledTreeView(self.schemaFrm,
                                         virtEvent="<<schema_table_selected>>")
        self.schemaTree.tree['columns'] = ("column", "type", "key")
        self.schemaTree.tree['displaycolumns'] = ("type", "key")
        self.schemaTree.tree.heading("#0", text="Table / Column", anchor='w')
        self.schemaTree.tree.heading("column", text="Column Name", anchor='w')
        self.schemaTree.tree.heading("type", text="Data Type", anchor='w')
        self.schemaTree.tree.heading("key", text="Primary Key", anchor='w')
        self.schemaTree.tree.column("#0", stretch=1)
        self.schemaTree.tree.column("type", stretch=0, width=100)
        self.schemaTree.tree.column("key", stretch=0, width=100)
        self.schemaTree.tree.configure(selectmode="browse")
        self.schemaTree.grid(row=1, column=0, padx=5, pady=3, sticky="NSEW")
        self.schemaFrm.columnconfigure(0, weight=1)
        self.schemaFrm.rowconfigure(1, weight=1)

        # Currently selected table
        self.curTabLab = tk.Label(self, anchor="nw", text="Selected Table:")
        self.curTabLab.grid(row=0, column=1, columnspan=1,
                            padx=5, pady=3, sticky="EW")
        self.curTab1Lab = tk.Label(self, anchor="nw",  background="white",
                                   relief="solid", borderwidth="1")
        self.curTab1Lab.grid(row=0, column=2, columnspan=2,
                            padx=5, pady=3, sticky="EW")
        sep = ttk.Separator(self, orient="horizontal")
        sep.grid(row=1, column=1, columnspan=3, padx=5, pady=5, sticky="NSEW")
        
        # View table
        self.titleViewLab = tk.Label(self, anchor="nw",
                                  text="View A Table in the Plot Window:")
        self.titleViewLab.grid(row=2, column=1, columnspan=3,
                            padx=5, pady=3, sticky="EW")
        self.nrows1Lab = tk.Label(self, anchor="nw", text="Number of Rows:")
        self.nrows1Lab.grid(row=3, column=1, padx=5, pady=3, sticky="EW")
        nRowLst = ["20", "50", "100", "500", "all"]
        self.nrows1Comb = ttk.Combobox(self, values=nRowLst, state="readonly")
        self.nrows1Comb.current(0)
        self.nrows1Comb.grid(row=3, column=2, padx=5, pady=3, sticky="EW")
        self.viewBtn = ttk.Button(self, width=10, text="Go", state="disabled",
                                  command=lambda action="view_table" :
                                  self._handlerGoButton(action))
        self.viewBtn.grid(row=3, column=3, rowspan=1, padx=5, pady=2,
                          sticky="NSEW")
        
        sep = ttk.Separator(self, orient="horizontal")
        sep.grid(row=4, column=1, columnspan=3,padx=5, pady=5, sticky="NSEW")
        
        # Export table(s)
        self.title3Lab = tk.Label(self, anchor="nw",
                                  text="Export a Table to Disk:")
        self.title3Lab.grid(row=5, column=1, padx=5, pady=3, sticky="EW")
        
        self.fmt1Lab = tk.Label(self, anchor="nw", text="Save Format:")        
        self.fmt1Lab.grid(row=6, column=1, padx=5, pady=3, sticky="EW")
        fmtLst = ["CSV"] #, "TSV", "VOT"]
        self.fmt1Comb = ttk.Combobox(self, values=fmtLst, state="readonly")
        self.fmt1Comb.current(0)
        self.fmt1Comb.grid(row=6, column=2, padx=5, pady=3, sticky="EW")
        self.direxp1Lab = tk.Label(self, justify="left", anchor="e",
                                  text="Export Directory:")
        self.direxp1Lab.grid(row=7, column=1, padx=5, pady=5, sticky="NW")
        self.browse1Btn = ttk.Button(self, text="Browse",
                                     command=lambda action="browse_export":
                                     self._handlerBrowseButton(action))
        self.browse1Btn.grid(row=7, column=3, padx=5, pady=2,sticky="NSEW" )
        self.exp1Dir = tk.StringVar()
        self.direxp1Ent = ttk.Entry(self, textvariable=self.exp1Dir)
        self.direxp1Ent.grid(row=8, column=1, columnspan=3, padx=5, pady=5,
                            sticky="EW")
        self.exp1Dir.set(os.getcwd())
        self.expBtn = ttk.Button(self, width=10, text="Go", state="disabled",
                                 command=lambda action="export_table":
                                 self._handlerGoButton(action))
        self.expBtn.grid(row=6, column=3, rowspan=1, padx=5, pady=2,
                         sticky="NSEW")
        
        # Blank row 9 is stretchable separator
        #sep = ttk.Separator(self, orient="horizontal")
        #sep.grid(row=9, column=0, columnspan=4, padx=5, pady=15,
        #sticky="NSEW")
        
        # SQL Query - TODO: IN DEVELOPMENT
        self.sqlFrm = tk.Frame(self)
        self.sqlFrm.grid(row=10, column=0, columnspan=6, padx=5, pady=3, 
                         sticky="NSEW")
        self.title4Lab = tk.Label(self.sqlFrm, anchor="nw",
                                  text="Run a SQL Query on the Database:")
        self.title4Lab.grid(row=0, column=0, columnspan=6, padx=5, pady=3,
                            sticky="EW")
        self.sqlTextBox = ScrolledText(self.sqlFrm, height=3)
        self.sqlTextBox.grid(row=1, column=0, columnspan=7, padx=5, pady=3, 
                             sticky="NSEW")
        self.sqlTextBox.frame.configure(borderwidth="0")
        sql = "SELECT * FROM sourceCat"
        self.sqlTextBox.insert("1.0", sql)
        self.queryBtn = ttk.Button(self.sqlFrm, text="Run Query", width=10,
                                   command=self._handlerQueryButton)
        self.queryBtn.grid(row=2, column=0, rowspan=3, padx=5, pady=5,
                           sticky="NSEW" )
        self.sqlOutChoice = tk.StringVar()
        self.sqlViewRad = ttk.Radiobutton(self.sqlFrm, 
                                      variable=self.sqlOutChoice, value="view",
                                 text="View result table in plotting window")
        self.sqlViewRad.grid(row=2, column=1, padx=5, pady=2, sticky="NW" )
        self.sqlDBRad = ttk.Radiobutton(self.sqlFrm, 
                                        variable=self.sqlOutChoice, value="db",
                                 text="Save result as new table in database")
        self.sqlDBRad.grid(row=3, column=1, padx=5, pady=2, sticky="NW" )
        self.sqlExpRad = ttk.Radiobutton(self.sqlFrm, 
                                       variable=self.sqlOutChoice, value="exp",
                                       text="Export result of query to disk")
        self.sqlExpRad.grid(row=4, column=1, padx=5, pady=2, sticky="NW")
        self.sqlOutChoice.set("view")

        sep = ttk.Separator(self.sqlFrm, orient="vertical")
        sep.grid(row=2, column=2, rowspan=3, padx=5, pady=5, sticky="NS")
        
        self.table3Lab = tk.Label(self.sqlFrm, anchor="nw",
                                  text="New Table:")
        self.table3Lab.grid(row=2, column=3, padx=5, pady=3, sticky="EW")
        self.table1Ent = ttk.Entry(self.sqlFrm)
        self.table1Ent.insert(0, "newTable")
        self.table1Ent.grid(row=2, column=4, padx=5, pady=5, sticky="EW")
        
        self.fmt2Lab = tk.Label(self.sqlFrm, anchor="nw", text="Save Format:")
        self.fmt2Lab.grid(row=2, column=5, padx=5, pady=3, sticky="EW")
        fmtLst = ["CSV"] #, "TSV", "VOT"]
        self.fmt2Comb = ttk.Combobox(self.sqlFrm, values=fmtLst, width=5,
                                     state="readonly")
        self.fmt2Comb.current(0)
        self.fmt2Comb.grid(row=2, column=6, padx=5, pady=3, sticky="EW")
        

        self.direxp2Lab = tk.Label(self.sqlFrm, justify="left", anchor="e",
                                  text="Export Directory:")
        self.direxp2Lab.grid(row=3, column=3, padx=5, pady=5, sticky="NW")
        self.browse2Btn = ttk.Button(self.sqlFrm, text="Browse",
                                     command=lambda action="browse_query":
                                     self._handlerBrowseButton(action))
        self.browse2Btn.grid(row=3, column=6, padx=5, pady=2,sticky="NSEW" )
        self.exp2Dir = tk.StringVar()
        self.direxp2Ent = ttk.Entry(self.sqlFrm, textvariable=self.exp2Dir)
        self.direxp2Ent.grid(row=4, column=3, columnspan=4, padx=5, pady=5,
                            sticky="EW")
        self.exp2Dir.set(os.getcwd())
        
        self.sqlFrm.rowconfigure(1, weight=1)
        self.sqlFrm.columnconfigure(2, weight=1)
        
        # Bindings
        self.schemaTree.bind("<<schema_table_selected>>",
                             self.on_tree_selected)

        # Set the expansion behaviour
        self.columnconfigure(0, weight=1)
        self.rowconfigure(10, weight=1)

    def _handlerGoButton(self, action):
        """Handle events for the view and export table buttons"""
        if self.tableName is None or self.tableName=="":
            errStr = "Please select a table in the box on the left."
            tkMessageBox.showinfo("Error", errStr)
            return
        
        if action=="view_table":
            self.sql = "SELECT * FROM %s" % self.tableName
            self.rowLimit =  self.nrows1Comb.get()
            try:
                limit = int(rowLimit)
                if limit>0:
                    self.sql += " LIMIT %d" % limit
            except Exception:
                pass
            self.event_generate("<<view_sql_table>>")
            
        if action=="export_table":
            self.sql = "SELECT * FROM %s" % self.tableName
            self.event_generate("<<export_sql_table>>")
        
    def _handlerBrowseButton(self, action):
        """Open the file selection dialog."""
        saveDir = filedialog.askdirectory(parent=self, initialdir=".",
             title="Please select a directory in which to save files")
        if saveDir!="":
            if action=="browse_export":
                self.exp1Dir.set(saveDir)
            if action=="browse_query":
                self.exp2Dir.set(saveDir)
            
    def _handlerQueryButton(self):
        """Run a custom SQL query on the data"""
        self.sql = self.sqlTextBox.get("1.0","end-1c")
        self.event_generate("<<run_custom_sql>>")

    def load(self, dataMan):
        """Load the tables from a datamanager instance"""
        self.sql = None
        self.dataMan = dataMan        
        schemaDict = self.dataMan.get_database_schema()
        self.clear_schematree()
        self.load_schematree(schemaDict)        
        self.viewBtn.configure(state="enabled")
        self.expBtn.configure(state="enabled")
        
    def on_tree_selected(self, event=None):        
        dummy, self.tableName = self.schemaTree.get_text_selected()
        self.curTab1Lab.configure(text=self.tableName)
        
    def clear_schematree(self):
        """Clear all the entries from the schema tree."""
        try:
            x = self.schemaTree.tree.get_children()
            for entry in x:
                self.schemaTree.tree.delete(entry)
        except Exception:
            pass

    def load_schematree(self, schemaDict):
        """Load the database description into the tree"""
        tabNameLst = schemaDict.keys()
        for tabName in tabNameLst:
            descArr = schemaDict[tabName]
            tabNode = self.schemaTree.tree.insert('', 'end', text=tabName)
            for e in descArr:
                colNode = self.schemaTree.tree.insert(tabNode, 'end',
                    text=e['name'],values=[e['name'], e["type"],
                                           "yes" if e["pk"]==1 else "no"])
    

#-----------------------------------------------------------------------------#
class PlotSctHstFrame(tk.Frame):
    """Frame presenting an interface to create a histogram or scatter plot
    driven by SQL queries. Up to four queries may be specified."""

    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        self.plotParm = None
    
        # Box in which to write the queries. 
        msg = "Enter up to 4 SQL queries on separate lines:"
        self.title1Lab = tk.Label(self, anchor="nw", text=msg)
        self.title1Lab.grid(row=0, column=0, columnspan=3, padx=5, pady=3,
                            sticky="EW")
        self.sqlTextBox = ScrolledText(self)
        #self.sqlTextBox.tag_configure("even", background="#e0e0e0")
        #self.sqlTextBox.tag_configure("odd", background="#ffffff")
        self.sqlTextBox.grid(row=1, column=0, columnspan=3, rowspan=7,
                             padx=5, pady=3, sticky="NSEW")
        self.sqlTextBox.frame.configure(borderwidth="0")

        # Buttons to create the plot
        self.title2Lab = tk.Label(self, justify="center", anchor="nw",
                                  text="Plotting Actions:")
        self.title2Lab.grid(row=0, column=3, columnspan=1, padx=5, pady=3,
                            sticky="W")
        sep = ttk.Separator(self, orient="horizontal")
        sep.grid(row=1, column=3, columnspan=1, padx=5, pady=5, sticky="NSEW")
        self.loadLab = tk.Label(self, text="Load parameters from file:")
        self.loadLab.grid(row=2, column=3, padx=5, pady=3, sticky="W")
        self.browseBtn = ttk.Button(self, text="Browse",
                                    command=self._handlerLoadQueryButton)
        self.browseBtn.grid(row=3, column=3, padx=5, pady=2, sticky="EW")

        self.saveLab = tk.Label(self, text="Save parameters to file:")
        self.saveLab.grid(row=4, column=3, padx=5, pady=3, sticky="W")
        self.saveBtn = ttk.Button(self, text="Browse",
                                  command=self._handlerSaveQueryButton)
        self.saveBtn.grid(row=5, column=3, padx=5, pady=2,sticky="EW")        
        self.plotBtn = ttk.Button(self, text="\nCREATE PLOT\n", width=10,
                                  command=self._handlerPlotButton)
        self.plotBtn.grid(row=7, column=3, padx=5, pady=2, sticky="SEW" )
        
        # Box for plot type and axis scaling
        self.axesFrm = ttk.LabelFrame(self, text=" Plot Type & Axis Scaling ")
        self.axesFrm.grid(column=0, row=8, padx=5, pady=2, sticky="NSEW")
        self.pltTypLab = tk.Label(self.axesFrm, anchor="nw", text="Plot Type:")
        self.pltTypLab.grid(row=0, column=0, padx=5, pady=3, sticky="E")
        self.pltTypLst = ["Histogram", "Scatter"]
        self.pltTypComb = ttk.Combobox(self.axesFrm, values=self.pltTypLst,
                                       width=10)
        self.pltTypComb.current(0)
        self.pltTypComb.grid(row=0, column=1, padx=5, pady=3, sticky="EW")
        self.nBinLab = ttk.Label(self.axesFrm, anchor="nw",
                                    text="Number of bins:")
        self.nBins = tk.StringVar()
        self.nBinLab.grid(row=1, column=0, padx=5, pady=3, sticky="E")
        self.nBinEnt = ttk.Entry(self.axesFrm, width=10)
        self.nBinEnt.grid(row=1, column=1,  padx=5, pady=3, sticky="EW")
        
        self.doLogX = tk.IntVar()
        self.logXChk = ttk.Checkbutton(self.axesFrm, text="Log X",
                                       variable=self.doLogX, onvalue=1,
                                       offvalue=0)
        self.logXChk.grid(row=2, column=0, columnspan=1, padx=5, pady=2,
                          sticky="E")
        self.doLogY = tk.IntVar()
        self.logYChk = ttk.Checkbutton(self.axesFrm, text="Log Y",
                                       variable=self.doLogY, onvalue=1,
                                       offvalue=0)
        self.logYChk.grid(row=2, column=1, columnspan=1, padx=5, pady=2,
                          sticky="E")
        self.zPowLab = tk.Label(self.axesFrm, anchor="nw", text="Z Power:")
        self.zPowLab.grid(row=3, column=0, padx=5, pady=3, sticky="E")
        zPowLst = ["0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9",
                   "1.0", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "1.8",
                   "1.9", "2.0"]
        self.zPower = tk.StringVar()
        self.zPowComb = ttk.Combobox(self.axesFrm, textvariable=self.zPower,
                                     values=zPowLst, width=5)
        self.zPowComb.current(8)
        self.zPowComb.grid(row=3, column=1, padx=5, pady=3, sticky="EW")
        
        # Axis labels
        self.titleFrm = ttk.LabelFrame(self, text="  Axis Labels  ")
        self.titleFrm.grid(column=1, row=8, padx=5, pady=2, sticky="NSEW")
        self.xLabelLab = tk.Label(self.titleFrm, justify="center", anchor="n",
                                  text="X-Axis:")
        self.xLabelLab.grid(row=0, column=0, padx=5, pady=3, sticky="E")
        self.xLabelEnt = ttk.Entry(self.titleFrm, width=20)
        self.xLabelEnt.grid(row=0, column=1,  padx=5, pady=3, sticky="EW")
        self.yLabelLab = tk.Label(self.titleFrm, justify="center", anchor="n",
                                  text="Y-Axis:")
        self.yLabelLab.grid(row=1, column=0, padx=5, pady=3, sticky="E")
        self.yLabelEnt = ttk.Entry(self.titleFrm, width=20)
        self.yLabelEnt.grid(row=1, column=1,  padx=5, pady=3, sticky="EW")     
        self.zLabelLab = tk.Label(self.titleFrm, justify="center", anchor="n",
                                  text="Z-Axis:")
        self.zLabelLab.grid(row=2, column=0, padx=5, pady=3, sticky="E")
        self.zLabelEnt = ttk.Entry(self.titleFrm, width=20)
        self.zLabelEnt.grid(row=2, column=1,  padx=5, pady=3, sticky="EW")
        self.tLabelLab = tk.Label(self.titleFrm, anchor="n",
                                  text="Title:")
        self.tLabelLab.grid(row=3, column=0, padx=5, pady=3, sticky="E")
        self.tLabelEnt = ttk.Entry(self.titleFrm, width=20)
        self.tLabelEnt.grid(row=3, column=1,  padx=5, pady=3, sticky="EW")
        self.titleFrm.columnconfigure(1, weight=1)
        
        # Legend labels
        self.legendFrm = ttk.LabelFrame(self, text="  Legend Labels  ")
        self.legendFrm.grid(column=2, row=8, padx=5, pady=2, sticky="NSEW")
        self.qEntLst = []
        self.q1Lab = tk.Label(self.legendFrm, justify="center", anchor="n",
                                  text="Query 1:")
        self.q1Lab.grid(row=0, column=0, padx=5, pady=3, sticky="W")
        self.qEntLst.append(ttk.Entry(self.legendFrm, width=20))
        self.qEntLst[-1].grid(row=0, column=1,  padx=5, pady=3, sticky="EW")
        self.q2Lab = tk.Label(self.legendFrm, justify="center", anchor="n",
                                  text="Query 2:")
        self.q2Lab.grid(row=1, column=0, padx=5, pady=3, sticky="W")
        self.qEntLst.append(ttk.Entry(self.legendFrm, width=20))
        self.qEntLst[-1].grid(row=1, column=1,  padx=5, pady=3, sticky="EW")
        self.q3Lab = tk.Label(self.legendFrm, justify="center", anchor="n",
                                  text="Query 3:")
        self.q3Lab.grid(row=2, column=0, padx=5, pady=3, sticky="W")
        self.qEntLst.append(ttk.Entry(self.legendFrm, width=20))
        self.qEntLst[-1].grid(row=2, column=1,  padx=5, pady=3, sticky="EW")
        self.q4Lab = tk.Label(self.legendFrm, justify="center", anchor="n",
                                  text="Query 4:")
        self.q4Lab.grid(row=3, column=0, padx=5, pady=3, sticky="W")
        self.qEntLst.append(ttk.Entry(self.legendFrm, width=20))
        self.qEntLst[-1].grid(row=3, column=1,  padx=5, pady=3, sticky="EW")
        self.legendFrm.columnconfigure(1, weight=1)
        
        # Data clipping, X, Y, Z
        self.clipFrm = ttk.LabelFrame(self, text="  Data Clipping  ")
        self.clipFrm.grid(column=3, row=8, padx=5, pady=2, sticky="NSEW")
        self.minLab = tk.Label(self.clipFrm, justify="center", anchor="n",
                                  text="Minimum")
        self.minLab.grid(row=0, column=1, padx=5, pady=3, sticky="EW")
        self.maxLab = tk.Label(self.clipFrm, justify="center", anchor="n",
                               text="Maximum")
        self.maxLab.grid(row=0, column=2, padx=5, pady=3, sticky="EW")
        self.clipXLab = tk.Label(self.clipFrm, justify="center", anchor="n",
                                  text="X:")
        self.clipXLab.grid(row=1, column=0, padx=5, pady=3, sticky="W")
        self.minXEnt = ttk.Entry(self.clipFrm, width=10)
        self.minXEnt.grid(row=1, column=1,  padx=5, pady=3)     
        self.maxXEnt = ttk.Entry(self.clipFrm, width=10)
        self.maxXEnt.grid(row=1, column=2,  padx=5, pady=3)
        self.clipYLab = tk.Label(self.clipFrm, justify="center", anchor="n",
                                  text="Y:")
        self.clipYLab.grid(row=2, column=0, padx=5, pady=3, sticky="W")
        self.minYEnt = ttk.Entry(self.clipFrm, width=10)
        self.minYEnt.grid(row=2, column=1,  padx=5, pady=3)     
        self.maxYEnt = ttk.Entry(self.clipFrm, width=10)
        self.maxYEnt.grid(row=2, column=2,  padx=5, pady=3)     
        self.clipZLab = tk.Label(self.clipFrm, justify="center", anchor="n",
                                  text="Z:")
        self.clipZLab.grid(row=3, column=0, padx=5, pady=3, sticky="W")
        self.minZEnt = ttk.Entry(self.clipFrm, width=10)
        self.minZEnt.grid(row=3, column=1,  padx=5, pady=3)     
        self.maxZEnt = ttk.Entry(self.clipFrm, width=10)
        self.maxZEnt.grid(row=3, column=2,  padx=5, pady=3)     

        # Set grid sections to stretch
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=1)
        self.rowconfigure(7, weight=1)

        # Set defaults
        self.reset_defaults()
        
    def _handlerLoadQueryButton(self):
        """Open the file selection dialog to load a query file."""
        queryFileName = filedialog.askopenfilename(parent=self,
                               title="Please select a parameter file to load:")
        if queryFileName!="":
            self.load_queryfile(queryFileName)
            
    def _handlerSaveQueryButton(self):
        """Open the file save dialogue to save a query file"""
        queryFileName = filedialog.asksaveasfilename(parent=self,
                           title="Save current plotting parameters to a file:")
        if queryFileName!="":
            self.save_queryfile(queryFileName)

    def _handlerPlotButton(self):
        self.event_generate("<<plot_query_scthst>>")

    def load_queryfile(self, queryFileName):
        """Load settings from a query file into the GUI"""
        self.plotParm = PlotParms(queryFileName)
        pDict = self.plotParm.configDict
        queryLst = self.plotParm.queryLst
        queryLabLst = self.plotParm.queryLabLst
        self.reset_defaults()
        
        # Query box
        #tag = "odd"
        for i in range(len(queryLst)):
            self.sqlTextBox.insert("1.0", queryLst[i] + "\n")#, tag)
            #tag = "even" if tag == "odd" else "odd"
            
        # Plot type
        valueLst = self.pltTypComb["values"]        
        pDict["TYPE"] = pDict.get("TYPE", valueLst[0])
        if not pDict["TYPE"] in valueLst:
            pDict["TYPE"] = valueLst[0]
        self.pltTypComb.current(valueLst.index(pDict["TYPE"]))
        
        # Number of bins
        nBins = pDict.get("NBINS", self.nBinEnt.get()) 
        self.nBinEnt.delete(0,tk.END)
        self.nBinEnt.insert(0, nBins)
        
        # Log X radiobutton
        pDict["DOLOGX"] = pDict.get("DOLOGX", "0")
        if pDict["DOLOGX"] in ["Y", "y", "1"]:
            self.doLogX.set(1)
        else:
            self.doLogX.set(0)
            
        # Log Y radiobutton
        pDict["DOLOGY"] = pDict.get("DOLOGY", "0")
        if pDict["DOLOGY"] in ["Y", "y", "1"]:
            self.doLogY.set(1)
        else:
            self.doLogY.set(0)
            
        # Z Power combobox (set nearest value)
        valueLst = [float(x) for x in self.zPowComb["values"]]
        cmdZpow = float(pDict.get("ZPOWER", "1.0"))
        closeZpow = min(valueLst, key=lambda x:abs(x-cmdZpow))
        self.zPowComb.current(valueLst.index(closeZpow))
        
        # Axis labels and title
        self.xLabelEnt.insert(0, pDict.get("XLABEL", ""))
        self.yLabelEnt.insert(0, pDict.get("YLABEL", ""))
        self.zLabelEnt.insert(0, pDict.get("ZLABEL", ""))
        self.tLabelEnt.insert(0, pDict.get("TITLE", ""))
        
        # Legend labels
        for i in range(len(queryLabLst)):
            self.qEntLst[i].delete(0,tk.END)
            self.qEntLst[i].insert(0, queryLabLst[i])
            
        # Data clipping
        self.minXEnt.insert(0, pDict.get("XDATAMIN", ""))
        self.maxXEnt.insert(0, pDict.get("XDATAMAX", ""))
        self.minYEnt.insert(0, pDict.get("YDATAMIN", ""))
        self.maxYEnt.insert(0, pDict.get("YDATAMAX", ""))
        self.minZEnt.insert(0, pDict.get("ZDATAMIN", ""))
        self.maxZEnt.insert(0, pDict.get("ZDATAMAX", ""))

    def save_queryfile(self, queryFileName):
        pass
            
    def reset_defaults(self):
        """Reset to default all plot settings"""
        self.sqlTextBox.delete("1.0",tk.END)
        self.pltTypComb.current(0)        
        self.nBinEnt.delete(0,tk.END)        
        self.nBinEnt.insert(0, "10")
        self.doLogX.set(0)
        self.doLogY.set(0)
        self.zPowComb.current(8)
        self.xLabelEnt.delete(0,tk.END)
        self.yLabelEnt.delete(0,tk.END)
        self.zLabelEnt.delete(0,tk.END)
        self.tLabelEnt.delete(0,tk.END)
        for i in range(len(self.qEntLst)):
            self.qEntLst[i].delete(0,tk.END)
            self.qEntLst[i].insert(0, "Query %s" % (i+1))
        self.minXEnt.delete(0,tk.END)
        self.maxXEnt.delete(0,tk.END)
        self.minYEnt.delete(0,tk.END)
        self.maxYEnt.delete(0,tk.END)
        self.minZEnt.delete(0,tk.END)
        self.maxZEnt.delete(0,tk.END)

    def update_plotparm_obj(self):
        """Push the plot settings in the GUI ro the PlotParm object."""
        textStr = self.sqlTextBox.get("1.0","end")
        textStr = cleanup_str_input(str(textStr))
        self.plotParm.queryLst = textStr.split("\n")
        labLst = []
        for i in range(len(self.qEntLst)):
            labLst.append(self.qEntLst[i].get())
        self.plotParm.queryLabLst = labLst
        self.plotParm.configDict["TYPE"] = self.pltTypComb.get()
        self.plotParm.configDict["NBINS"] = self.nBinEnt.get()
        self.plotParm.configDict["DOLOGX"] = self.doLogX.get()
        self.plotParm.configDict["DOLOGY"] = self.doLogY.get()
        self.plotParm.configDict["ZPOWER"] = self.zPower.get()
        self.plotParm.configDict["XDATAMIN"] = self.minXEnt.get()
        self.plotParm.configDict["XDATAMAX"] = self.maxXEnt.get()
        self.plotParm.configDict["YDATAMIN"] = self.minYEnt.get()
        self.plotParm.configDict["YDATAMAX"] = self.maxYEnt.get()
        self.plotParm.configDict["ZDATAMIN"] = self.minZEnt.get()
        self.plotParm.configDict["ZDATAMAX"] = self.maxZEnt.get()
        self.plotParm.configDict["TITLE"] = self.tLabelEnt.get()
        self.plotParm.configDict["XLABEL"] = self.xLabelEnt.get()
        self.plotParm.configDict["YLABEL"] = self.yLabelEnt.get()
        self.plotParm.configDict["ZLABEL"] = self.zLabelEnt.get()

        
#-----------------------------------------------------------------------------#
class SingleFigFrame(tk.Frame):
    
    def __init__(self, parent, fig, titleStr=""):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        self.frame = tk.Frame(self.parent)
        self.frame.grid(row=0, column=0, padx=5, pady=5, sticky="NSEW")
        
    
        # Layout the figure, title and button row 
        self.titleLab = tk.Label(self.frame, justify="left", anchor="nw",
                                 text=titleStr)
        self.titleLab.grid(row=0, column=0, padx=5, pady=3, sticky="NW")
        figCanvas = FigureCanvasTkAgg(fig, master=self.frame)
        figCanvas.show()
        self.myCan = figCanvas.get_tk_widget()
        self.myCan.grid(row=1, column=0, padx=5, pady=5, sticky="NSEW")
        self.toolbarFrame = tk.Frame(self.frame)
        toolbar = NavigationToolbar2TkAgg(figCanvas, self.toolbarFrame)
        toolbar.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.toolbarFrame.grid(row=2, column=0, padx=5, pady=5, sticky="EW")
        self.frame.columnconfigure(0, weight=1)
        self.frame.rowconfigure(1, weight=1)

        
#-----------------------------------------------------------------------------#
class SingleTabFrame(tk.Frame):
    """Frame presenting a single table from a recArray"""

    def __init__(self, parent, title="", footer=""):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        
        # Table 
        self.titleLab = tk.Label(self, justify="left", anchor="nw",
                                  text=title)
        self.titleLab.grid(row=0, column=0, padx=5, pady=3, sticky="EW")
        self.table = ScrolledTreeTab(self)
        self.table.grid(column=0, row=1, padx=5, pady=3, sticky="NSEW")
        self.footerLab = tk.Label(self, justify="left", anchor="nw",
                                  text=footer)
        self.footerLab.grid(column=0, row=2,padx=5, pady=3, sticky="EW")
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        

#-----------------------------------------------------------------------------#
if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    

