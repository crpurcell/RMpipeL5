#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_tk.py                                                        #
#                                                                             #
# PURPOSE:  Functions and classes for TK graphical elements.                  #
#                                                                             #
# MODIFIED: 29-May-2015 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
#  ScrolledTreeTab     ... scrolled Treeview widget to mimicing a table       #
#  ScrolledTreeView    ... standard scrolled TreeView widget                  #
#  ScrolledCanvasFrame ... frame embedded in a scrolling canvas               #
#                                                                             #
#=============================================================================#

import Tkinter as tk
import ttk
import tkFont


#-----------------------------------------------------------------------------#
class ScrolledTreeTab(tk.Frame):
    """
    Use a ttk.Treeview as a multicolumn ListBox with scrollbars.
    """
    
    def __init__(self, parent, virtEvent="<<tab_row_selected>>", strPad=20,
                 *args, **kw):
        tk.Frame.__init__(self, parent, *args, **kw)
        self.parent = parent
        self.rowSelected = None
        self.virtEvent = virtEvent
        self.strPad = strPad
        
        # Create the treeview and the scrollbars
        self.tree = ttk.Treeview(self, show="headings")
        vsb = ttk.Scrollbar(self, orient="vertical",
                            command=self.tree.yview)
        hsb = ttk.Scrollbar(self, orient="horizontal",
                            command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        # Grid the tree and scrollbars within the container frame
        self.tree.grid(column=0, row=0, sticky="NSWE")
        vsb.grid(column=1, row=0, sticky="NS")
        hsb.grid(column=0, row=1, sticky="WE")
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=0)
        self.rowconfigure(0, weight=1)
        
        # Limit to single row selections
        self.tree.configure(selectmode="browse")
        self.tree.bind("<ButtonRelease-1>", self.on_row_select)
        
    def insert_rows(self, rows, colNames=None):
        """Insert rows from a 2D iterable object."""

        # Populate the headers and set the sort function
        if colNames is None:
            colNames = ["Row "+ str(x+1) for x in range(len(rows[0]))]
        self.tree['columns'] = colNames
        for col in colNames:
            self.tree.heading(col, text=col, command=lambda c=col: \
                                  self._sortby(self.tree, c, 0))
        
            # Set the column width to the width of the header string
            strWidth = tkFont.Font().measure(col)
            self.tree.column(col, width=strWidth + self.strPad)
            self.tree.column(col, minwidth=strWidth + self.strPad)
            
        # Populate the rows
        rowIndx = 0
        for row in rows:
            row = [str(x) for x in row]   # Convert to plain strings
            self.tree.insert('', 'end', values=row, text=str(rowIndx))
            rowIndx += 1
            
            # Adjust the column width (& minwidth) to fit each value
            for i, val in enumerate(row):
                strWidth = tkFont.Font().measure(val)
                if self.tree.column(colNames[i], width=None)<\
                       (strWidth + self.strPad):
                    self.tree.column(colNames[i], width=strWidth +
                                     self.strPad)
                    self.tree.column(colNames[i], minwidth=strWidth +
                                     self.strPad)
        
    def on_row_select(self, event=None):
        """Store the index of the selected row and generate an event.
        The original index of the column is stored in the 'text' property
        of the TreeView widget"""
        item =  event.widget.identify("item", event.x, event.y)
        if not item=="":
            indx = event.widget.item(item, "text")
            self.rowSelected = int(indx)
            self.event_generate(self.virtEvent)

    def get_indx_selected(self):
        """Return the index of the last row selected."""        
        if self.rowSelected is None:
            return None
        else:
            return int(self.rowSelected)
        
    def insert_recarray(self, arr):
        """Insert a numpy.recarray into the treeview widget"""
        colNames = arr.dtype.names
        self.insert_rows(arr, colNames)

    def clear_entries(self):
        """Clear all the entries from the table."""
        try:
            x = self.tree.get_children() 
            for entry in x:
                self.tree.delete(entry)
        except Exception:
            pass
            
    def _sortby(self, tree, col, descending):
        """Sort tree contents when a column header is clicked."""

        # Fetch column IDs and data values to sort on
        data = [(tree.set(child, col), child)
                for child in tree.get_children('')]

        # If the data to be sorted is numeric change to float
        data = self._change_numeric_onestep(data)
        
        # now sort the data in place
        data.sort(reverse=descending)
        for i, item in enumerate(data):
            tree.move(item[1], '', i)
            
        # Toggle the sort function
        tree.heading(col,
                     command=lambda col=col: \
                     self._sortby(tree, col, int(not descending)))
        
    def _change_numeric_onestep(self, data):
        """If the data to be sorted is numeric change to float."""
        newData = []
        try:
            for child, col in data:
                if child=="None":
                    child = "-inf"   # Regard text "None" as -infinity
                newData.append((float(child), col))
            print "Sorting column as numeric data."            
            return newData
        except Exception:
            print "Sorting column as ascii data."
            return data


#-----------------------------------------------------------------------------#
class ScrolledTreeView(tk.Frame):

    def __init__(self, parent, virtEvent="<<tree_selected>>", *args, **kw):
        tk.Frame.__init__(self, parent, *args, **kw)
        self.parent = parent
        self.textSelected = None
        self.textRootSelected = None
        self.virtEvent = virtEvent
    
        # Create the treeview and the scrollbars
        self.tree = ttk.Treeview(self)
        vsb = ttk.Scrollbar(self, orient="vertical",
                            command=self.tree.yview)
        hsb = ttk.Scrollbar(self, orient="horizontal",
                            command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        
        # Grid the tree and scrollbars within the container frame
        self.tree.grid(column=0, row=0, sticky="NSEW")
        vsb.grid(column=1, row=0, sticky="NS")
        hsb.grid(column=0, row=1, sticky="WE")
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=0)
        self.rowconfigure(0, weight=1)

        # Binding to handle row selections
        self.tree.bind("<ButtonRelease-1>", self.on_row_select)
        
    def on_row_select(self, event=None):
        item =  event.widget.identify("item", event.x, event.y)
        if not item=="":
            self.textSelected = event.widget.item(item, "text")
            self.textRootSelected = event.widget.item(item, "text")
            while True:                
                parentItem =  event.widget.parent(item)
                if parentItem=="":
                    break
                self.textRootSelected = event.widget.item(parentItem, "text")
                item = parentItem
            self.event_generate(self.virtEvent)

    def get_text_selected(self):        
        return self.textSelected, self.textRootSelected


#-----------------------------------------------------------------------------#
class ScrolledCanvasFrame(tk.Frame):
    
    def __init__(self, parent, *args, **kw):
        tk.Frame.__init__(self, parent, bg="blue", *args, **kw)
        self.parent = parent

        # Create the canvas and the scrollbars
        self.canvas = tk.Canvas(self, bg="white", border=0)
        vsb = ttk.Scrollbar(self, orient="vertical",
                            command=self.canvas.yview)
        hsb = ttk.Scrollbar(self, orient="horizontal",
                            command=self.canvas.xview)
        self.canvas.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        
        self.canvas.xview_moveto(0)
        self.canvas.yview_moveto(0)
        
        # Grid the canvas and scrollbars within the container frame
        self.canvas.grid(column=0, row=0, sticky="NSEW")
        vsb.grid(column=1, row=0, sticky="NS")
        hsb.grid(column=0, row=1, sticky="WE")
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=0)
        self.rowconfigure(0, weight=1)
        
        # Now create a frame within the canvas
        self.interior = tk.Frame(self.canvas, bg="red")
        self.winID = self.canvas.create_window((0,0), window=self.interior,
                                         anchor="nw", tags="self.interior")
        
        self.interior.bind('<Configure>', self._configure_interior)
        self.canvas.bind('<Configure>', self._configure_canvas)

    def _configure_interior(self, event):
        size = (self.interior.winfo_reqwidth(),
                self.interior.winfo_reqheight())
        self.canvas.config(scrollregion="0 0 %s %s" % size)
        #self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        if self.interior.winfo_reqwidth() != self.canvas.winfo_width():
            self.canvas.config(width=self.interior.winfo_reqwidth())

    def _configure_canvas(self, event):
        if self.interior.winfo_reqwidth() != self.canvas.winfo_width():
            self.canvas.itemconfigure(self.winID,
                                      width=self.canvas.winfo_width())
        
