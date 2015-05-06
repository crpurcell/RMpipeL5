#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_tk.py                                                        #
#                                                                             #
# PURPOSE:  Functions and classes for TK graphical elements.                  #
#                                                                             #
# MODIFIED: 21-Apr-2015 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
#  ScrolledTreeTab          ... scrolled Treeview widget to mimic a table     #
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
    
    def __init__(self, parent):
        tk.Frame.__init__(self, parent=None)
        self.parent = parent
        self.frame = tk.Frame(self.parent)
        
        # Create the treeview and the scrollbars
        self.tree = ttk.Treeview(self.frame, show="headings")
        vsb = ttk.Scrollbar(self.frame, orient="vertical",
                            command=self.tree.yview)
        hsb = ttk.Scrollbar(self.frame, orient="horizontal",
                            command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        # Grid the tree and scrollbars within the container frame
        self.tree.grid(column=0, row=0, sticky="NSWE")
        vsb.grid(column=1, row=0, sticky="NS")
        hsb.grid(column=0, row=1, sticky="WE")
        self.frame.columnconfigure(0, weight=1)
        self.frame.columnconfigure(1, weight=0)
        self.frame.rowconfigure(0, weight=1)
        
        # Limit to single row selections
        self.tree.configure(selectmode="browse")
        
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
            self.tree.column(col, width=strWidth)
            self.tree.column(col, minwidth=strWidth)
            
        # Populate the rows
        for row in rows:
            row = [str(x) for x in row]   # Convert to plain strings
            self.tree.insert('', 'end', values=row)
            
            # Adjust the column width (& minwidth) to fit each value
            for i, val in enumerate(row):
                strWidth = tkFont.Font().measure(val)
                if self.tree.column(colNames[i], width=None)<strWidth:
                    self.tree.column(colNames[i], width=strWidth)
                    self.tree.column(colNames[i], minwidth=strWidth)
        
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
        data = self._change_numeric(data)
        
        # now sort the data in place
        data.sort(reverse=descending)
        for i, item in enumerate(data):
            tree.move(item[1], '', i)
            
        # Toggle the sort function
        tree.heading(col,
                     command=lambda col=col: \
                     self._sortby(tree, col, int(not descending)))
        
    def _change_numeric(self, data):
        """If the data to be sorted is numeric change to float."""
        newData = []
        try:
            for child, col in data:            
                newData.append((float(child), col))
            return newData
        except Exception:            
            return data



#-----------------------------------------------------------------------------#
class TableFrame(tk.Frame):
    """Frame presenting the a table to the user."""
    
    def __init__(self, parent, title="Table1:", 
                 virtEvent="<<tab_row_selected>>",
                 foot="Click on a column header to sort up or down."):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        self.rowSelected = None
        self.virtEvent = virtEvent
        
        # Layout title, table & footnotes
        self.labTitle = tk.Label(self, justify="left", anchor="nw",
                                 text=title, font=("Helvatica", 10))
        self.labTitle.grid(row=0, column=0,  padx=5, pady=3, sticky="NW")
        self.table = ScrolledTreeTab(self)
        self.table.frame.grid(row=1, column=0, padx=5, pady=5, sticky="NSEW")
        self.labNote = tk.Label(self, justify="left", anchor="sw", text=foot)
        self.labNote.grid(row=2, column=0,  padx=5, pady=3, sticky="NW")
        
        # Allow only the table to expand
        self.rowconfigure(1, weight=1)
        self.columnconfigure(0, weight=1)
        
        # Binding to handle row selections
        self.table.tree.bind("<ButtonRelease-1>", self.on_row_select)
        
    def on_row_select(self, event=None):
        """Store the index of the selected row and generate an event."""
        
        item =  event.widget.identify("item", event.x, event.y)
        rowID = event.widget.identify_row(event.y)
        #rowIDdecimal =  int(rowID[1:], 16)
        indx = event.widget.index(rowID)
        if not item=="":
            self.rowSelected = indx
            self.event_generate(self.virtEvent)

    def get_indx_selected(self):
        """Return the index of the last row selected."""
        
        if self.rowSelected is None:
            return None
        else:
            return int(self.rowSelected)
            
    def insert_table(self, arr):
        """Insert a numpy.recarray into the treetable widget"""

        colNames = arr.dtype.names
        self.table.insert_rows(arr, colNames)

    def clear_entries(self):
        """Clear all entries from the widget"""
        
        self.table.clear_entries()
