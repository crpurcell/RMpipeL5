# Database file path
DBFILE = atlas.sqlite

# Type of plot
TYPE = Histogram

# Queries to execute (<=4)
QUERY = SELECT phiPeakPIfit_rm2 FROM dirtyFDFparms WHERE phiPeakPIfit_rm2 IS NOT NULL
QUERY = SELECT phiPeakPIfit_rm2 FROM dirtyFDFparms WHERE phiPeakPIfit_rm2 IS NOT NULL
QUERY = SELECT phiPeakPIfit_rm2 FROM dirtyFDFparms WHERE phiPeakPIfit_rm2 IS NOT NULL
QUERY = SELECT phiPeakPIfit_rm2 FROM dirtyFDFparms WHERE phiPeakPIfit_rm2 IS NOT NULL

# Query labels (shown as legends on plot)
QLABEL = All

DOLOGX = n
DOLOGY = n
#ZPOWER = 1.1
# Data cutoffs to exclude outliers
#XDATAMIN = 2.5
#XDATAMAX = 10

# Plotting controls
NBINS = 20

# Axis Labels
TITLE = Faraday Depth of Dirty Peak
YLABEL = Number of Sources
XLABEL = Faraday Depth (rad/m^2)
ZLABEL = 