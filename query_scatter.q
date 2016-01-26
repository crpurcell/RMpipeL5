# Database file path
DBFILE = atlas.sqlite

# Type of plot
TYPE = Scatter

# Queries to execute (<=4)
QUERY = SELECT spectraParms.fluxMedI_Jybm, dirtyFDFparms.ampPeakPIfit_Jybm, phiPeakPIfit_rm2 FROM spectraParms INNER JOIN dirtyFDFparms ON spectraParms.uniqueName = dirtyFDFparms.uniqueName

# Query labels (shown as legends on plot)
QLABEL = All

# Plotting controls
DOLOGX = n
DOLOGY = n
ZPOWER = 1.0

# Data cutoffs to exclude outliers
#XDATAMIN = 10
#XDATAMAX = 100
#YDATAMIN = .01
#YDATAMAX = 10
#ZDATAMIN = 0.0
#ZDATAMAX = 1

# Axis Labels
XLABEL = Median Flux (Jy/beam)
YLABEL = Polarised Intensity (Jy/beam)
ZLABEL = Faraday Depth of Peak (rad/m^2)
