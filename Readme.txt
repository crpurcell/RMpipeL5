# ABOUT:
# This is the prototype RM-synthesis pipeline intended for use on Level-5 
# ASKAP data. See POSSUM Report 62 for a full description.


#-----------------------------------------------------------------------------#
# USAGE INSTRUCTIONS:
# Execute the following Python scripts in order to create a test dataset and
# run the Level-5 RM-pipeline on the data.

./0_mk_test_image_data.py testData/
./1_verify_image_data.py testData/
./2_create_image_session.py -o testData/ testSession/ testData/testCat.dat testData/testCatDesc.sql

#> Edit the file 'testSession/inputs.config' to modify default pipeline inputs.

./3_extract_spectra.py testSession/
./4_do_RM-synthesis.py testSession/
./5_do_RM-clean.py testSession/
./rmPipeViewer.py

# Note: a '-h' argument after most scripts will print help & usage information.

#-----------------------------------------------------------------------------#
# TODO PPC:

# Extractor
* Re-write the spectral extraction module to save a sparse FITS file of
  spectra containing valid emission. This can then be used with the
  other PPC schemes.
* Also allow the extraction module create a mini-cube on disk with
  area equivalent to the noise box.
* Set flag when source is near the edge of image or box contains NaNs.

# Pipeline RM and measurements
* Set flag when RM is detected near edge of spectrum.
* Create a best-fit RM-thin model and subtract to get residual.
* Implement complexity measurements based on residuals.
* Show clean cutoff and complexity thresholds on plots.


# GUI & interface
* Re-write pipeline inputs tab to show driving file on left and
  derived parmaters on right.
* Show the result values in the plotting window.
* Add a name filter to the table
* Annotate main result values on plots.
* Subclass the plot control bar and add a button to hide the legends
* Pipeline script to batch create publication plots and export results tables.

# LONG TERM TODO:
* Measure noise from an annulus instead of a box. Even using MADFM the
Stokes I RMS is biased.
* Implement other methods of source extraction in PPC vote.
* Turn on -o flag in RM-synthesis script?
* re-write feedback for RM-clean and implement useful logging.

#-----------------------------------------------------------------------------#

TODO PVACAT:

* New scripts to run qu-fitting and model comparison
