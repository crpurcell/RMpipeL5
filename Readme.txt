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
# TODO:
* allow user to specify a uniqueName instead of creating a
position-dependent one. PRIORITY
* Set flag when source is near the edge of image or contains NaNs.
* Set flag when RM is detected near edge of spectrum.
* Implement complexity measurements.
* Update the spectral plots with best fits for RM-synth. IN PROG
* Re-write & simplify the GUI to view plots and export all tables. IN PROG
* Add a filter to the table
* Annotate main results on plots

# LONG TERM TODO:
* Measure noise from an annulus instead of a box. Even using MADFM the
Stokes I RMS is biased.
* Implement other methods of source extraction in PPC vote.
* Turn on -o flag in RM-synthesis script? 
* re-write feedback for RM-clean and implement useful logging 
