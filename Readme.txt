# ABOUT:
# This is the prototype RM-synthesis pipeline intended for use on Level-5 
# ASKAP data. See POSSUM Report 62 for a full description.


#-----------------------------------------------------------------------------#
# USAGE INSTRUCTIONS:
# Execute the following Python scripts in order to create a test dataset and
# run the Level-5 RM-pipeline on the data.


./0_mk_test_ascii_data.py -n PAF_MKII_Tsys.dat -f 1.0e9,1.3e9 catalogue.csv testASCIIData/
./1_verify_ascii_data.py testASCIIData/
./2_create_session.py -o testSessionASCII/ testASCIIData/ testASCIIData/testCat.txt testASCIIData/testCatDesc.sql
#> Edit the file 'testSessionASCII/inputs.config' to modify default pipeline inputs.
./3_extract_spectra.py testSessionASCII/
./4_do_RM-synthesis.py testSessionASCII/
./5_do_RM-clean.py testSessionASCII/
./6_measure_complexity.py testSessionASCII/
./rmPipeViewer.py

# OR:

./0_mk_test_image_data.py -n PAF_MKII_Tsys.dat -f 1.0e9,1.3e9  catalogue.csv testImageData/
./1_verify_image_data.py testImageData/
./2_create_session.py -o testSessionImage/ testImageData/ testImageData/testCat.txt testImageData/testCatDesc.sql
#> Edit the file 'testSessionImage/inputs.config' to modify default pipeline inputs.
./3_extract_spectra.py testSessionImage/
./4_do_RM-synthesis.py testSessionImage/
./5_do_RM-clean.py testSessionImage/
./6_measure_complexity.py testSessionImage/
./rmPipeViewer.py

# Note: a '-h' argument after most scripts will print help & usage information.
x
#-----------------------------------------------------------------------------#
# TODO PPC:

# Extractor
* Set flag when source is near the edge of image or box contains NaNs.

# Pipeline RM and measurements
* Set flag when RM is detected near edge of spectrum.
* Create a best-fit RM-thin model and subtract to get residual.
TO BE PLOTTED
* Implement complexity measurements based on residuals. 
IN PROGRESS - IN TESTING


# GUI & interface
* Re-write pipeline inputs tab to show driving file on left and
  derived parmaters on right.
* Show the result values in the plotting window.
* (Add a name filter to the table)
* Annotate main result values on plots.
* Subclass the plot control bar and add a button to hide the legends
* Pipeline script to batch create publication plots and export results tables.

# LONG TERM TODO:
* re-write feedback for RM-clean and implement useful logging.

#-----------------------------------------------------------------------------#

TODO PVACAT:

* Integrate new scripts to run qu-fitting and model comparison
