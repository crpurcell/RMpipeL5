
./0_mk_test_data.py
./1_verify_image_data.py testData/
./2_create_image_session.py -o testData/ testSession/ testData/testCat.dat testData/testCatDesc.sql
#> Edit the file 'testSession/inputs.config' (optional)
./3_extract_spectra.py testSession/
./4_do_RM-synthesis.py testSession/
./5_do_RM-clean.py testSession/
./rmPipeViewer.py
