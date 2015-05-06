./1_verify_image_data.py DATA/C628/
./2_create_image_session.py DATA/C628/ SESSION/testSession catAegean.txt catAegeanDesc.sql -o
./3_extract_spectra.py SESSION/testSession/
./4_do_RM-synthesis.py SESSION/testSession/
./5_do_RM-clean.py SESSION/testSession/
./rmPipeViewer.py
