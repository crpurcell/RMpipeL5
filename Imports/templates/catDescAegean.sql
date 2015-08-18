#-----------------------------------------------------------------------------#
# SQL Description of the input ASCII catalogue.
#-----------------------------------------------------------------------------#
# The table must be called sourceCat and have at least two columns:
#    x_deg double         # World coordinate X-position
#    y_deg double         # World coordinate Y-position
#
# A column named 'uniqueName' will be created by the pipeline script to serve
# as a primary key, so this column name is reserved. Any additional columns
# will be parsed and inserted into the database. All columns in the ASCII file
# MUST be described. Valid types are 'double', 'varchar(n)', and 'int'.
#
# This sample file describes the format of the catalogue produced by the 
# Aegean source finder (http://adsabs.harvard.edu/abs/2012MNRAS.422.1812H)
# https://www.sites.google.com/site/mrpaulhancock/data-and-code/aegean
#
#-----------------------------------------------------------------------------#

CREATE TABLE sourceCat (
inName varchar(20),
bkg double,
rms double,
RAstr varchar(20),
DECstr varchar(20),
x_deg double,
dX_deg double,
y_deg double,
dY_deg double,
I_jybm double,
dI_jybm double,
Sint_jy double,
dSint_jy double,
majAxis_asec double,
dMajAxis_asec double,
minAxis_asec double,
dMinAxis_asec double,
PA_deg double,
dPA_deg double,
flags varchar(100)
);
