
CREATE TABLE spectraParms (
uniqueName varchar(20) PRIMARY KEY,
boxScale_pix int(10),
fluxMedI_Jybm double,
fluxMedQ_Jybm double,
fluxMedU_Jybm double,
rmsMedI_Jybm double,
rmsMedQ_Jybm double,
rmsMedU_Jybm double,
rmsMedQUAvg_Jybm double,
nFreqChan int(10),
deltaFreqChan_Hz double,
coeffPolyIspec varchar(500),
bmaj_deg double,
bmin_deg double,
bpa_deg double,
pixscale_deg double,
fNormSumbox double,
fitIredChiSq double,
fitIstatus int(2),
extractStatus int(2)
);

CREATE TABLE dirtyFDFparms (
uniqueName varchar(20) PRIMARY KEY,
lam0Sq_m2 double,
freq0_Hz double,
nPhiChan int(10),
deltaPhiChan_rm2 double,
phiCentre_rm2 double,
weightType varchar(20),
fwhmRMSF double,
phiPeakPIchan_rm2 double,
dPhiPeakPIchan_rm2 double,
ampPeakPIchan_Jybm double,
ampPeakPIchanEff_Jybm double,
dAmpPeakPIchan_Jybm double,
snrPIchan double,
indxPeakPIchan double,
peakFDFimagChan double,
peakFDFrealChan double,
polAngleChan_deg double,
dPolAngleChan_deg double,
polAngle0Chan_deg double,
dPolAngle0Chan_deg double,
phiPeakPIfit_rm2 double,
dPhiPeakPIfit_rm2 double,
ampPeakPIfit_Jybm double,
ampPeakPIfitEff_Jybm double,
dAmpPeakPIfit_Jybm double,
snrPIfit double,
indxPeakPIfit double,
peakFDFimagFit double,
peakFDFrealFit double,
polAngleFit_deg double,
dPolAngleFit_deg double,
polAngle0Fit_deg double,
dPolAngle0Fit_deg double,
thresholdSignalPI double,
detectionF int(2),
status int(2)
);

CREATE TABLE cleanFDFparms (
uniqueName varchar(20) PRIMARY KEY,
phiPeakPIchan_rm2 double,
dPhiPeakPIchan_rm2 double,
ampPeakPIchan_Jybm double,
ampPeakPIchanEff_Jybm double,
dAmpPeakPIchan_Jybm double,
snrPIchan double,
indxPeakPIchan double,
peakFDFimagChan double,
peakFDFrealChan double,
polAngleChan_deg double,
dPolAngleChan_deg double,
polAngle0Chan_deg double,
dPolAngle0Chan_deg double,
phiPeakPIfit_rm2 double,
dPhiPeakPIfit_rm2 double,
ampPeakPIfit_Jybm double,
ampPeakPIfitEff_Jybm double,
dAmpPeakPIfit_Jybm double,
snrPIfit double,
indxPeakPIfit double,
peakFDFimagFit double,
peakFDFrealFit double,
polAngleFit_deg double,
dPolAngleFit_deg double,
polAngle0Fit_deg double,
dPolAngle0Fit_deg double,
thresholdSignalPI double,
detectionF int(2),
nIterDone int(100),
cleanCutoff_sigma double,
cleanCutoff_Jybm double,
status int(2)
);

CREATE TABLE complexMeasures (
uniqueName varchar(20) PRIMARY KEY,
complexM1 double,
complexM2 double,
complexM3 double,
status int(2)
);
