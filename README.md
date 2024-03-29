# Untargeted Region of Interest Selection for GC-MS Data using a Pseudo F-Ratio Moving Window (ψFRMV)

## Ryland T. Giebelhaus, Michael D. Sorochan Armstrong, A. Paulina de la Mata, and James J. Harynuk*

### University of Alberta, Edmonton, Canada
#### james.harynuk@ualberta.ca

## 1.0 About
This is a MatLab algorithm developed to locate regions of interest (ROI) within 1D GC-MS datasets. The intention is that the ROIs are selected then subjected to further chemometrics, especially approaches where a low number of components is desired.

The current release (v1.3) is the version described in the paper we recently published in the Journal of Chromatography A. See Referencing section below for the paper and how to cite this work.

## 2.0 Use
Users input data as an **M** by **N** matrix into the function, where **M** is individual scans and **N** is the _m/z_ of each ion, and the matrix is ion intensities. 

Download our latest release (v1.3), open in the most recent release of MatLab, and enjoy. 

### 2.1 Inputs
* **data**: M x N array of ion intensities, where M is scans and N is ion m/z.
* **wndw**: moving window size. Ideally this should be about the approximate width of the average peak in your separation. 10 is a good starting point if unsure.
* **CutOff**: probability cutoff. 0.7 is a good place to start.

### 2.2 Outputs
* **pv**: probability values for each scan. Is overlaid on output plot by default.
* **modPVans**: modified probability values, above the CutOff threshold.
* **ticData**: Total Ion Chromatogram of the data. Plotted by default.
* **noiseDroppedTIC**: TIC with non ROIs dropped.
* **noiseDropped**: all ions and scans with noise (non ROIs) dropped. This is intended for use in further analysis.
* **boolCutOff**: 0 = scan not in ROI, 1 = scan in ROI.

## 3.0 Referencing
Please cite the following paper when using this work.

* Ryland T. Giebelhaus, Michael D. Sorochan Armstrong, A. Paulina de la Mata, and James J. Harynuk. Untargeted region of interest selection for gas chromatography - mass spectrometry data using a pseudo F-ratio moving window, Journal of Chromatography A, 2022, 1682, 463499; 10.1016/j.chroma.2022.463499
