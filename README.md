# Materials for "Improving Infrared Radiance Ensemble Data Assimilation though Mitigating Deleterious Non-Gaussian Artifacts"
Author: Man-Yau Chan, Department of Geography, The Ohio State University, Columbus, Ohio, United States of America

&nbsp; &nbsp;



## Abstract

Ensemble Kalman Filters (EnKFs) with quality control (QC) are frequently used to assimilate satellite-sensed infrared radiance (IR) observations. Unfortunately, the forecast ensembles used in such data assimilation (DA) often possess non-Gaussian statistics that violate the EnKFs' and QCs' assumption of Gaussian forecast statistics. These violations likely generate statistical artifacts in the analysis ensemble. Such artifacts may limit the beneficial impacts of EnKFs-based IR DA. However, it is not clear what those artifacts are and how they limit the performance of EnKFs-based IR BT DA. This study addresses that knowledge gap using perfect model observing system simulation experiments (OSSEs) of EnKFs-based IR DA. Two artefacts resulting from those assumption violations are identified: 1) the creation of clouds in subsaturated humidity conditions and 2) a general depletion of specific humidity. The first artifact results from the EnKFs' Gaussian forecast assumption and the latter arises from applying QC on left-skewed IR innovation distributions. As is demonstrated here, these artifacts cause severe cold temperature biases that annihilate the beneficial impacts of IR radiance DA. This study treats those artifacts by 1) disabling the hydrometeor analysis increments and 2) reducing the impacts of quality control. OSSEs demonstrate that these treatments noticeably improve the performance of IR DA. The findings from this study not only advances the community's understanding of how non-Gaussian statistics degrade the performance of DA, but also has the potential to improve the performance EnKFs-based IR DA under more realistic situations. That potential deserves future investigation.

&nbsp; &nbsp;


## Repository Description
This repository contains the following materials for Chan's original research article manuscript ("Improving Infrared Radiance Ensemble Data Assimilation though Mitigating Deleterious Non-Gaussian Artifacts"). 
1) `DART_modified`: A modified version of the NCAR Data Assimilation Research Testbed (DART)
2) `WRF_v4.5.1`: Weather Research and Forecast model version 4.5.1
3) `Expt_Performance_Diagnostics`: Files containing the performance diagnostics used to generate the figures in the manuscript
4) `Plot_Diagnostics`: Some of the Python scripts used to generate the figures in the manuscript.

Due to the sheer volume of raw simulation data produced in this study (~11 terabytes), raw simulation data is not archived on this repository.

If you have questions, please email Man-Yau (Joseph) Chan at the Ohio State University's Department of Geography.