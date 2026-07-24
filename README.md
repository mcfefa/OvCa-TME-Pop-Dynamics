# OvCa-TME-Pop-Dynamics
Population dynamics modeling of sensitive and resistant ovarian cancer populations along with fibroblast and macrophages

This repository contains the files and data needed to parameterize existing parameters and new parameters explain how macrophages and fibroblasts affect sensitive and resistant population growth dynamics for two ovarian cancer cell lines, A2780 and Tyk-nu. The repository is divided into two folders one for M2TAM (macrophages) and one for Fibroblasts. The analysis for the two cell lines is similar, even though results are different. 

## Macrophage effects 
To replicate the code, download the data folder and each of the cell lines folders (A2780withM2TAM and TykwithM2TAM) inside the M2TAM folder. 
To estimate the growth rates for sensitive and resistant populations under M2TAM conditioned media, run A2780withM2TAM_in_co_culture.m and TykwithM2TAM_in_co_culture.m MATLAB scripts for A2780 and Tyk-nu cell populations, respectively. 

## Fibroblast effects
To replicate the code, download the data, ACAFandTCAF_parameterization, and each of the cell lines folders (A2780withACAF and TykwithTCAF) inside the fibroblast folder. To add the fibroblast population to the sensitive-resistant growth dynamics run the scripts in the following order:
- ACAFandTCAF_parameterEstimation.m MATLAB script to estimate for the logistic growth parameters for the ACAF and TCAF populations
- ACAF_inco_culture.m and TCAF_inco_culture.m MATLAB sript to estimate parameters and BIC for models A and B for ACAFs and TCAFs effects on A2780 and Tyk-nu co-cultures

