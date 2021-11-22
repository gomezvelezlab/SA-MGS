# Spectral Analysis of Multiscale Groundwater Systems (SA-MGS)

This repository contains the Matlab codes for the estimation of groundwater fluxes using spectral analysis with dominant frequencies. See Perez et al., (2021) (Under review in WRR) "Identification of Characteristic Spatial Scales to Improve the Performance of Analytical Spectral Solutions to the Groundwater Flow Equation".

The Matlab functions included in this repository are:
- SS_DomFreq.m: Selection of Dominant Frequencies to be used on SS_Top.m
- SS_Top.m: Estimation of the Fourier surface to be used on SS_Vel.m
- SS_Vel.m: Exact Solution to groundwater flow in 3D based on dominant frequencies
- SS_Example.m: Examples to execute SS_DomFreq, SS_Top ans SS_Vel codes 
