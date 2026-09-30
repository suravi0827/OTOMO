clear;
close all;
clc;

csvFile = fullfile(pwd,'OTOMO_Gamma_Sensitivity','OTOMO_Gamma_Sensitivity_Summary.csv');
figureFolder = fullfile(pwd,'OTOMO_Gamma_Sensitivity','Sensitivity_Figures');
datasetOrder = ["NASA","RELINK","SOFTLAB","AEEEM"];

plot_sensitivity_figures(csvFile,figureFolder,datasetOrder);
