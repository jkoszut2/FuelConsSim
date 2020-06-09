% May 28, 2020
% Joe Koszut
% Directions:
% Choose one of the log files below and run script

clc
clear all
close all

% run params
format shortG

% NOTE: First FSAEL log file has formatting issue and is corrected with
% the file below it
LogFSAEM = 'FSAEM_Endurance_20190511-1260803.csv';
% Log = 'FSAEL_Endurance_20190622-1260800.csv';
LogFSAEL = 'FSAEL_Endurance_20190622-1260800_MATLABfix.csv';

% fuelSim(<Log Name>, <Entire Race>, <Analyze Error>, <Coupled Mode>)
[newTimeFSAEM, newFuelFSAEM, straightsFSAEM] = fuelSim(LogFSAEM,0,0,1);
[FSAEMendpts,FSAEMeffpts] = calcPointsFSAEM(newTimeFSAEM,newFuelFSAEM);
[newTimeFSAEL, newFuelFSAEL, straightsFSAEL] = fuelSim(LogFSAEL,0,0,1);
[FSAELendpts,FSAELeffpts] = calcPointsFSAEL(newTimeFSAEL,newFuelFSAEL);

