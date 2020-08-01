% June 7, 2020
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

t5 = tic();
fprintf('Running parfor loop\n')

% torque = MBT*(-2.65*lambdaDev.^2 - 0.0437*lambdaDev + 1);
factor1 = 0:0.1:0.4; % lambdaDev
FSAEMpts = zeros(length(factor1),2);
FSAELpts = zeros(length(factor1),2);
FSAEMsimData = cell(length(factor1),1);
FSAELsimData = cell(length(factor1),1);
parfor k = 1:length(factor1)
    fprintf('=======================================\n')
    fprintf('------------- %0.0f out of %0.0f -------------\n', k, length(factor1))
    % FSAEM
    [newTimeFSAEM, newFuelFSAEM, straightsFSAEM, accelSimDataFSAEM] = fuelSim(LogFSAEM,0,0,1,[factor1(k)]);
    [FSAEMendpts,FSAEMeffpts] = calcPointsFSAEM(newTimeFSAEM,newFuelFSAEM);
    FSAEMpts(k,:) = [FSAEMendpts FSAEMeffpts];
    FSAEMsimData{k} = accelSimDataFSAEM;
    % FSAEL
    [newTimeFSAEL, newFuelFSAEL, straightsFSAEL, accelSimDataFSAEL] = fuelSim(LogFSAEL,0,0,1,[factor1(k)]);
    [FSAELendpts,FSAELeffpts] = calcPointsFSAEL(newTimeFSAEL,newFuelFSAEL);
    FSAELpts(k,:) = [FSAELendpts FSAELeffpts];
    FSAELsimData{k} = accelSimDataFSAEL;
end

dt5 = toc(t5);
fprintf('Time to run parfor loop: %0.2f sec\n', dt5)
