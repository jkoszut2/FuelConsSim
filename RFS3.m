% June 7, 2020
% Joe Koszut
% Script used for sweeping vehicle parameters in fuelSim.

% DIRECTIONS:
% 1. Choose one of the log files below
% 2. Define a sweep
% 3. Go to accelSim and set up the vehice parameters of interest so that 
%    they can be changed during the sweep
% 3. Set up data mapping accordingly in the for loop
% 4. Run script

clc
clear all
close all

format shortG

% CHOOSE LOG FILE
% Can use just one log if only performing sweep for either FSAEM or FSAEL. 
% Can also use both log files, simply set up the sweep accordingly.
% NOTE: First FSAEL log file has a formatting issue and is corrected with
% the file below it
LogFSAEM = 'FSAEM_Endurance_20190511-1260803.csv';
% LogFSAEL = 'FSAEL_Endurance_20190622-1260800.csv';
LogFSAEL = 'FSAEL_Endurance_20190622-1260800_MATLABfix.csv';

% SWEEP DEFINITION
% Use factor1 to pass in a modifier to accelSim and change a vehicle
% parameter
% Example: Pass in how much lambda should change and set up an
% equation in accelSim to recalculate torque. Equation example:
% torque = MBT*(-2.65*lambdaDev.^2 - 0.0437*lambdaDev + 1);
factor1 = 0:0.01:0.4; % lambdaDev

% INITIALIZE ARRAYS TO STORE SWEEP RESULTS
% Arrays below for storing total endurace lap time and fuel consumed
FSAEMresults = zeros(length(factor1),2);
FSAELresults = zeros(length(factor1),2);
% Arrays below for storing endurace and effciency points
FSAEMpts = zeros(length(factor1),2);
FSAELpts = zeros(length(factor1),2);
% Cells below for storing parameters of each straight and each accelSim run
FSAEMsimData = cell(length(factor1),1);
FSAELsimData = cell(length(factor1),1);
% Time loop
t5 = tic();
fprintf('Running parfor loop\n')
parfor k = 1:length(factor1)
    fprintf('=======================================\n')
    fprintf('------------- %0.0f out of %0.0f -------------\n', k, length(factor1))
    % FSAEM
    [newTimeFSAEM, newFuelFSAEM, straightsFSAEM, accelSimDataFSAEM] = fuelSim(LogFSAEM,0,0,1,[factor1(k)]);
    [FSAEMendpts,FSAEMeffpts] = calcPointsFSAEM(newTimeFSAEM,newFuelFSAEM);
    FSAEMresults(k,:) = [newTimeFSAEM newFuelFSAEM];
    FSAEMpts(k,:) = [FSAEMendpts FSAEMeffpts];
    FSAEMsimData{k} = accelSimDataFSAEM;
    % FSAEL
    [newTimeFSAEL, newFuelFSAEL, straightsFSAEL, accelSimDataFSAEL] = fuelSim(LogFSAEL,0,0,1,[factor1(k)]);
    [FSAELendpts,FSAELeffpts] = calcPointsFSAEL(newTimeFSAEL,newFuelFSAEL);
    FSAELresults(k,:) = [newTimeFSAEL newFuelFSAEL];
    FSAELpts(k,:) = [FSAELendpts FSAELeffpts];
    FSAELsimData{k} = accelSimDataFSAEL;
end

dt5 = toc(t5);
fprintf('Time to run parfor loop: %0.2f sec\n', dt5)


% Run SweepPlotGeneral to produce a general plot of the sweep.
run SweepPlotGeneral
