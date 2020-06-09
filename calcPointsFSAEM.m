function [deltaFSAEMend,deltaFSAEMeff] = PointsCalcFSAEM(newFSAEMtime, ...
                                            newFSAEMfuel)
% Endurance + Efficiency Points Calculator
% Returns points gain for endurance event and efficiency score
% UNFILTERED DECOUPLED Front/Rear Wheel Speed Original Times and Fuel Consumption
% straightsFSAEMacc = 100.74; % Original predicted accel time [s] for straights
% straightsFSAEMfuel = 704.11; % Original predicted fuel consumption [cc] for straights
% DECOUPLED Front/Rear Wheel Speed Original Times and Fuel Consumption
% straightsFSAEMacc = 84.38; % Original predicted accel time [s] for straights
% straightsFSAEMfuel = 606.126157040443; % Original predicted fuel consumption [cc] for straights
% COUPLED Front/Rear Wheel Speed Original Times and Fuel Consumption
straightsFSAEMacc = 81.94; % Original predicted accel time [s] for straights
straightsFSAEMfuel = 577.904557504403; % Original predicted fuel consumption [cc] for straights
fprintf('=======================================\n')
fprintf('---------------- FSAEM ----------------\n')
fprintf('=======================================\n')
% FSAEM - Endurance
numlaps = 11;
Tyour = 1533.958/numlaps;
Tmin = 1267.742/numlaps;
Tmax = 1838.226/numlaps; % 145% of Tmin
ptsFSAEM = 250*((Tmax/Tyour)-1)/((Tmax/Tmin)-1)+25;
fprintf('Original Endurance Points: %0.2f pts\n', ptsFSAEM)

% FSAEM - Predicted Endurance
deltaFSAEMend = newFSAEMtime - straightsFSAEMacc;
numlaps = 11;
Tyour = (1533.958+deltaFSAEMend)/numlaps;
Tmin = 1267.742/numlaps;
Tmax = 1838.226/numlaps; % 145% of Tmin
ptsFSAEMnew = 250*((Tmax/Tyour)-1)/((Tmax/Tmin)-1)+25;
fprintf('New Endurance Points: %0.2f pts\n', ptsFSAEMnew)
deltaFSAEMend = ptsFSAEMnew-ptsFSAEM;
fprintf('Endurance Points Gain: %0.2f pts\n', deltaFSAEMend)

% FSAEM - Efficiency
numlaps = 11;
origFSAEMfuel = 4.112; % liters
factorCO2 = 2.31; % kg of CO2 per liter of gasoline
adjCO2 = origFSAEMfuel*factorCO2/numlaps; % your adjusted CO2 per lap
minCO2 = 2.650*2.31/11; % best adjusted CO2 per lap out of eligible teams
Effmin = 0.319;
Effmax = 0.989;
Effyour = Tmin/Tyour*minCO2/adjCO2;
ptsFSAEMeff = 100*((Effmin/Effyour)-1)/((Effmin/Effmax)-1);
% Result is actually too high by 0.1 compared to PDF due to rounding errors
fprintf('Original Efficiency Points: %0.2f pts\n', ptsFSAEMeff)

% FSAEM - Predicted Efficiency
numlaps = 11;
deltaFSAEMfuel = (newFSAEMfuel - straightsFSAEMfuel)/1000; % cc to liters
adjCO2 = (origFSAEMfuel+deltaFSAEMfuel)*factorCO2/numlaps; % your adjusted CO2 per lap
Effyour = Tmin/Tyour*minCO2/adjCO2;
ptsFSAEMeffnew = 100*((Effmin/Effyour)-1)/((Effmin/Effmax)-1);
fprintf('New Efficiency Points: %0.2f pts\n', ptsFSAEMeffnew)
deltaFSAEMeff = ptsFSAEMeffnew-ptsFSAEMeff;
fprintf('Efficiency Points Gain: %0.2f pts\n', deltaFSAEMeff)


fprintf('=======================================\n')
fprintf('--- FSAEM Net Points Gain: %0.2f pts ---\n', deltaFSAEMend + deltaFSAEMeff)
fprintf('=======================================\n')

