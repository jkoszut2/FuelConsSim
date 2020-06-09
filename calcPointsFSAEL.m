function [deltaFSAELend,deltaFSAELeff] = PointsCalcFSAEL(newFSAELtime, ...
                                            newFSAELfuel)
% Endurance + Efficiency Points Calculator
% Returns points gain for endurance event and efficiency score
% UNFILTERED DECOUPLED Front/Rear Wheel Speed Original Times and Fuel Consumption
% straightsFSAELacc = 163.33; % Original predicted accel time [s] for straights
% straightsFSAELfuel = 1022.00; % Original predicted fuel consumption [cc] for straights
% DECOUPLED Front/Rear Wheel Speed Original Times and Fuel Consumption
% straightsFSAELacc = 127.13; % Original predicted accel time [s] for straights
% straightsFSAELfuel = 820.812189911384; % Original predicted fuel consumption [cc] for straights
% COUPLED Front/Rear Wheel Speed Original Times and Fuel Consumption
straightsFSAELacc = 123.48; % Original predicted accel time [s] for straights
straightsFSAELfuel = 777.401686972868; % Original predicted fuel consumption [cc] for straights
fprintf('=======================================\n')
fprintf('---------------- FSAEL ----------------\n')
fprintf('=======================================\n')
% FSAEL - Endurance
numlaps = 15;
penalties = 30;
Tyour = (1294.611+penalties)/(numlaps-2);
Tmin = 1364.599/numlaps; 
Tmax = 1978.669/numlaps;
ptsFSAEL = 250*((Tmax/Tyour)-1)/((Tmax/Tmin)-1)+25;
fprintf('Original Endurance Points: %0.2f pts\n', ptsFSAEL)

% FSAEL - Predicted Endurance
deltaFSAELend = newFSAELtime - straightsFSAELacc;
numlaps = 15;
penalties = 30;
Tyour = (1294.611+penalties+deltaFSAELend)/(numlaps-2);
Tmin = 1364.599/numlaps; 
Tmax = 1978.669/numlaps;
ptsFSAELnew = 250*((Tmax/Tyour)-1)/((Tmax/Tmin)-1)+25;
fprintf('New Endurance Points: %0.2f pts\n', ptsFSAELnew)
deltaFSAELend = ptsFSAELnew-ptsFSAEL;
fprintf('Endurance Points Gain: %0.2f pts\n', deltaFSAELend)

% FSAEL - Efficiency
numlaps = 15;
origFSAELfuel = 3.465; % liters
factorCO2 = 2.31; % kg of CO2 per liter of gasoline
adjCO2 = origFSAELfuel*factorCO2/(numlaps-2); % your adjusted CO2 per lap
minCO2 = 2.018*2.31/15; % best adjusted CO2 per lap out of eligible teams
Effmin = 0.243;
Effmax = 0.904;
Effyour = Tmin/Tyour*minCO2/adjCO2;
ptsFSAELeff = 100*((Effmin/Effyour)-1)/((Effmin/Effmax)-1);
% Result is actually too high by 0.1 compared to PDF due to rounding errors
fprintf('Original Efficiency Points: %0.2f pts\n', ptsFSAELeff)

% FSAEL - Predicted Efficiency
numlaps = 15;
deltaFSAELfuel = (newFSAELfuel - straightsFSAELfuel)/1000; % cc to liters
adjCO2 = (origFSAELfuel+deltaFSAELfuel)*factorCO2/(numlaps-2); % your adjusted CO2 per lap
Effyour = Tmin/Tyour*minCO2/adjCO2;
ptsFSAELeffnew = 100*((Effmin/Effyour)-1)/((Effmin/Effmax)-1);
fprintf('New Efficiency Points: %0.2f pts\n', ptsFSAELeffnew)
deltaFSAELeff = ptsFSAELeffnew-ptsFSAELeff;
fprintf('Efficiency Points Gain: %0.2f pts\n', deltaFSAELeff)


fprintf('=======================================\n')
fprintf('--- FSAEL Net Points Gain: %0.2f pts ---\n', deltaFSAELend + deltaFSAELeff)
fprintf('=======================================\n')

