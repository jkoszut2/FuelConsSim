function [outputTime, outputFuel, ...
        straightsStateData, accStateData] = fuelSim(Log, ...
                                        testEntireRace,analyzeOrigError, ...
                                        modeCoupled, modifierArray)
r = 9.875*0.0254;
gears = 2.073*[31/12 32/16 30/18 26/18 27/21 23/20]; % crank:sprocket
if (strcmp(Log,'FSAEM_Endurance_20190511-1260803.csv'))
    fdr = 36/11;
    % 2.6.1.10e
    % All FEPW values multipled by 1.015 to account for average air temp 
    % compensation during FSAEM 2019
    % All values also multipled by 1.014 to account for average lambda 
    % closed loop control effort during FSAEM 2019
    FEPW = [0 2.4192,2.5564,3.10415,3.7058,3.7513,3.9263,3.9697,4.0796, ...
        4.0922,4.0957,4.01625,4.08975,4.2679,4.4856,4.7614,5.319475, ...
        5.327816,5.19843,5.07126,4.98575,4.956,5.01025,4.74985,4.45573, ...
        4.421025]*1.015*1.014; % ms
    % 0.3 lbft loss assumed due to changes from 2.6.1.8e to 2.6.1.10e
    Torque_lbft = [0.3 15.17 15.93 22.67 26.13 27.57 29.07 28.83 30.07 31.17 ...
    30.97 31.6 31.53 31.83 33 35.33 38.1 37.67 36.67 35 33.67 31.83 31 ...
    29.13 27.66 25.5]-0.3; % lb-ft
elseif (strcmp(Log,'FSAEL_Endurance_20190622-1260800_MATLABfix.csv'))
    fdr = 34/11;
    % 2.6.1.15e
    % All FEPW values multipled by 0.992 to account for average air temp 
    % compensation during FSAEL 2019
    % All FEPW values also multipled by 0.984 to account for average lambda 
    % closed loop control effort during FSAEL 2019
    FEPW = [0 2.4192 2.5564 3.10415 3.6288 3.7513 3.84335 3.9305 4.0229 ...
        4.0194 4.0145 3.92875 4.02325 4.2308 4.4548 4.7124 5.2346 5.215 ...
        5.13846 4.93686 4.81425 4.83 4.95308 4.7047 4.393783 4.34665]*0.992*0.984; % ms
    % 0.2 lbft loss assumed due to changes from 2.6.1.10e to 2.6.1.15e
    % 0.3 lbft has to be accounted for 2.6.1.8e --> 2.6.1.10e, so really
    % have to assume loss of 0.5 lbft from 2.6.1.08e --> 2.6.1.15e
    Torque_lbft = [0.5 15.17 15.93 22.67 26.13 27.57 29.07 28.83 30.07 31.17 ...
    30.97 31.6 31.53 31.83 33 35.33 38.1 37.67 36.67 35 33.67 31.83 31 ...
    29.13 27.66 25.5]-0.5; % lb-ft
end

% Other fuel maps
% 2.6.1.8e
% FEPW = [0 2.578 2.535 3.163 3.661 3.811 4.085 3.981 4.145 4.176 4.244 ...
%     4.208 4.31 4.45	4.765 5.032	5.58 5.651 5.665 5.365 5.47	5.354 5.376 ...
%     5.144 4.915]; % ms
% redline = 13500; % rpm

% 2.6.1.8p
% FEPW = [0 2.742 2.607 3.201 3.797 3.953 4.245 4.13 4.309 4.363 4.391 ...
%     4.174 4.269	4.406 4.707	4.881 5.524	5.629 5.614	5.322 5.401	5.271 ...
%     5.202 5.053	4.811]; % ms
% Torque_lbft = [0 16.26666667 16.93333333 22.66666667 26.43333333 ...
%     27.63333333	29.96666667	28.63333333	30.13333333	31.76666667	 ...
%     30.93333333	31.96666667	31.8	32.23333333	33.36666667	35.66666667	...
%     38.8 37.96666667 36.9 35.23333333 34 31.93333333 30.66666667 ...
%     28.96666667	27.06666667]; % lb-ft

fprintf('=======================================\n')
fprintf('------- Start of new simulation -------\n')
fprintf('=======================================\n')

t0 = tic();
fprintf('Running fuelSim\n')
% If desire to calculate fuel consumption over an entire log file, turn on
% the 'testEntireRace' boolean in the fuelSim function call.
% If testEntireRace is on, analyzeOrigError is not used in the function, so
% it can be set to 0 or any other value as well
if ~exist('testEntireRace','var') || isempty(testEntireRace)
    testEntireRace = 0;
end
% Use analyzeOrigError to compare the error between the predicted results
% (calculated using accelSim) and the original results (calculated 
% using euler's approximation with the raw logged data)
if ~exist('analyzeOrigError','var') || isempty(analyzeOrigError)
    analyzeOrigError = 0;
end

if ~exist('modifierArray','var') || isempty(modifierArray)
    modifierArray = [1 1];
end

% modeCoupled determines if the front and rear wheels are coupled in the
% acceleration sim. Decoupling the front and rear wheels improves the
% fuel consumption estimation by matching RPM to the RPM in the logged
% data. However, this should be turned off when running sweeps for points
% in order to use a more theoretically correct assumption which is no wheel
% slip for straight line acceleration at speeds > ~30 mph.
if ~exist('modeCoupled','var') || isempty(modeCoupled)
    modeCoupled = 1;
end

% Set value of fuelConsIndex to the index that is used for storing
% cumulative fuel consumption in accelSim
fuelConsIndex = 10;

% Process csv
t1 = tic();
fprintf('Running processLog\n')
[~,LoggedData2,~,index_RPM,index_TP, ...
    index_FEPW,index_WSSFL,index_WSSRL,fs] = processLog(Log);
dt1 = toc(t1);
fprintf('Time to run processLog: %0.2f sec\n', dt1)

% Process throttle position data
% Set constraints for minimum throttle position and duration at that
% position in order to be considered a straight
tpThresh = 90; % percent
if testEntireRace
    tpThresh = 0; % percent
end
timeThresh = 1; % seconds
sampleThresh = timeThresh*fs; % samples

% Find straights
t2 = tic();
fprintf('Running findStraights\n')
straights = findStraights(Log, LoggedData2, index_TP, index_WSSFL, ...
                        index_WSSRL, index_RPM, tpThresh, sampleThresh, ...
                        testEntireRace);
dt2 = toc(t2);
fprintf('Time to run findStraights: %0.2f sec\n', dt2)

% Plot all TP data and overlay TP during straights
% Uncomment plot commands in for loop below to get straights overlay
%{
time = LoggedData2(:,index_Time);
plot(time,tpData) % All TP data
hold on
for i=1:length(straights(:,1))
    plot(time(straights(i,1):straights(i,2)),tpData(straights(i,1):straights(i,2)),'r')
    hold on
end
%}

% Calculate rear wheel speed based off rpm
% NOTE: This is currently not used for anything. It was used for analyzing 
% the error between the logged rear wheel speed reading and the theoretical
% rear wheel speed based off RPM and a gear estimate. This is not perfect 
% since drivers do not always perform optimal shifts and thus there are
% occasional errors in predicted gear versus actual gear. Gear position was
% not logged and thus exact gears are not known.
WSSRL_RPM = zeros(length(LoggedData2(:,1)),1);
for o = 1:length(LoggedData2)
    gear = findGear3(LoggedData2(o,index_RPM),LoggedData2(o,index_WSSFL),13500);
    WSSRL_RPM(o) = LoggedData2(o,index_RPM)*60*2*pi*r/fdr/gears(gear)/1609.4;
end

% Calculate distance traveled and fuel consumed during each straight from
% logged data
accelTimes = zeros(length(straights(:,1)),1); % miles
distances = zeros(length(straights(:,1)),1); % miles
fuelCons = zeros(length(straights(:,1)),1); % cc
% Injector flow rate below from Ford injector characterization data
% Ford injector data provided in mass flow rate, so needs to be converted
% to volumetric flow rate. Also needs to be adjusted based off Bernoulli's
% principle as Ford testing was performed at 270 kPa
% Average of gauge fuel pressure from FSAEM and FSAEL is 44.03psi
% 44.03psi --> 303.57616 kPa
injFlowRate = 2.58*60/0.738*sqrt(303.57616/270); % cc/min
for p=1:length(straights(:,1))
    sampleStart = straights(p,1);
    sampleEnd = straights(p,2);
    currDistTraveled = 0;
    currFuelCons = 0;
    cylCount = 4;
    for sampIn = sampleStart:(sampleEnd-1)
        % Calculate Distance Using Midpoint Sum
        WSSFL = LoggedData2(sampIn,index_WSSFL); % vehicle velocity in mph
        WSSRL = LoggedData2(sampIn,index_WSSRL); % vehicle velocity in mph
        WSSRL2 = WSSRL_RPM(sampIn); % vehicle velocity in mph
        instDistTraveled = 0.5*(LoggedData2(sampIn,index_WSSFL)+LoggedData2(sampIn+1,index_WSSFL))/3600*(1/fs); % instantaneous distance in miles
        currDistTraveled = currDistTraveled + instDistTraveled; % distance traveled in miles
        % Calculate Fuel Consumption Using Midpoint Sum
        currFPW = LoggedData2(sampIn,index_FEPW)/1000; % injector pulse width in seconds
        currFPWmax = 1/(LoggedData2(sampIn,index_RPM)/60/2); % pulse width correspond to 100% duty cycle
        currFDC = currFPW/currFPWmax*100; % injector duty cycle
        nextFPW = LoggedData2(sampIn+1,index_FEPW)/1000; % injector pulse width in seconds
        nextFPWmax = 1/(LoggedData2(sampIn+1,index_RPM)/60/2); % pulse width correspond to 100% duty cycle
        nextFDC = nextFPW/nextFPWmax*100; % injector duty cycle
        instFuelCons = 0.5*(currFDC+nextFDC)/100*injFlowRate/60*(1/fs)*cylCount; % instantaneous fuel consumed in cc
        % Check to make sure over-run fuel cut (ORFC) isn't active
        % NOTE: These should NOT be getting used by the script when
        % testEntireRace is off
        if (strcmp(Log,'FSAEM_Endurance_20190511-1260803.csv'))
            ORFC_RPM = 5500; % RPM above which ORFC is active
        elseif (strcmp(Log,'FSAEL_Endurance_20190622-1260800_MATLABfix.csv'))
            ORFC_RPM = 6000; % RPM above which ORFC is active
        end
        ORFC_TP = 20; % Throttle position below which ORFC is active
        if ((LoggedData2(sampIn,index_TP) < ORFC_TP) && ...
                (LoggedData2(sampIn,index_RPM)>ORFC_RPM))
            % Over-run fuel cut is active
            instFuelCons = 0;
        end
        currFuelCons = currFuelCons+instFuelCons; % fuel consumed in cc
    end
    accelTimes(p) = (sampleEnd-sampleStart)/fs;
    distances(p) = currDistTraveled;
    fuelCons(p) = currFuelCons;
    % Uncomment below to get straights overlay for throttle position plot
    % plot(time(straights(i,1):straights(i,2)),tpData(straights(i,1):straights(i,2)),'r')
    % hold on
end

% Overall distance traveled and fuel consumed
% Not required, only used as a method of double checking calculations above
dt2 = 0; % distance in miles
fc2 = 0; % fuel consumption in cc
for q = 1:length(LoggedData2)
    WSSFL = LoggedData2(q,index_WSSFL); % miles per hour
    dt2 = dt2+WSSFL/3600*(1/fs);
    FPW = LoggedData2(q,index_FEPW)/1000; % injector pulse width in seconds
    FPWmax = 1/(LoggedData2(q,index_RPM)/60/2); % pulse width correspond to 100% duty cycle
    FDC = FPW/FPWmax*100; % injector duty cycle
    fc2 = fc2+FDC/100*injFlowRate/60*(1/fs)*cylCount; % fuel consumed in cc
end

% Find initial conditions for each straight
initRPM = zeros(length(straights(:,1)),1);
initWSSFL = zeros(length(straights(:,1)),1);
initWSSRL = zeros(length(straights(:,1)),1);
for v = 1:length(straights(:,1))
    initRPM(v) = LoggedData2(straights(v,1),index_RPM);
    initWSSFL(v) = LoggedData2(straights(v,1),index_WSSFL);
    initWSSRL(v) = LoggedData2(straights(v,1),index_WSSRL);
    % Correct for wrong rear wheel speed calibration
    initWSSRL(v) = initWSSRL(v);
end

% Create a table with the characteristics of each straight
% Purely for convenience of viewing details for each straight
straightsStateData = table(straights,accelTimes,distances, ...
    fuelCons,initRPM,initWSSFL,initWSSRL);

% Create cell to store states of each acceleration simulation
accStateData = cell(length(straights(:,1)),1);
% NOTE: For header below, copy from accelSim
header = {'time' 'pos' 'vel' 'acc' 'gear' 'rpm' 'wssfl' 'wssrl', ...
 'fuelcons', 'cumfuelcons' 'Torque_Nm' 'wheel_Torque_Nm' 'slip' , ...
'Weight_Transfer_Acc' 'Fx' 'MaxFx' 'Drag' 'Downforce' 'FEPW' 'IDC'};

if ~testEntireRace
    % Run accel sim over each straight to get predicted state behavior
    % Accel sim takes initial velocity in meters/sec, so input to function must
    % be converted from mph
    predAccelTimes = zeros(length(straights(:,1)),1); % msec
    predFuelCons = zeros(length(straights(:,1)),1); % cc
    t3 = tic();
    fprintf('Running accelSim for each straight\n')
    imax = length(straights);
    prog = 0;
    fprintf(1,'Computation Progress:  %3d%%\n',prog);
    for i = 1:length(straights(:,1))
        prog = ( 100*(i/imax) );
        fprintf(1,'\b\b\b\b\b\b%3.0f %%\n%',prog);
        pause(0.1); % Deleting 4 characters (The three digits and the % symbol)
        makePlots=0;
        % Uncomment code below for debugging help
        %{
        if i == 22
            makePlots = 1;
            fprintf('Max RPM: %0.1f\n', max(LoggedData2(straightsStateData{i,1}(1):straightsStateData{i,1}(2),index_RPM)))
        end
        %}
        maxRPM = max(LoggedData2(straightsStateData{i,1}(1):straightsStateData{i,1}(2),index_RPM));
        if modeCoupled
            [accelTime, stateData] = accelSim(distances(i)*1609.4, ...
                initWSSFL(i)*1609.4/3600, initWSSFL(i)*1609.4/3600, ...
                findGear3(initRPM(i), initWSSRL(i), maxRPM), ...
                13500, ...
                fdr,FEPW,Torque_lbft,modifierArray,makePlots);
        elseif ~modeCoupled
            [accelTime, stateData] = accelSim(distances(i)*1609.4, ...
                initWSSFL(i)*1609.4/3600, findVehV(initRPM(i), findGear3(initRPM(i),initWSSRL(i),maxRPM))*1609.4/3600, ...
                findGear3(initRPM(i), initWSSRL(i), maxRPM), ...
                maxRPM, ...
                fdr,FEPW,Torque_lbft,modifierArray,makePlots);
        end
        
%         [accelTime, stateData] = accelSim(distances(i)*1609.4, ...
%             initWSSFL(i)*1609.4/3600,findGear3(initRPM(i), initWSSFL(i)), ...
%             modifierArray);
        %{
        gear = findGear3(initRPM(i), initWSSRL(i));
        [accelTime, stateData] = accelSim(distances(i)*1609.4, ...
            findVehV(initRPM(i),gear)*1609.4/3600,findGear3(initRPM(i), ...
            initWSSFL(i)), modifierArray);
         %}
%         [accelTime, stateData] = accelSim(distances(i)*1609.4,initWSSFL(i)*1609.4/3600,findGear2(initWSSFL(i)));
%         [accelTime, stateData] = accelSim(distances(i)*1609.4,findVehV(initRPM(i),1)*1609.4/3600,1);  
%         [accelTime, stateData] = accelSim(distances(i)*1609.4,initWSSFL(i)*1609.4/3600,2);
        predAccelTimes(i) = accelTime;
        predFuelCons(i) = max(stateData{:,fuelConsIndex});
        accStateData{i} = cell2table(num2cell(table2array(stateData)),'VariableNames', header);
    end
    dt3 = toc(t3);
    fprintf('Time to run accelSim: %0.2f sec\n', dt3)
    
    % Originally, results were compared against actual acceleration
    % times and fuel consumption from the csv log. After verifying that the
    % acceleration model is relatively accurate, but NOT perfect, all 
    % results are now compared to the predicted acceleration times and
    % fuel consumption from the original log file.
    % Print statistics for comparison between actual and predicted results
    fprintf('========================================\n')
    fprintf('---------------- RESULTS ---------------\n')
    fprintf('========================================\n')
    % If analyzeOrigError on, then compare with actual results from CSV
    if analyzeOrigError
        fprintf('Total actual accel time: %0.2f seconds\n', sum(accelTimes));
        fprintf('Total actual fuel cons: %0.2f cc\n', sum(fuelCons));
    end
    if (strcmp(Log,'FSAEM_Endurance_20190511-1260803.csv'))
        % NOTE: Current time and fuel consumption values are obtained from
        % coupled simulation. Refer to PointsCalcFSAEM for original and
        % decoupled values.
        fprintf('Total predicted original accel time: %0.2f seconds\n', 81.94);
        fprintf('Total predicted original fuel cons: %0.2f cc\n', 577.904557504403);
        % Below if statement allows for analyzing sample-wise error of
        % actual raw accel times and fuel consumption from log file
        % compared to predicted results
        if ~analyzeOrigError
            % Here, set the baselines for the predictions of the original
            % results
            accelTimes = 81.94; %
            fuelCons = 577.904557504403; % 
        end
    elseif (strcmp(Log,'FSAEL_Endurance_20190622-1260800_MATLABfix.csv'))
        % NOTE: Current time and fuel consumption values are obtained from
        % coupled simulation. Refer to PointsCalcFSAEL for original and
        % decoupled values.
        fprintf('Total predicted accel time: %0.2f seconds\n',  123.48);
        fprintf('Total predicted fuel cons: %0.2f cc\n', 777.401686972868);
        % Below if statement allows for analyzing sample-wise error of
        % actual raw accel times and fuel consumption from log file
        % compared to predicted results
        if ~analyzeOrigError
            %       Accel Time [s]   Fuel Cons [cc]
            % FSAEL:    132.80           830.76             
            %       Accel Time [s]   Fuel Cons [cc]   Script Runtime [s]
            %   0.1:    134.5            863.67             22.37
            %  0.01:    129.39           826.25             28.34
            % 0.001:    128.89           822.55             81.75
            % Here, set the baseline for the predictions of the original
            % results
            accelTimes =  123.48; % Compare results against this value
            fuelCons = 777.401686972868; % Compare results against this value
        end
    end
    fprintf('Total predicted new accel time: %0.2f seconds\n', sum(predAccelTimes));
    fprintf('Total predicted new fuel cons: %0.2f cc\n', sum(predFuelCons));
    % Error
    fprintf('----------------------------------------\n')
    fprintf('Difference in accel time: %0.2f seconds\n', ...
        sum(sum(predAccelTimes) - sum(accelTimes)));
    fprintf('Difference in fuel cons: %0.2f cc\n', ...
        sum(sum(predFuelCons) - sum(fuelCons)));
    % Print and plot error mean and standard deviation
    
    if analyzeOrigError
        fprintf('----------------------------------------\n')
        fprintf('Accel time error mean: %0.3f seconds\n', ...
            mean( -((accelTimes)-(predAccelTimes)) ));
        fprintf('Accel time error stdev: %0.3f seconds\n', ...
            std( -((accelTimes)-(predAccelTimes)) ));
        fprintf('Accel time error mean: %0.3f %%\n', ...
            mean( -((accelTimes)-(predAccelTimes)) ./accelTimes*100));
        fprintf('Accel time error stdev: %0.3f %%\n', ...
            std( -((accelTimes)-(predAccelTimes)) ./accelTimes*100 ));
         fprintf('Fuel cons error mean: %0.3f cc\n', ...
            mean( -((fuelCons)-(predFuelCons)) ));
        fprintf('Fuel cons error stdev: %0.3f cc\n', ...
            std( -((fuelCons)-(predFuelCons)) ));
        fprintf('Fuel cons error mean: %0.3f %%\n', ...
            mean( -((fuelCons)-(predFuelCons)) ./fuelCons*100 ));
        fprintf('Fuel cons error stdev: %0.3f %%\n', ...
            std( -((fuelCons)-(predFuelCons)) ./fuelCons*100 ));
        % Plot actual versus predicted results along with respective error
        % Acceleration Times
        figure
        subplot(2,1,1)
        plot(0:1:length(straights(:,1))-1,accelTimes,0:1:length(straights(:,1))-1,predAccelTimes)
        predAccErr = zeros(length(straights(:,1)),1); % percent error
        title('Actual vs Predicted Accel Times')
        xlabel('Sample Number')
        ylabel('Time [sec]')
        legend('Actual', 'Predicted', 'location', 'best')
        grid on
        subplot(2,1,2)
        for i = 1:length(straights(:,1))
            predAccErr(i) = -(accelTimes(i)-predAccelTimes(i))/accelTimes(i)*100;
        end
        plot(0:1:length(straights(:,1))-1,predAccErr)
        hold on
        yline(0, 'k--');
        ylim([-10 10])
        title('Acceleration Time Error (Predicted-Actual)/Actual')
        xlabel('Sample Number')
        ylabel('Percent Error [%]')
        legend('Percent Error', 'Target', 'location', 'best')
        grid on
        if (strcmp(Log,'FSAEM_Endurance_20190511-1260803.csv'))
            sgt = sgtitle('FSAEM 2019');
            sgt.FontSize = 14;
        elseif (strcmp(Log,'FSAEL_Endurance_20190622-1260800_MATLABfix.csv'))
            sgt = sgtitle('FSAEL 2019');
            sgt.FontSize = 14;
        end
        % Fuel Consumption
        figure
        subplot(2,1,1)
        plot(0:1:length(straights(:,1))-1,fuelCons,0:1:length(straights(:,1))-1,predFuelCons)
        hold on
        predFuelErr = zeros(length(straights(:,1)),1); % percent error
        title('Actual vs Predicted Fuel Consumption')
        xlabel('Sample Number')
        ylabel('Time [sec]')
        legend('Actual', 'Predicted', 'location', 'best')
        grid on
        subplot(2,1,2)
        for i = 1:length(straights(:,1))
            predFuelErr(i) = -(fuelCons(i)-predFuelCons(i))/fuelCons(i)*100;
        end
        plot(0:1:length(straights(:,1))-1,predFuelErr)
        hold on
        yline(0, 'k--');
        ylim([-10 10])
        title('Fuel Consumption Error (Predicted-Actual)/Actual')
        xlabel('Sample Number')
        ylabel('Percent Error [%]')
        legend('Percent Error', 'Target', 'location', 'best')
        grid on
        if (string(Log) == 'FSAEM_Endurance_20190511-1260803.csv')
            sgt = sgtitle('FSAEM 2019');
            sgt.FontSize = 14;
        elseif (string(Log) == 'FSAEL_Endurance_20190622-1260800_MATLABfix.csv')
            sgt = sgtitle('FSAEL 2019');
            sgt.FontSize = 14;
        end
    end
    
    outputTime = sum(predAccelTimes);
    outputFuel = sum(predFuelCons);
else
    fprintf('Accel times for entire race:\n')
    disp(accelTimes)
    fprintf('Total: %0.1f\n', sum(accelTimes))
    fprintf('Fuel consumption for entire race:\n')
    disp(fuelCons)
    fprintf('Total: %0.1f\n', sum(fuelCons))
    outputTime = sum(accelTimes);
    outputFuel = sum(fuelCons);
end
dt0 = toc(t0);
fprintf('Time to run fuelSim: %0.2f sec\n', dt0)

% Functions
function gear = findGear3(currRPM, currWSSrear, maxRPM)
% Same as findGear2 but two major changes:
% 1. Prioritizes RPM over vehicle velocity to find gear
% 2. No optimal gear selection
% Use this function to initialize the accel sim
% currRPM is current RPM
% currWSSrear is current rear wheel speed reading in miles per hour
n = 6;
err = 13500*ones(6,1); % Initialize an error array
% Loop over each possible gear and determine the error between actual
% RPM and theoretical RPM for that gear. If a gear will cause an RPM
% greater than redline, do not calculate the error, simply leave the error
% to be as redline to ensure it does not become a minimum
redlineLoc = maxRPM;
while ((n>=1) && (currWSSrear/60*1609.4/(2*pi*r)*fdr*gears(n) < redlineLoc))
    predRPM = currWSSrear/60*1609.4/(2*pi*r)*fdr*gears(n);
    err(n) = abs(predRPM-currRPM);
    n = n-1;
end
% Use the index corresponding to the minimum of the error array as the gear
[~,minIndex] = min(err);
gear = minIndex;
end

function VehV = findVehV(currRPM, currGear)
% Finds vehicle velocity in miles per hour
% global r fdr gears
VehV = currRPM*60/1609.4*(2*pi*r)/gears(currGear)/fdr;
end

function straights = findStraights(Log, LoggedData2,index_TP, index_WSSFL, ...
                                index_WSSRL, index_RPM, tpThresh, sampleThresh, ...
                                testEntireRace)
% This function finds all the "straights" in a log file
% The sample points corresponding to the starts and ends of the straights
% are stored in the straights array 

tpData = LoggedData2(:,index_TP);
straights = []; % initialize array 
i = 1;
while ( (1<=i) && (i<=length(tpData)) )
    tstpnt = i;
    % Set a limit on how much higher rear wheel speeds can be than fronts
    % to avoid issues with high wheel slip.
    if ~testEntireRace
        if ( (tpData(tstpnt) >= tpThresh) && ...
                (LoggedData2(tstpnt,index_WSSFL) > 30) && ...
                ((LoggedData2(tstpnt,index_WSSRL)-LoggedData2(tstpnt,index_WSSFL)) < 15))
            start = tstpnt;
            while ( (tstpnt < length(tpData)-1) && (tpData(tstpnt+1) >= tpThresh) )
                tstpnt = tstpnt+1;
            end
            stop = tstpnt;
            % If the straight meets the minimum time requirement, 
            % then append the straight start stop points to the 
            % straights array.
            if ((stop-start >= sampleThresh) && ...
                    mean(LoggedData2(start:stop,index_WSSRL)-LoggedData2(start:stop,index_WSSFL)) < 10 )
                    straights = [straights; [start stop]];
            end
        end
    else
        % To calculate fuel over an entire endurance run, use all samples
        % occuring after the predefined sample point
        if (strcmp(Log,'FSAEM_Endurance_20190511-1260803.csv'))
            strt = 2350;
        elseif (strcmp(Log,'FSAEL_Endurance_20190622-1260800_MATLABfix.csv'))
            strt = 3470;
        end
        if ((tstpnt-1) >= strt)
            if tpData(tstpnt) >= tpThresh
                start = tstpnt;
                % RPM constraint is to stop during driver change (engine off)
                % and after race (engine off again). Reason for this is
                % that the ECU continues to log for a short bit after the
                % engine turns off but the battery supply to the ECU is
                % still active (since the car is still on).
                while ((tstpnt < length(tpData)-1) && (tpData(tstpnt+1) >= tpThresh) ...
                        && (LoggedData2(tstpnt,index_RPM) > 10) )
                    tstpnt = tstpnt+1;
                end
                stop = tstpnt;
                % If the straight meets the minimum time requirement, 
                % then append the straight start stop points to the 
                % straights array.
                if ( (stop-start >= sampleThresh) )
                        straights = [straights; [start stop]];
                end
            end
        end
    end
    i = tstpnt+1;
end
end

end
