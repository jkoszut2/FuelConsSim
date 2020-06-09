function [LoggedData,LoggedData2,index_Time,index_RPM,index_TP, ...
    index_FEPW,index_WSSFL,index_WSSRL,fs] = processLog(logFileName)
% Import Log File and Convert to Table and Double Arrays
% Table used for easy referencing, double used for data processing

% Below warning removes the following error caused by unideal header names
% in the csv files: "Warning: Table variable names were modified to make 
% them valid MATLAB identifiers. The original names are saved in the 
% VariableDescriptions property."
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
logfilename = logFileName;
LoggedData = readtable(logfilename); % Converts logged data to table
try
    LoggedData2 = readmatrix(logfilename); % Converts logged data to double array
catch
    fprintf(['Could not find function "readmatrix".\n' ...
        'Likely using MATLAB version pre-2019a.\n' ...
        'Current MATLAB version: %s\n' ...
        'Will proceed using preprocessed data log.\n'], version)
    if (strcmp(logfilename,'FSAEM_Endurance_20190511-1260803.csv'))
        preProcessedData = load('LD2FSAEM.mat');
        LoggedData2 = preProcessedData.LD2FSAEM;
        clear preProcessedData;
    elseif (strcmp(logfilename,'FSAEL_Endurance_20190622-1260800_MATLABfix.csv'))
        preProcessedData = load('LD2FSAEL.mat');
        LoggedData2 = preProcessedData.LD2FSAEL;
        clear preProcessedData;
    end
end
LoggedData2 = LoggedData2(1:end-1,:); % Get rid of last row which is NaN
% NOTE: LoggedData2 array will have samples indexed one row higher than
% LoggedData table due to second row of LoggedData table showing units
% All calculations below done with respect to index of LoggedData2 array
% Therefore, a straight corresponding to indexes 2200 to 2225 in
% LoggedData2 array corresponds to indexes 2201 to 2226 in LoggedData table
% NOTE2: Lincoln csv has a weird issue with the header in the csv and
% MATLAB interprets it as 4 rows as opposed to 2 rows with FSAEM data. As a
% result the LoggedData2 array will have samples indexed THREE rows higher
% than the LoggedData table

% Find Header Names of Logged Channels and Store Respective Column Numbers
Names = LoggedData.Properties.VariableNames; % Header names of logged data
cols_LoggedData = size(LoggedData2,2);
rows_LoggedData = size(LoggedData2,1);

% Find Column Number from Logged Data
index_Time = find(strcmpi(LoggedData.Properties.VariableNames,'Time'));
index_RPM = find(strcmpi(LoggedData.Properties.VariableNames,'EngineRPM'));
index_MAP = find(strcmpi(LoggedData.Properties.VariableNames,'ManifoldPres'));
index_TP = find(strcmpi(LoggedData.Properties.VariableNames,'ThrottlePos'));
index_FBPW = find(strcmpi(LoggedData.Properties.VariableNames,'FuelBasePW'));
index_FEPW = find(strcmpi(LoggedData.Properties.VariableNames,'FuelEffectivePW'));
index_Lambda = find(strcmpi(LoggedData.Properties.VariableNames,'Lambda1'));
index_LambdaAim = find(strcmpi(LoggedData.Properties.VariableNames,'Lambda1Aim'));
index_ET = find(strcmpi(LoggedData.Properties.VariableNames,'EngineTemp'));
index_FT = find(strcmpi(LoggedData.Properties.VariableNames,'FuelTemp'));
index_OT = find(strcmpi(LoggedData.Properties.VariableNames,'EngOilTemp'));
index_FP_abs = find(strcmpi(LoggedData.Properties.VariableNames,'FuelPressureSense'));
index_AE = find(strcmpi(LoggedData.Properties.VariableNames,'FuelAccel'));
index_DE = find(strcmpi(LoggedData.Properties.VariableNames,'FuelDecel'));
index_StartComp = find(strcmpi(LoggedData.Properties.VariableNames, 'FuelStartingComp'));
index_La1ST = find(strcmpi(LoggedData.Properties.VariableNames,'La1ShortTrim'));
index_La1LT = find(strcmpi(LoggedData.Properties.VariableNames,'La1LongTrim'));
index_WSSFL = find(strcmpi(LoggedData.Properties.VariableNames,'GroundSpeedLeft'));
index_WSSRL = find(strcmpi(LoggedData.Properties.VariableNames,'DriveSpeedLeft'));

% Adjust for incorrect WSS calbration
LoggedData2(:,index_WSSFL)= LoggedData2(:,index_WSSFL)*2361/2045.51;

% Find sampling frequency
fs = 1/(LoggedData2(2,index_Time)-LoggedData2(1,index_Time)); % [Hz]

% Filter RPM
originalRPM = LoggedData2(:,index_RPM);
LoggedData2(:,index_RPM) = channelFilt(originalRPM,12,0.4); 

% Filter WSSFL
originalWSSFL = LoggedData2(:,index_WSSFL);
LoggedData2(:,index_WSSFL) = channelFilt(originalWSSFL,12,0.14); 

% Filter WSSRL
originalWSSRL = LoggedData2(:,index_WSSRL);
LoggedData2(:,index_WSSRL) = channelFilt(originalWSSRL,12,0.14); 
