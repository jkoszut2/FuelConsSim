clc
clear
close all

load Fuel26115e.mat
fuelTable = Fuel26115;

RPM = [0 500 1000 1250 1500 1630 1750 1880 2000 2250 2500 2750 3000:500:6500 6750 7000 7250 7500 7750 8000:500:15000];
MAP = [0.1:0.1:0.5 0.55:0.05:0.85 0.87 0.9 0.92 0.93 0.94 0.95 0.97 0.98 0.99 1];

MAP_WOT_2618e = [0.988 0.987 0.986 0.983 0.982 0.977 0.977 0.973 0.972 0.968 0.967 0.964 0.959 0.95 0.946 0.917 0.903 0.902 0.887 0.872 0.875 0.885 0.882 0.874 0.867 0.846 0.818];
MAP_WOT_2618p = [0.988 0.987 0.983 0.983 0.981 0.98 0.977 0.973 0.972 0.968 0.963 0.961 0.955 0.946 0.926 0.912 0.896 0.896 0.881 0.863 0.865 0.88 0.875 0.867 0.86 0.841 0.814];
MAP_WOT_2618 = 0.5*MAP_WOT_2618e + 0.5*MAP_WOT_2618p;

fuelSimRPM = 2000:500:14000;
FEPW = zeros(length(MAP_WOT_2618),1);
for i = 1:length(fuelSimRPM)
    currRPM = fuelSimRPM(i);
    currMAP = MAP_WOT_2618(i);
    tableRPMindex = find(RPM==currRPM);
    tableMAPindex = find(MAP>currMAP,1);
    tableMAPindex = tableMAPindex - 1;
    FEPWvals = fuelTable(tableMAPindex:end,tableRPMindex);
    FEPW(i) = interp1(MAP(tableMAPindex:end),FEPWvals,MAP_WOT_2618(i));
end

FEPW = FEPW/100*7;