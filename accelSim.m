function [accel_time,stateData] = ...
    accelSim(input_Dist_meters,input_Vel_meterspersec, ...
              input_RearVel_meterspersec,input_Gear, ...
              input_redline,input_fdr,input_FEPW,input_Torque_lbft, ...
              modifierArray,producePlots)
          
% Simulation for straight line non-traction limited acceleration

% DESCRIPTION:
% This tool is made to be used for conditions under which the following
% assumptions can be confidently made: 1. No rear wheel slip, 2. Wide-open
% throttle (ie. throttle position >= ~95%), 3. An accurate engine torque
% map is available, 4. An accurate engine fuel map is available
% This tool should NOT be used for conditions under which the following is
% true: 1. Rear wheel slip is present 2. Throttle cannot be approximated to
% be at 100%

% DIRECTIONS:
% Required inputs are: 1. Desired distance to travel, 2. initial velocity, 
% 3. initial rear wheel velocity, 4. initial gear, 5. redline, 6. fdr, 
% 7. 100% throttle fuel map, 8. 100% throttle torque map,
% and lastly two optional inputs:
% 9. A modifier array that can be used for running sweeps on parameters. If
% this feature is not desired, simply use a "0" in the function call.
% 10. A boolean to turn on plots which show various
% characteristics of the simulated acceleration (Torque, RPM, Tractive
% Force Fx, Fuel Consumption, Gear, etc.)
% **** In order to use this tool as standalone (e.g. not from within the
% fuelSim simulation), take the following steps.
% 1. Copy and paste a fuel map and torque map from the fuelSim.m file into
% your workspace and name them input_FEPW and input_Torque_lbft,
% respectively.
% 2. For whatever value is passed in as "input_Vel_meterspersec", make sure
% the same exact value is passed into input_RearVel_meterspersec
% 3. After taking the above steps, the function can be called. For example,
% accelSim(65, 10, 10, 1, 12500, 36/11, input_FEPW, input_Torque_lbft, 0, 1)
% runs accelSim for a distance of 65 meters with a starting vehicle
% velocity of 10 meters/sec. The car starts in 1st gear, has a redline of
% 12500, a final drive ratio of 36/11. The fuel and torque maps are pulled
% from the local workspace variables "input_FEPW" and "input_Torque_lbft".
% Lastly, the modifierArray feature is ignored, but plots are turned on.

if ~exist('producePlots','var') || isempty(producePlots)
    producePlots = 0;
end

if ~exist('modifierArray','var') || isempty(modifierArray)
    modifierArray = 0;
end

% Vehicle Parameters
global dt dteff fdr g gears m r redline shifttime shifttimer ...
    RPM FEPW Torque_Nm LC_wss wdist wb cgh cp uk
dt = 0.001; % secs
m = 450+175; % lbs
m = m/2.20462; % kg
wdist = 0.42; % weight distrubution F/R
r = 9.875*0.0254; % meters
wb = 5.325*12*2.54/100; % meters
cgh = .655833333*12*2.54/100; % meters
uk = 1.4; % tire coefficient of friction
cp = 0.45; % aerodynamic downforce center of pressure
g = 9.81; % m/s^2
redline = input_redline;
gears = 2.073*[31/12 32/16 30/18 26/18 27/21 23/20]; % crank:sprocket
fdr = input_fdr; % ratio
RPM = [0 2000:500:13500 14000]; % rpm
FEPW = input_FEPW;
Torque_lbft = input_Torque_lbft;
Torque_Nm = Torque_lbft/2.20462*g*12*2.54/100; % Nm
Lambda = [1 1.046 0.95 0.9667 0.982 0.9314 0.9608	0.9296 0.9347 0.9348 ...
    0.9357 0.9305 0.9203 0.9215 0.9274 0.919 0.9171	0.9181 0.9003 0.8841 ...
    0.9235 0.9071 0.8763 0.8519 0.8711 0.8509]; % Reference lambda
FEPW = FEPW./((Lambda+modifierArray(1))./Lambda);
Torque_Nm = Torque_Nm*(-2.65*modifierArray(1)^2 - 0.0437*modifierArray(1) + 1);
% LC_wss = [0 1 2 3 4 5 7.5 10 15 20 25];
shifttime = 0.125; % sec
shifttime = 0.125; % sec - SET TO 0 ONLY FOR RFS4
shifttimer = 0; % sec
shiftflag = 0; % boolean to signify paddle upshift
dteff = 0.85; % drivetrain efficiency
if redline>max(RPM)
    fprintf(">> Error! Engine data is not sufficiently populated!\n")
    fprintf(">> Script may not run correctly if time step is too large.\n")
end

% NOTE: No warning if rpm input to findTorque > redline
% Result for above is 0 and NOT nan, thus need another case to detect

distMax = input_Dist_meters;

shifttimer = 0; % ensure timer is at 0 in case params were not run

% This is used to preallocate a time array assuming the acceleration event
% will take no longer than tMax seconds
% Once the acceleration run has been simulated, only the nonzero data
% entries are used for the function's output
tMax = 6; % seconds

% Initialize all states and set their indices
% NOTE: If adding a new state, it must also be added to the header variable
% around line 340
% Time
index_Time = 1;
init_Time = 0; % seconds
% Position
index_Pos = 2;
init_Pos = 0; % meters
% Velocity
index_Vel = 3;
init_Vel = input_Vel_meterspersec; % meters/sec
% Gear
index_Gear = 5;
init_Gear = input_Gear; % gear 
% Wheel Speed (Front)
index_WSSFL = 7;
init_WSSFL = init_Vel*3600/1609.4; % miles/hr
% Wheel Speed (Rear)
index_WSSRL = 8;
init_WSSRL = input_RearVel_meterspersec*3600/1609.4; % miles/hr
% Drag
index_Drag = 17;
init_Drag = findDrag(init_WSSFL);
% Downforce
index_Downforce = 18;
init_Downforce = findDownforce(init_WSSFL);
% RPM
index_RPM = 6;
% init_RPM = init_Vel/(pi*2*r)*60*fdr*gears(init_Gear); % rpm
init_RPM = input_RearVel_meterspersec/(pi*2*r)*60*fdr*gears(init_Gear); % rpm
% Acceleration
index_Acc = 4;
init_Acc = ((findTorque(init_RPM)*fdr*gears(init_Gear)) ...
    /r-init_Drag)/m/g*dteff; % g's
% Instantatneous Fuel Consumption
index_FuelCons = 9;
init_FuelCons = 0; % cc
% Cumulative Fuel Consumption
index_CumFuelCons = 10;
init_CumFuelCons = 0; % cc
% Engine Torque
index_Torque_Nm = 11;
init_Torque_Nm = findTorque(init_RPM); % Nm
% Wheel Torque
index_Wheel_Torque_Nm = 12;
init_Wheel_Torque_Nm = init_Torque_Nm; % Nm
% Slip
index_Slip = 13;
init_Slip = 0; % boolean
% Weight transfer due to longitudinal accel
index_WeightTransferAcc = 14;
init_WeightTransferAcc = m*init_Acc*(cgh/wb)*g; % N
% Fx of rear tires
index_Fx = 15;
init_Fx = init_Torque_Nm*fdr*gears(init_Gear)*dteff/r; % N
% Max Fx
index_MaxFx = 16;
init_MaxFx = ((1-wdist)*m*9.8 + init_Downforce*cp ...
    + init_WeightTransferAcc)*uk; % N
% FEPW
index_FEPW = 19;
init_FEPW = interp1(RPM,FEPW,init_RPM);
% Injector Duty Cycle
index_IDC = 20;
init_IDC = init_FEPW/(1/(init_RPM/60/2)*1000)*100;

tmp1 = who('index*');
tmp3 = 0;
for i = 1:length(who('index*'))
    tmp2 = eval(tmp1{i});
    tmp3 = [tmp3 tmp2];
end
x = [zeros(tMax/dt,max(tmp3))];
clear tmp*

A = [0 1;
     0 0];
B = [0 0;
     1/m -1/m];

states = who('index*');
stateinits = who('init*');
for i = 1:length(who('index*'))
    x(1,eval(states{i})) = eval(stateinits{i});
end


t = 0:dt:tMax;
i = 2; % Start calculations the next sample after t=0 sec
while ((x(i-1,index_Pos) < distMax) && i<=(tMax/dt+1))
%     if producePlots
%         fprintf('Time: %0.5f\n', t(i))
%     end
    % Using previous solution, find wheel speed, gear, and max Fx
    prev_WSSFL = x(i-1,index_WSSFL); % front wheel speed [mph]
    prev_WSSRL = x(i-1,index_WSSRL); % rear wheel speed [mph]
    prev_Gear = findGear(prev_WSSRL, x(i-1,index_Gear));
    prev_Drag = findDrag(prev_WSSFL); % N
    prev_Downforce = findDownforce(prev_WSSFL); % N
    prev_WeightTransferAcc = m*x(i-1,index_Acc)*(cgh/wb)*g; % N
    prev_WeightRear = (1-wdist)*m*9.8 + prev_Downforce*cp ...
        + prev_WeightTransferAcc; % N
    maxFx = prev_WeightRear*uk; % N

    % Find engine speed and torque output
    prev_rpm = findRPM(prev_WSSRL,prev_Gear);
    % Find engine torque output
    prev_engine_Torque_Nm = findTorque(prev_rpm);
    prev_wheel_Torque_Nm = x(i-1,index_Wheel_Torque_Nm);

    % Find previous longitudinal force at rear wheels assuming no slip
    prev_Fx = prev_engine_Torque_Nm*fdr*gears(prev_Gear)*dteff/r;

    % Solve state space model
    u = [prev_Fx;
         prev_Drag];
    dxdt = A*x(i-1,index_Pos:index_Vel)' + B*u;
    x(i,index_Time) = t(i);
    x(i,index_Pos) = x(i-1,index_Pos) + dxdt(1)*dt;
    x(i,index_Vel) = x(i-1,index_Vel) + dxdt(2)*dt;
    x(i,index_Acc) = dxdt(2)/g;
    x(i,index_Gear) = prev_Gear;
    x(i,index_WSSFL) = x(i,index_Vel)*3600/1609.4;
    x(i,index_WSSRL) = x(i-1,index_WSSRL) + dxdt(2)*dt*3600/1609.4;
    curr_rpm = findRPM(x(i,index_WSSRL),prev_Gear);
    x(i,index_RPM) = curr_rpm;
    currFuelCons = findFuelCons(curr_rpm);
    x(i,index_FuelCons) = currFuelCons;
    x(i,index_CumFuelCons) = x(i-1,index_CumFuelCons) + currFuelCons;
    x(i,index_Torque_Nm) = prev_engine_Torque_Nm;
    x(i,index_Wheel_Torque_Nm) = prev_wheel_Torque_Nm;
    x(i,index_WeightTransferAcc) = prev_WeightTransferAcc;
    x(i,index_Fx) = prev_Fx;
    x(i,index_MaxFx) = maxFx;
    x(i,index_Drag) = prev_Drag;
    x(i,index_Downforce) = prev_Downforce;
    x(i,index_FEPW) = interp1(RPM,FEPW,curr_rpm);
    x(i,index_IDC) = x(i,index_FEPW)/(1/(curr_rpm/60/2)*1000)*100;

    if prev_Fx > maxFx
        slip = 1;
    else
        slip = 0;
    end
    x(i,index_Slip) = slip;
    i = i+1;
end

lastindex = i-1; % Avoid plotting zeros located at end of state array

accel_time = t(find(x(:,index_Pos)>=distMax,1)); % seconds
% Should not get in here but double check anyways to simplify debugging in
% the event that tMax is not large enough
if isempty(accel_time)
error(['Simulation ended before target distance was covered. ', ...
        'Increase tMax to provide simulation with enough time to ', ...
        'cover distance set by function input.\n', ...
        'Current tMax = %0.0f sec\nDistance input = %0.3f meters', ...
        '\nDistance traveled in simulation = %0.3f meters', ...
        '\nInitial RPM = %0.3f', ...
        '\nInitial WSSFL = %0.3f'], ...
        tMax, distMax, x(lastindex,index_Pos),x(1,index_RPM),x(1,index_WSSFL))
end
if producePlots
    fprintf('Accel time: %0.1f msec\n', accel_time*1000)
    fprintf('Fuel consumed: %0.4f mcc\n', x(i-1,index_CumFuelCons)*1000)
end

if producePlots
    figure
    subplot(2,2,1)
    plot(t(1:lastindex),x(1:lastindex,index_Pos))
    hold on
    plot(t(1:lastindex),x(1:lastindex,index_WSSFL))
    hold on
    yline(75);
    hold on
    grid on
    legend('{Position [m]}', 'WSS\_front [mph]', '75m')
    xlim([0 t(lastindex)])

    subplot(2,2,2)
    plot(t(1:lastindex),x(1:lastindex,index_Acc))
    hold on
    yyaxis right
    plot(t(1:lastindex),x(1:lastindex,index_Torque_Nm))
    hold on
    plot(t(1:lastindex),x(1:lastindex,index_Wheel_Torque_Nm), 'k-')
    hold on
    grid on
    legend('Acceleration [g]', 'Engine Torque [Nm]', 'Wheel Torque [Nm]')
    xlim([0 t(lastindex)])

    subplot(2,2,3)
    plot(t(1:lastindex),x(1:lastindex,index_Fx), 'linewidth', 3)
    hold on
    plot(t(1:lastindex),x(1:lastindex,index_MaxFx), '--', 'linewidth', 2)
    hold on
    yyaxis right
    plot(t(1:lastindex),x(1:lastindex,index_Gear))
    hold on
    area(t(1:lastindex),x(1:lastindex,index_Slip))
    ylim([0 max(x(1:lastindex,index_Gear))+0.5])
    hold on
    grid on
    yyaxis left
    legend('Fx [N]', 'Max Fx [N]', 'Gear', 'Slip')
    xlim([0 t(lastindex)])

    subplot(2,2,4)
    plot(t(1:lastindex),x(1:lastindex,index_RPM))
    hold on
    yyaxis right
    plot(t(1:lastindex),x(1:lastindex,index_FuelCons))
    hold on
    plot(t(1:lastindex),x(1:lastindex,index_CumFuelCons))
    hold on
    grid on
    legend('RPM', 'Inst Fuel Cons [cc]', 'Fuel Cons [cc]')
    grid on
    xlim([0 t(lastindex)])
end
%{
% Time
index_Time = 1;
% Position
index_Pos = 2;
% Velocity
index_Vel = 3;
% Gear
index_Gear = 5;
% RPM
index_RPM = 6;
% Acceleration
index_Acc = 4;
% Wheel Speed
index_WSSFL = 7;
% Instantatneous Fuel Consumption
index_FuelCons = 8;
% Cumulative Fuel Consumption
index_CumFuelCons = 9;
% Engine Torque
index_Torque_Nm = 10;
% Wheel Torque
index_Wheel_Torque_Nm = 11;
% Slip
index_Slip = 12;
% Weight Transfer due to longitudinal accel
index_WeightTransferAcc = 13;
% Fx of rear tires
index_Fx = 14;
% Max Fx
index_MaxFx = 15;
% Drag
index_Drag = 16;
% Downforce
index_Downforce = 17;
%}

header = {'time' 'pos' 'vel' 'acc' 'gear' 'rpm' 'wssfl' 'wssrl', ...
     'fuelcons', 'cumfuelcons' 'Torque_Nm' 'wheel_Torque_Nm' 'slip' , ...
    'Weight_Transfer_Acc' 'Fx' 'MaxFx' 'Drag' 'Downforce' 'FEPW' 'IDC'};
stateData = cell2table(num2cell(x(1:lastindex,:)), 'VariableNames', header);

function fuelcons = findFuelCons(curr_rpm)
% global dt FEPW RPM redline
% Injector flow rate below from Ford injector characterization data
% Ford injector data provided in mass flow rate, so needs to be converted
% to volumetric flow rate. Also needs to be adjusted based off Bernoulli's
% principle as Ford testing was performed at 270 kPa
% Average of gauge fuel pressure from FSAEM and FSAEL is 44.03psi
% 44.03psi --> 303.57616 kPa
inj_flowrate = 2.58*60/0.738*sqrt(303.57616/270); % cc/min
% Calculating volume for now so don't need fuel density
% fuel_density = 0.73; % kg/L

if (curr_rpm < redline)
    curr_fepw = interp1(RPM,FEPW,curr_rpm);
else
    % Should only happen right before a shift
    % Engine still supplies fuel under shift so use constant extrapolation
    % to assume FEPW right above redline is equal to FEPW at redline
    curr_fepw = FEPW(end);
end
inj_max_ms = 1/(curr_rpm/60/2)*1000; % max pulsewidth per cycle [ms]
inj_duty = curr_fepw/inj_max_ms*100; % duty cycle in percent
fuelcons = inj_duty/100*dt*inj_flowrate/60*4; % cc
end   

function downforce = findDownforce(velocity)
% Fit obtained from 2019 FCA wind tunnel data
%{
x = [35, 45, 55]'; % speed [mph]
f = [324.6733759, 519.6826238, 774.942625]'; % downforce [N]
M=fit(x, f, 'poly2', 'upper', [inf, 0, 0], 'lower', [0, 0, 0]);
%}
downforce = 0.257206872977041*velocity^2;
end

function drag = findDrag(curr_WSSf_mph)
% Fit obtained from 2019 FCA wind tunnel data
%{
x = [35, 45, 55]'; % speed [mph]
f = [197.960913, 317.1329839, 469.0276004]'; % drag [N]
M=fit(x, f, 'poly2', 'upper', [inf, 0, 0], 'lower', [0, 0, 0]);
%}
drag = 0.15614997429366*curr_WSSf_mph^2;
end

function gear = findGear(currVehV, prevGear)
% currVehV is current vehicle velocity in miles per hour
% prevGear is previous gear (min=1, max = 6) and is used to determine if
% a shift should occur
% global fdr gears redline shiftflag r
gear = 0;
j = 6;
while ((j>=1) && (currVehV/60*1609.4/(2*pi*r)*fdr*gears(j) < redline))
    gear = j;
    j = j-1;
end
% Avoid downshifting
if (prevGear > gear)
    gear = prevGear;
end
% Check if upshift present
if (gear ~= prevGear)
    shiftflag = 1;
else
    shiftflag = 0;
end
end

function gear = findGear2(currVehV)
% Same as findGear but one minor change:
% 1. No shift flag
% Use this function to initialize the accel sim

% currVehV is current vehicle velocity in miles per hour
% global fdr gears redline r
gear = 0;
j = 6;
while ((j>=1) && (currVehV/60*1609.4/(2*pi*r)*fdr*gears(j) < redline))
    gear = j;
    j = j-1;
end
end

function rpm_launch = findLCrpm(curr_wheelspeed_mph)
% global LC_rpm LC_wss
rpm_launch = interp1(LC_wss,LC_rpm,curr_wheelspeed_mph);
end

function rpm = findRPM(currVehV, currGear)
% currVehV is current vehicle velocity in miles per hour
% currGear is current gear (min=1, max = 6)
% global r fdr gears
rpm = currVehV/60*1609.4/(2*pi*r)*fdr*gears(currGear);
end

function torque = findTorque(curr_rpm)
% global dt RPM Torque_Nm redline ...
%     shiftflag shifttime shifttimer
if (curr_rpm <= redline)
    torque = interp1(RPM,Torque_Nm,curr_rpm);
else
%     torque = 0;
    torque = interp1(RPM,Torque_Nm,redline);
end
% Check to make sure not shifting gears
if (shifttimer - dt > 0)
    shifttimer  = shifttimer - dt;
    torque = 0;
elseif (shifttimer - dt < 0)
        shifttimer = 0;
end
if (shiftflag)
    torque = 0;
    shiftflag = 0;
    shifttimer  = shifttimer + shifttime;
end
end

function VehV = findVehV(currRPM, currGear)
% Finds vehicle velocity in miles per hour
% global r fdr gears
VehV = currRPM*60/1609.4*(2*pi*r)/gears(currGear)/fdr;
end
end