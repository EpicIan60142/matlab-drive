
clc;
clear all;
close all;


recuv_tempest;

Va_trim = 22;
h_trim = 2438.5;

wind_inertial = [0;0;0];

trim_definition = [Va_trim; h_trim];


%%% Determine trim
[trim_variables, fval] = CalculateTrimVariables(trim_definition, aircraft_parameters);
[trim_state, trim_input]= TrimStateAndInput(trim_variables, trim_definition);


%%% Linear matrices
[Alon, Blon, Alat, Blat] = AircraftLinearModel(trim_definition, trim_variables, aircraft_parameters);

% Pull out eigenvalues and eigenvectors
lon_eigVal = eig(Alon);
[lon_eigVec, ~] = eig(Alon);

lat_eigVal = eig(Alat);
[lat_eigVec, ~] = eig(Alat);

% Process short period mode
    % Isolate short period poles and mode shapes
    % - Mode shapes organized [u; w; q; theta; x; z]
spPoles = lon_eigVal(3:4);
spModes = lon_eigVec(:,3:4);
    % Preserve dimensionalized mode shapes for sim initialization
spModes_dim = spModes;
    % Non-dimensionalize u, w, q for phasor plots
spModes(1,:) = spModes(1,:)/Va_trim;
spModes(2,:) = spModes(2,:)/Va_trim;
spModes(3,:) = spModes(3,:)/(2*Va_trim/aircraft_parameters.c);
    % Normalize by theta
spModes = spModes./spModes(4,1:2);
spModes_dim = spModes_dim./spModes_dim(4,1:2);
    % Pull out real and imaginary components for phasor plots
realSp = real(spModes);
imagSp = imag(spModes);
    % Pull out dimensionalized mode shape real components for sim initialization
realSp_dim = real(spModes_dim);
    % Calculate natural frequency and damping
spNatFreq = norm(spPoles(1));
[~, spDamp] = damp(spPoles(1));


% Process phugoid mode
    % Isolate phugoid poles and mode shapes
    % - Mode shapes organized [u; w; q; theta; x; z]
phPoles = lon_eigVal(5:6);
phModes = lon_eigVec(:,5:6);
    % Preserve dimensionalized mode shapes for sim initialization
phModes_dim = phModes;
    % Non-dimensionalize u, w, q for phasor plots
phModes(1,:) = phModes(1,:)/Va_trim;
phModes(2,:) = phModes(2,:)/Va_trim;
phModes(3,:) = phModes(3,:)/(2*Va_trim/aircraft_parameters.c);
    % Normalize by theta
phModes = phModes./phModes(4,1:2);
phModes_dim = phModes_dim./phModes_dim(4,1:2);
    % Pull out real and imaginary components for phasor plots
realPh = real(phModes);
imagPh = imag(phModes);
    % Pull out dimensionalized mode shape real components for sim initialization
realPh_dim = real(phModes_dim);
    % Calculate natural frequency and damping
phNatFreq = norm(phPoles(1));
[~, phDamp] = damp(phPoles(1));


% Process dutch roll mode
    % Isolate dutch roll poles and mode shapes
    % - Mode shapes organized [v; p; r; phi; psi; y]
drPoles = lat_eigVal(4:5);
drModes = lat_eigVec(:,4:5);
    % Preserve dimensionalized mode shapes for sim initialization
drModes_dim = drModes;
    % Non-dimensionalize v, p, r for phasor plots
drModes(1,:) = drModes(1,:)/Va_trim;
% drModes(2,:) = drModes(2,:)/(2*Va_trim/aircraft_parameters.b);
% drModes(3,:) = drModes(3,:)/(2*Va_trim/aircraft_parameters.b);
    % Normalize by phi
drModes = drModes./drModes(4,1:2);
drModes_dim = drModes_dim./drModes_dim(4,1:2);
     % Pull out real and imaginary components for phasor plots
realDr = real(drModes);
imagDr = imag(drModes);
    % Pull out dimensionalized mode shape real components for sim initialization
realDr_dim = real(drModes_dim);


% Process spiral mode
    % Isolate spiral poles and mode shapes
    % - Mode shapes organized [v; p; r; phi; psi; y]
sPoles = lat_eigVal(6);
sMode = lat_eigVec(:,6);

sMode_dim = sMode;

sMode_dim = sMode_dim./sMode_dim(4);

realS_dim = real(sMode_dim);


% Process roll mode
    % Isolate roll poles and mode shapes
    % - Mode shapes organized [v; p; r; phi; psi; y]
rPoles = lat_eigVal(3);
rMode = lat_eigVec(:,3);

rMode_dim = rMode;

rMode_dim = rMode_dim./rMode_dim(4);

realR_dim = real(rMode_dim);


%% Phasor plots

% Setup
size = 25;  % Size of points
col = [[1 0 0]; [0 1 0]; [0 0 1]; [1 0 1]];%; [0 0 0]; [0 1 1]]; % One color for each mode component
variable = ["u", "w", "q", "\theta", "x^E", "z^E"]; % Vector of mode components

figure(7)
hold on
grid on
title("Short Period Mode Phasor Plot")
xline(0, 'k--')
yline(0, 'k--')
scatter(realSp(1:4,1), imagSp(1:4,1), size, col, 'filled') % Plot real and imaginary components
for k = 1:4%length(spModes) % Problem only wants first 4 components
    shortPlot(k) = plot([0,realSp(k,1)], [0,imagSp(k,1)], 'Color', col(k,:)); % Draw line to each point
    if k == 4
        label(k) = sprintf("$\\Delta %s = 1$", variable(k));
    elseif k == 3
        label(k) = sprintf("$\\Delta %s$", variable(k));
    else
        label(k) = sprintf("$\\hat{%s}$", variable(k));
    end
end

subset = shortPlot; % Only create legend for the lines
legend(subset, label, 'Location', 'best', 'Interpreter', 'latex') % Make sure legend doesn't generate on top of phasor plot

figure(8)
hold on
grid on
title("Phugoid Mode Phasor Plot")
xline(0, 'k--')
yline(0, 'k--')
scatter(realPh(1:4,1), imagPh(1:4,1), size, col, 'filled')
for k = 1:4%length(phModes)
    phugoidPlot(k) = plot([0,realPh(k,1)], [0,imagPh(k,1)], 'Color', col(k,:));
    if k == 4
        label(k) = sprintf("$\\Delta %s = 1$", variable(k));
    elseif k == 3
        label(k) = sprintf("$\\Delta %s$", variable(k));
    else
        label(k) = sprintf("$\\hat{%s}$", variable(k));
    end
end

subset = phugoidPlot;
legend(subset, label, 'Location', 'best', 'Interpreter', 'latex')

col = [[1 0 0]; [0 1 0]; [0 0 1]; [1 0 1]; [0 0 0]];% [0 1 1]]; % One color for each mode component
variable = ["v", "p", "r", "\phi", "\psi", "y^E"]; % Vector of mode components

figure(9)
hold on
grid on
title("Dutch Roll Mode Phasor Plot")
xline(0, 'k--')
yline(0, 'k--')
scatter(realDr(1:5,1), imagDr(1:5,1), size, col, 'filled')
for k = 1:5%length(phModes)
    dutchPlot(k) = plot([0,realDr(k,1)], [0,imagDr(k,1)], 'Color', col(k,:));
    if k == 4
        label(k) = sprintf("$\\Delta %s = 1$", variable(k));
    elseif k == 1
        label(k) = sprintf("$\\hat{%s}$", variable(k));
    else
        label(k) = sprintf("$\\Delta %s$", variable(k));
    end
end

subset = dutchPlot;
legend(subset, label, 'Location', 'best', 'Interpreter', 'latex');

return;
% close all;

%%%%% Set initial condition
% STUDENTS COMPLETE
% Setup
deltaTheta = 0; % deg, pitch perturbation outlined in the problem statement
deltaPhi = 2; % deg, roll perturbation outlined in the problem statement
mode = realR_dim(:,1); % mode to simulate: trim (zeros(6,1)), short period (realSp_dim(:,1)), phugoid (realPh_dim(:,1)), dutch roll (realDr_dim(:,1)), spiral (realS_dim), roll (realR_dim)
tfinal = 150; % Change to value specified in problem

% Pull out trim variables
alpha = trim_variables(1);
del_e = trim_variables(2);
del_t = trim_variables(3);

% Calculate trim (can actually pull out of "trim_state" calculated on line
% 19, didn't see that until now -_-)
x_trim = [0; 0; -h_trim];
a_trim = [0; alpha; 0];
v_trim = [Va_trim*cos(alpha); 0; Va_trim*sin(alpha)];
w_trim = [0; 0; 0];
aircraft_state_trim = [x_trim; a_trim; v_trim; w_trim];

% Define the initial perturbation, corresponds to the real components of
% each mode shape
deltaX_lon = deg2rad(deltaTheta)*mode;
deltaX_lat = deg2rad(deltaPhi)*mode;

% Define initial state and control vectors, accounting for trim where necessary
x0 = deltaX_lon(5) + x_trim(1);
y0 = deltaX_lat(6) + x_trim(2);
z0 = deltaX_lon(6) + x_trim(3);
phi0 = deltaX_lat(4);
theta0 = deltaX_lon(4) + a_trim(2);
psi0 = deltaX_lat(5);
u0 = deltaX_lon(1) + v_trim(1);
v0 = deltaX_lat(1);
w0 = deltaX_lon(2) + v_trim(3);
p0 = deltaX_lat(2);
q0 = deltaX_lon(3);
r0 = deltaX_lat(3);

aircraft_state_0 = [x0; y0; z0; phi0; theta0; psi0; u0; v0; w0; p0; q0; r0];

control_input0 = [del_e; 0; 0; del_t];

%%% Full sim in ode45
% Nonlinear simulation
TSPAN = [0 tfinal];
[TOUT,YOUT] = ode45(@(t,y) AircraftEOM(t,y,control_input0,wind_inertial,aircraft_parameters),TSPAN,aircraft_state_0,[]);


for i=1:length(TOUT)
    UOUT(i,:) = control_input0';
end

% Display response
PlotAircraftSim(TOUT,YOUT,UOUT,'b')


%%% Linear simulation
% STUDENTS COMPLETE

% Specify a consant time vector so linear outputs are the same size
time = 0:0.01:tfinal;

trimLon = [v_trim(1), v_trim(3), w_trim(2), a_trim(2), x_trim(1), x_trim(3)]; % u, w, q, theta, x, z
sysLon = ss(Alon, Blon, eye(6,6), []); % Create state space model based on Alon
% X0Lon = zeros(1,6); % delta: u, w, q, theta, x, z
X0Lon = deg2rad(deltaTheta)*mode; % Specify initial linear perturbation (same as "deltaX_lon" above)
[deltXLon, linTLon, ~] = initial(sysLon, X0Lon, time); % Run linear simulation
linRespLon = trimLon + deltXLon; % Add perturbation output to trim

trimLat = [v_trim(2), w_trim(1), w_trim(3), a_trim(1), a_trim(3), x_trim(2)]; % v, p, r, phi, psi, y
sysLat = ss(Alat, Blat, eye(6,6), []);
% X0Lat = zeros(1,6); % delta: v, p, r, phi, psi, y
X0Lat = deg2rad(deltaPhi)*mode;
[deltXLat, linTLat, ~] = initial(sysLat, X0Lat, time);
linRespLat = trimLat + deltXLat;

% Parse linear response, linear time, and linear control vectors from longitudinal/lateral mode simulations
x = linRespLon(:,5) + Va_trim*time'; % Account for lack of xdot from "initial" output, otherwise plane doesn't move forward
y = linRespLat(:,6);
z = linRespLon(:,6);
phi = linRespLat(:,4);
theta = linRespLon(:,4);
psi = linRespLat(:,5);
u = linRespLon(:,1);
v = linRespLat(:,1);
w = linRespLon(:,2);
p = linRespLat(:,2);
q = linRespLon(:,3);
r = linRespLat(:,3);

linX = [x, y, z, phi, theta, psi, u, v, w, p, q, r];

linT = time';

for i=1:length(linT)
    linU(i,:) = control_input0';
end

% Display response
PlotAircraftSim(linT, linX, linU, 'r')

