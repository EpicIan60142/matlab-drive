

clear all;
close all;


recuv_tempest;
ap = aircraft_parameters;
Va_trim = 22;
h_trim = 2438.5;

wind_inertial = [0;0;0];

trim_definition = [Va_trim; h_trim];


%%% Determine trim
[trim_variables, fval] = CalculateTrimVariables(trim_definition, aircraft_parameters);
[trim_state, trim_input]= TrimStateAndInput(trim_variables, trim_definition);


%%% Linear matrices
[Alon, Blon, Alat, Blat] = AircraftLinearModel(trim_definition, trim_variables, aircraft_parameters);

lat_eigen = eig(Alat);
[lat_eigVec, ~] = eig(Alat);
% Alat eigen, 1 and 2 correspond to Delta Yaw and Delta Y
% Alat eigen 3 is roll rate
% Alat eigen 4 and 5 are dutch roll
% Alat eigen 6 is spiral mode


lat_eigen_spiral_mode = lat_eigen(6);
lat_dr_wn = sqrt(lat_eigen(4)*lat_eigen(5));
lat_dr_zeta = -(lat_eigen(4)+lat_eigen(5))/(2*lat_dr_wn);
lat_spiral_tc = -1/lat_eigen(6);
lat_roll_tc = -1/lat_eigen(3);

dr_mode = lat_eigen(4:5);


drModes = lat_eigVec(1:5,4:5);
drModes(1,:) = drModes(1,:)./Va_trim;
drModes(2,:) = drModes(2,:)./(2*Va_trim/ap.b);
drModes(3,:) = drModes(3,:)./(2*Va_trim/ap.b);

drModes = drModes./(drModes(4,1:2));

realdr = real(drModes);
imagdr = imag(drModes);

size = 25;
col = [[1 0 0]; [0 1 0]; [0 0 1]; [1 0 1]; [0 0 0]];% [0 1 1]];
variable = ["v", "p", "r", "\phi", "\psi", "y^E"];

figure
hold on
grid on
title("Dutch Roll Mode Phasor Plot")
xline(0, 'k--')
yline(0, 'k--')
scatter(realdr(:,1), imagdr(:,1), size, col, 'filled')
for k = 1:length(drModes)
    dutchPlot(k) = plot([0,realdr(k)], [0,imagdr(k)], 'Color', col(k,:));
    if k == 4
        label(k) = sprintf("\\Delta %s = 1", variable(k));
    else
        label(k) = sprintf("\\Delta %s", variable(k));
    end
end
subsetlat = dutchPlot;

legend(subsetlat, label, 'Location', 'best')


return;

%%%%% Set initial condition
% STUDENTS COMPLETE

%%% Full sim in ode45
tfinal = 150;
TSPAN = [0 tfinal];
[TOUT2,YOUT2] = ode45(@(t,y) AircraftEOM(t,y,control_input0,wind_inertial,aircraft_parameters),TSPAN,aircraft_state0,[]);


for i=1:length(TOUT2)
    UOUT2(i,:) = control_input0';
end

PlotAircraftSim(TOUT,YOUT,UOUT,'b')


%%% Linear simulation
% STUDENTS COMPLETE


