
clc;
clear all;
close all;


recuv_tempest;

Va_trim = 21;
h_trim = 1800;

wind_inertial = [0;0;0];

trim_definition = [Va_trim; h_trim];


%%% Use full minimization to determine trim
[trim_variables, fval] = CalculateTrimVariables(trim_definition, aircraft_parameters);
[trim_state, trim_input]= TrimStateAndInput(trim_variables, trim_definition);
[Alon, Blon, Alat, Blat] = AircraftLinearModel(trim_definition, trim_variables, aircraft_parameters);

[num, den] = ss2tf(Alat, Blat(:,2), [0 0 1 0 0 0], 0);

% Yaw Damper Root Locus
figure(1)
hold on
rlocus(num, den)
title("Positive k_r Root Locus")

figure(2)
hold on
rlocus(-num, den)
title("Negative k_r Root Locus")

eigValsLat = eig(Alat);
[eigVecLat, ~] = eig(Alat);

drPolesLat = eigValsLat(4);
[drWn, drDamp] = damp(drPolesLat);

kr = -6.83;
K_yaw = [0 0 kr 0 0 0];

Ayaw = (Alat - Blat(:,2)*K_yaw);

eigValsGain = eig(Ayaw);
[eigVecGain, ~] = eig(Ayaw);

drPolesGain = eigValsGain(4);
[drWnGain, drDampGain] = damp(drPolesGain);

% return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Transfer function from rudder to yaw rate
[num_rud2r, den_rud2r] = ss2tf(Alat, Blat(:,2), [0 0 1 0 0 0], 0);

%%% Controller
kr = -6.83;

%%% Closed loop transfer function
G_yaw = tf(num_rud2r, den_rud2r);
G_ol = tf(num, den);
C = tf(kr,1);

rud_noDamp = feedback(G_ol, 1);
rud_ol = feedback(G_yaw,1);
rud_cl = feedback(C*G_yaw, 1);

% Check eigenvalues
[rud_num, rud_den] = tfdata(rud_cl, 'v');
eigValsTF = roots(rud_den);

% Impulse rudder response
figure(3)
hold on
impulse(rud_noDamp)
impulse(rud_cl)
title("Yaw Damper Impulse Response")
ylabel("Yaw rate [rad/s]")
xlim([0 100])
legend("No Yaw Damper", "Yaw Damper")

% return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper code for Problem 3. Look inside the AircraftEOMControl function
% for hints on how to set up controllers for Problem 2 and Problem 4.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Transfer function from elevator to pitch angle
[num_elev2pitch, den_elev2pitch] = ss2tf(Alon, Blon(:,1), [0 0 0 1 0 0],0);

%%% Controller
kq = -1; %NOT REASONABLE VALUES
kth = -4.55; %NOT REASONABLE VALUES
num_c = [kq kth];
den_c = 1;

%%% Closed loop transfer function
G_pitch = tf(num_elev2pitch, den_elev2pitch);

pitch_cl = feedback(tf(conv(num_c, num_elev2pitch), conv(den_c, den_elev2pitch)),1);
% pitch_ol = feedback(G, 1);

[num_cl, den_cl] = tfdata(pitch_cl,'v');

%%% Poles of the closed loop (linear) system. Now do the same with the
%%% state stpace model.
eigValsLon = roots(den_elev2pitch)
eigValsPitch = roots(den_cl)

% State Space Recovery
K_pitch = [0 0 kq kth 0 0];

A_pitch = Alon - Blon(:,1)*K_pitch;

% eigValsPitchMat = eig(A_pitch)

[num_ail2p, den_ail2p] = ss2tf(Alat, Blat(:,1), [0 1 0 0 0 0], 0)

figure(4)
hold on
step(G_pitch)
step(pitch_cl)
title("Pitch Control Step Response")
ylabel("Pitch [rad]")
xlim([0 100])
legend("No Pitch Control", "Pitch Control")

%%% Controller inner
ka = -0.5;
C = tf(ka, 1);

G_roll_inner = tf(num_ail2p, den_ail2p);

roll_cl_inner = feedback(C*G_roll_inner,1)

[num_roll_cl, den_roll_cl] = tfdata(roll_cl_inner, 'v');

eigvalsRoll = roots(den_ail2p)
eigValsRollGain = roots(den_roll_cl)

figure(5)
hold on
impulse(G_roll_inner)
impulse(roll_cl_inner)
ylabel("Roll rate [rad/s]")
title("Roll Control Inner Impulse Response")
legend("No Roll Control", "Roll Control")

%%% Controller outer
[num_phi2p, den_phi2p] = ss2tf(Alat, Blat(:,1), [0 0 0 1 0 0], 0);

G_roll_outer = tf(num_phi2p, den_phi2p);

kp = 20;
C = tf(kp, 1);

roll_cl_outer = feedback(C*roll_cl_inner, 1, 1 ,1);

[num_roll_cl_outer, den_roll_cl_outer] = tfdata(roll_cl_outer, 'v');

eigValsRollOuter = roots(den_phi2p)
eigValsRollOuterGain = roots(den_roll_cl_outer)

figure(6)
hold on
step(G_roll_outer)
step(roll_cl_outer)
ylabel("Roll rate [rad/s]")
title("Roll Control Outer Step Response")
xlim([0 100])
legend("No Roll Control", "Roll Control")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return;

%%% Full sim in ode45

close all;

aircraft_state0 = trim_state;
% aircraft_state0(12) = 0.01; % Initial yaw rate perturbation for yaw
% damper
control_input0 = trim_input;

tfinal = 400;
TSPAN = [0 tfinal];

% No Yaw Damper
% [TOUT, YOUT] = ode45(@(t,y) AircraftEOM(t,y,control_input0,wind_inertial,aircraft_parameters),TSPAN,aircraft_state0,[]);

% Yaw Damper
[TOUT2,YOUT2] = ode45(@(t,y) AircraftEOMControl(t,y,control_input0,wind_inertial,aircraft_parameters),TSPAN,aircraft_state0,[]);

% for i=1:length(TOUT)
%     UOUT(i,:) = control_input0';
% end

for i=1:length(TOUT2)
    theta_c = 5*pi/180;

    r_c = 0;

    omega = 0.021;
    u0 = 21;
    theta = 0;
    g = 9.8062;
    phi_c = atan2(omega*u0,g*cos(theta));
    
    ail_perturb = rollControlFull(phi_c, YOUT2(4), YOUT2(10), ka, kp);

    elev_perturb = 0*PitchAttitudeControl(theta_c, YOUT2(i,5), YOUT2(i,11), kth, kq); 
    rud_perturb = YawDamperControl(0, YOUT2(i,12), kr);
    UOUT2(i,:) = control_input0 + [elev_perturb; ail_perturb; rud_perturb; 0];
end

% PlotAircraftSim(TOUT, YOUT, UOUT, wind_inertial, 'b');
traj = PlotAircraftSim(TOUT2,YOUT2,UOUT2,wind_inertial, 'b');

figure(traj.Number);
title("3D trajectory")
zlim([0 1801])
xlim([-1200 1200])
ylim([0 2400])

% legend("No Roll Control", "Roll Control")





