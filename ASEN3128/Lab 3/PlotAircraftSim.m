function PlotAircraftSim(time, aircraft_state_array, control_input_array, fig, col)
% Section 011 - Gabriel Agostine, Felix Evrard, Ian Faber, Blake Hogen
%   Function to plot aircraft state and control inputs with given
%   formatting parameters
% Inputs:
%           time: Length n vector spanning simulation interval
%           aircraft_state_array: 12xn vector array of aircraft states
%           control_input_array: 4xn vector array of control inputs ([Zc,
%                                Lc, Mc, Nc]')
%           fig: 6x1 vector of figure numbers to plot over
%           col: String defining format for each plot (e.g. 'b-')

x = aircraft_state_array(1,:);
y = aircraft_state_array(2,:);
z = aircraft_state_array(3,:);

phi = aircraft_state_array(4,:);
theta = aircraft_state_array(5,:);
psi = aircraft_state_array(6,:);

u = aircraft_state_array(7,:);
v = aircraft_state_array(8,:);
w = aircraft_state_array(9,:);

p = aircraft_state_array(10,:);
q = aircraft_state_array(11,:);
r = aircraft_state_array(12,:);

Zc = control_input_array(1,:);
Lc = control_input_array(2,:);
Mc = control_input_array(3,:);
Nc = control_input_array(4,:);

% Inertial position
figure(fig(1))
sgtitle("Inertial Position of Drone")
subplot(3,1,1)
title("X position vs. time")
plot(time, x, col); hold on;
ylabel("X position (m)");
xlabel("time (sec)")

subplot(3,1,2)
title("Y position vs. time")
plot(time, y, col); hold on;
ylabel("Y position (m)");
xlabel("time (sec)")

subplot(3,1,3)
title("Z position vs. time")
plot(time, z, col); hold on;
ylabel("Z position (m)");
xlabel("time (sec)")

% Euler angles
figure(fig(2))
sgtitle("Euler Angles of Drone")
subplot(3,1,1)
title("\phi vs. time")
plot(time, phi, col); hold on;
ylabel("\phi (rad)");
xlabel("time (sec)")

subplot(3,1,2)
title("\theta vs. time")
plot(time, theta, col); hold on;
ylabel("\theta (rad)");
xlabel("time (sec)")

subplot(3,1,3)
title("\psi vs. time")
plot(time, psi, col); hold on;
ylabel("\psi (rad)");
xlabel("time (sec)")

% Inertial velocity
figure(fig(3))
sgtitle("Inertial Velocity of Drone")
subplot(3,1,1)
title("u^E vs. time")
plot(time, u, col); hold on;
ylabel("u^E (m/s)");
xlabel("time (sec)")

subplot(3,1,2)
title("v^E vs. time")
plot(time, v, col); hold on;
ylabel("v^E (m/s)");
xlabel("time (sec)")

subplot(3,1,3)
title("w^E vs. time")
plot(time, w, col); hold on;
ylabel("w^E (m/s)");
xlabel("time (sec)")

% Angular rates
figure(fig(4))
sgtitle("Angular Rates of Drone")
subplot(3,1,1)
title("p vs. time")
plot(time, p, col); hold on;
ylabel("p (rad/s)")
xlabel("time (sec)")

subplot(3,1,2)
title("q vs. time")
plot(time, q, col); hold on;
ylabel("q (rad/s)")
xlabel("time (sec)")

subplot(3,1,3)
title("r vs. time")
plot(time, r, col); hold on;
ylabel("r (rad/s)")
xlabel("time (sec)")

% Control inputs
figure(fig(5))
sgtitle("Control Inputs of Drone")
subplot(4,1,1)
title("Zc vs. time")
plot(time, Zc, col); hold on;
ylabel("Zc (N)")
xlabel("time (sec)")

subplot(4,1,2)
title("Lc vs. time")
plot(time, Lc, col); hold on;
ylabel("Lc (Nm)")
xlabel("time (sec)")

subplot(4,1,3)
title("Mc vs. time")
plot(time, Mc, col); hold on;
ylabel("Mc (Nm)")
xlabel("time (sec)")

subplot(4,1,4)
title("Nc vs. time")
plot(time, Nc, col); hold on;
ylabel("Nc (Nm)")
xlabel("time (sec)")

% 3D trajectory
figure(fig(6))
title("3D Trajectory of Drone")
plot3(x, y, z, col); hold on;
start = plot3(x(1), y(1), z(1), 'g.', 'markerSize', 15);
stop = plot3(x(end), y(end), z(end), 'r.', 'MarkerSize', 15);
set(gca, 'YDir', 'reverse', 'ZDir', 'reverse')
view([-70 40]);
xlabel("X position (m)");
ylabel("Y position (m)");
zlabel("Z position (m)");

subset = [start, stop];
labels = ["Trajectory start", "Trajectory end"];

legend(subset, labels)

end