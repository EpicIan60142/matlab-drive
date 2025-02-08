%% OEMP 1 Problem 5: Max Point Load
% By: Ian Faber, 09/15/21

clc;clear;close all;

data = load('OEMP1_data.mat');
flight_envelope = data.flight_envelope;

Cd = 1.2;
A = 2*0.01859; %0.2 ft^2
v = flight_envelope.mach.*flight_envelope.('speed_of_sound (m/s)');
rho = flight_envelope.('air_density (kg/m^3)');

F = 0.5*Cd*A*rho.*(v.^2); %N

maxF = max(F)
index = find(F == maxF);
altitude = flight_envelope.('altitude (m)')(index)
mach = flight_envelope.mach(index)

figure

sgtitle('OEMP 1 data')

subplot(3,1,1)
hold on;
yyaxis left
plot(flight_envelope.('altitude (m)'),v,'r');
plot(flight_envelope.('altitude (m)'),flight_envelope.('speed_of_sound (m/s)'),'b');
ylabel('air speed (m/s)')
yyaxis right
plot(flight_envelope.('altitude (m)'),flight_envelope.mach,'b.')
ylabel('mach number')
subtitle('Altitude vs. air speed and mach number')
xlabel('Altitude (m)')
grid on;
grid minor;

subplot(3,1,2)
plot(flight_envelope.('altitude (m)'),flight_envelope.('air_density (kg/m^3)'));
subtitle('Altitude vs. air density')
xlabel('Altitude (m)')
ylabel('air density (kg/m^3)')
grid on;
grid minor;

subplot(3,1,3)
plot(flight_envelope.('altitude (m)'),F);
subtitle('Altitude vs. RAT Point Load')
xlabel('Altitude (m)')
ylabel('Point force (N)')
grid on;
grid minor;





