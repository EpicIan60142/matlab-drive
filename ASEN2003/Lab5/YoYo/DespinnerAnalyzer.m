%% YoYo Despinner Exploratory Exercise
% Group 27: Ian Faber, Aaron orshan, Owen Mcclung, and Chris Kong

clc; clear; close all;

% Derived radial length: 6.904 to 7.048 in

cases.case{1} = load('New Data\5pt5_inches');
cases.case{2} = load('New Data\6_inches');
cases.case{3} = load('New Data\6pt5_inches');
cases.case{4} = load('New Data\7_inches');
cases.case{5} = load('New Data\7pt5_inches');
cases.case{6} = load('New Data\8_inches');
cases.case{7} = load('New Data\8pt5_inches');
cases.case{8} = load('New Data\9_inches');
cases.friction = load('New Data\no_masses1');

closestIndex = 6;

if closestIndex <= 1
   closestIndex = 2;
elseif closestIndex >= 8
    closestIndex = 7;
end

closest.time = cases.case{closestIndex}(:,1);
closest.rpm = cases.case{closestIndex}(:,2);

below.time = cases.case{closestIndex-1}(:,1);
below.rpm = cases.case{closestIndex-1}(:,2);

above.time = cases.case{closestIndex+1}(:,1);
above.rpm = cases.case{closestIndex+1}(:,2);

friction.time = cases.friction(:,1);
friction.rpm = cases.friction(:,2);

figure
hold on;
title("Ideal Cord Length Comparison")
plot(closest.time/1000, closest.rpm)
plot(above.time/1000, above.rpm)
plot(below.time/1000, below.rpm)
xlabel("Time (sec)")
ylabel("Angular velocity (rpm)")
legend("Closest Length", "Above Length", "Below length")
hold off;

figure
hold on;
title("Effect of Friction on Spin Rate")
plot(friction.time/1000, friction.rpm)
xlabel("Time (sec)")
ylabel("Angular Velocity (rpm)")
hold off;


