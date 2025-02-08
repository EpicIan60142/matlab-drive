%% Housekeeping
clc; clear; close all;

%% Constants (problem 4a)
Alat =  [
            -0.0561,     0,          -775.0000,      32.2000,   0,   0
            -0.0038,     -0.4138,    0.4231,         0,         0,   0
            0.0011,      -0.0064,    -0.1456,        0,         0,   0
            0,           1.0000,     0,              0,         0,   0
            0,           0,          1,              0,         0,   0             
            1,           0,          0,              0,         775, 0   
        ];

Blat = [
            0,          5.601
            -0.14120,   0.1201
            0.00375,    -0.4898
            0,          0
            0,          0
            0,          0
       ];

u0 = -Alat(1,3);

%% Problem 4b

eigVals = eig(Alat);
[eigVec, ~] = eig(Alat);

%% Problem 4d

kr = -2;

K = [
        0, 0, 0,  0, 0, 0
        0, 0, kr, 0, 0, 0
    ];

%% Problem 4e

Anew = Alat - Blat*K;

eigValsNew = eig(Anew);
[eigVecNew, ~] = eig(Anew);

%% Problem 4f

% Process modes
drModes = eigVecNew(:,3:4);
drModes(1,:) = drModes(1,:)/u0; % Convert v to beta
drModes(6,:) = drModes(6,:)/100; % Units of hundreds of feet
drModes = drModes./(drModes(5,:));
drModes = drModes*0.1;

realDr = real(drModes(:,1));
imagDr = imag(drModes(:,1));

phModes = eigVecNew(:,5:6);
phModes(1,:) = phModes(1,:)/u0; % Convert v to beta
phModes(6,:) = phModes(6,:)/100; % Units of hundreds of feet
phModes = phModes./(phModes(5,:));
phModes = phModes*0.1;

realPh = real(phModes(:,1));
imagPh = imag(phModes(:,1));

% Plot phasors

size = 25;
col = [[1 0 0]; [0 1 0]; [0 0 1]; [1 0 1]; [0 0 0]; [0 1 1]]; % One color for each mode component
variable = ["\beta", "p", "r", "\phi", "\psi", "y"]; % Vector of mode components

figure(9)
hold on
grid on
title("Dutch Roll Mode Phasor Plot")
xline(0, 'k--')
yline(0, 'k--')
scatter(realDr, imagDr, size, col, 'filled')
for k = 1:length(realDr)
    dutchPlot(k) = plot([0,realDr(k)], [0,imagDr(k)], 'Color', col(k,:));
    if k == 5
        label(k) = sprintf("$%s = 0.1$", variable(k));
%     elseif k == 1
%         label(k) = sprintf("$\\hat{%s}$", variable(k));
    elseif k == 6
        label(k) = sprintf("$\\Delta %s$", variable(k));
    else
        label(k) = sprintf("$%s$", variable(k));
    end
end

subset = dutchPlot;
legend(subset, label, 'Location', 'best', 'Interpreter', 'latex');

figure(10)
hold on
grid on
title("Lateral Phugoid Mode Phasor Plot")
xline(0, 'k--')
yline(0, 'k--')
scatter(realPh, imagPh, size, col, 'filled')
for k = 1:length(realPh)
    phugoidPlot(k) = plot([0,realPh(k)], [0,imagPh(k)], 'Color', col(k,:));
    if k == 5
        label(k) = sprintf("$%s = 0.1$", variable(k));
%     elseif k == 1
%         label(k) = sprintf("$\\hat{%s}$", variable(k));
    elseif k == 6
        label(k) = sprintf("$\\Delta %s$", variable(k));
    else
        label(k) = sprintf("$%s$", variable(k));
    end
end

subset = phugoidPlot;
legend(subset, label, 'Location', 'best', 'Interpreter', 'latex');

figure(11)
hold on
grid on
title("Modes Comparison Phasor Plot")
xline(0, 'k--')
yline(0, 'k--')
scatter(realDr, imagDr, size, col, 'filled')
for k = 1:length(realDr)
    dutchPlot2(k) = plot([0,realDr(k)], [0,imagDr(k)], 'Color', col(k,:));
    if k == 5
        label(k) = sprintf("$%s = 0.1$", variable(k));
%     elseif k == 1
%         label(k) = sprintf("$\\hat{%s}$", variable(k));
    elseif k == 6
        label(k) = sprintf("$\\Delta %s$", variable(k));
    else
        label(k) = sprintf("$%s$", variable(k));
    end
end

scatter(realPh, imagPh, size, col, 'filled')
for k = 1:length(realPh)
    phugoidPlot2(k) = plot([0,realPh(k)], [0,imagPh(k)], 'Color', col(k,:), 'LineStyle', '--');
    if k == 5
        label(k) = sprintf("$%s = 0.1$", variable(k));
%     elseif k == 1
%         label(k) = sprintf("$\\hat{%s}$", variable(k));
    elseif k == 6
        label(k) = sprintf("$\\Delta %s$", variable(k));
    else
        label(k) = sprintf("$%s$", variable(k));
    end
end

subset = [dutchPlot2(1), phugoidPlot2(1)];
legend(subset, ["Dutch Roll", "Lateral Phugoid"], 'Location', 'best');

