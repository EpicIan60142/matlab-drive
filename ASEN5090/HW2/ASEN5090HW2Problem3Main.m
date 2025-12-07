%% ASEN 5090 HW 2 Problem 3 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

addpath("HW2_Code\")
addpath("HW2_Data\")
addpath("..\utilities\")

%% Setup
    % User position
NIST_ECEF = [-1288398.567; -4721696.932; 4078625.35]; % m

    % SP3 data
sp3Data = read_sp3("IGS0OPSFIN_20252230000_01D_15M_ORB.SP3");

    % Example GPS satellites to make measurements to
PRNs = [5, 16];

%% Part a/b. Write a range-azimuth-elevation calculation function and
%            calculate measurements for the same satellite as 1b both when 
%            visible and not, and repeat for a second satellite with 
%            different visibility

    % Pull out satellite data of interest
idxSats = sp3Data(:,3) == PRNs(1);
if length(PRNs) > 1
    for k = 2:length(PRNs)
        idxSats = idxSats + (sp3Data(:,3) == PRNs(k));
    end
end
idxSats = logical(idxSats);

t_Meas = unique(sp3Data(idxSats, 2));
t_Meas = (t_Meas - t_Meas(1))/3600;
satECEF = sp3Data(idxSats, 4:6)'*1000;

    % Call az-el-range function
[satAz, satEl, satRange] = compute_azelrange(NIST_ECEF, satECEF);

    % Pull out satellite measurements into individual PRN vectors
satMeasurements.Az = [];
satMeasurements.El = [];
satMeasurements.Range = [];
satMeasurements.visible = [];
for k = 1:length(PRNs)
    satMeasurements.Az = [satMeasurements.Az; satAz(k:length(PRNs):length(satAz))];
    satMeasurements.El = [satMeasurements.El; satEl(k:length(PRNs):length(satEl))];
    satMeasurements.Range = [satMeasurements.Range; satRange(k:length(PRNs):length(satRange))];
    satMeasurements.visible = [satMeasurements.visible; satEl(k:length(PRNs):length(satEl)) > 0];
end

for k = 1:length(PRNs)
    for kk = 1:length(satMeasurements.visible(k,:))
        if satMeasurements.visible(k, kk) == 0
            satMeasurements.visible(k,kk) = NaN;
        end
    end
end

%% Part c/d. Plot measurements for each satellite, both ignoring and taking
% into account visibility

    % Plot ignoring visibility
fig = figure; tl = tiledlayout(3,1); ax = [];

titleText = "Measurements from PRN";
for k = 1:length(PRNs)-1
    if length(PRNs) == 2
        titleText = titleText + sprintf(" %.f", PRNs(k));
    else
        titleText = titleText + sprintf(" %.f,", PRNs(k));
    end
end
titleText = titleText + sprintf(" and %.f, ignoring visibility", PRNs(end));
title(tl, titleText)

nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    for k = 1:length(PRNs)
        plot(t_Meas, satMeasurements.Az(k, :), '.', 'Displayname', sprintf("PRN %.f", PRNs(k)));
    end
    xlabel("Time [hr]"); ylabel("Azimuth angle [deg]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    for k = 1:length(PRNs)
        plot(t_Meas, satMeasurements.El(k, :), '.', 'Displayname', sprintf("PRN %.f", PRNs(k)));
    end
    xlabel("Time [hr]"); ylabel("Elevation angle [deg]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    for k = 1:length(PRNs)
        plot(t_Meas, satMeasurements.Range(k, :), '.', 'Displayname', sprintf("PRN %.f", PRNs(k)));
    end
    xlabel("Time [hr]"); ylabel("Range [m]")
    legend('Location', 'bestoutside')
linkaxes(ax, 'x');

    % Plot taking visibility into account
fig = figure; tl = tiledlayout(4,1); ax = [];

titleText = "Measurements from PRN";
for k = 1:length(PRNs)-1
    if length(PRNs) == 2
        titleText = titleText + sprintf(" %.f", PRNs(k));
    else
        titleText = titleText + sprintf(" %.f,", PRNs(k));
    end
end
titleText = titleText + sprintf(" and %.f, obeying visibility", PRNs(end));
title(tl, titleText)

nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    for k = 1:length(PRNs)
        plot(t_Meas, satMeasurements.Az(k,:).*satMeasurements.visible(k,:), '.', 'Displayname', sprintf("PRN %.f", PRNs(k)));
    end
    xlabel("Time [hr]"); ylabel("Azimuth angle [deg]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    for k = 1:length(PRNs)
        plot(t_Meas, satMeasurements.El(k,:).*satMeasurements.visible(k,:), '.', 'Displayname', sprintf("PRN %.f", PRNs(k)));
    end
    xlabel("Time [hr]"); ylabel("Elevation angle [deg]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    for k = 1:length(PRNs)
        plot(t_Meas, satMeasurements.Range(k,:).*satMeasurements.visible(k,:), '.', 'Displayname', sprintf("PRN %.f", PRNs(k)));
    end
    xlabel("Time [hr]"); ylabel("Range [m]")
    legend('Location', 'bestoutside')
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    for k = 1:length(PRNs)
        plot(t_Meas, satMeasurements.visible(k,:), '.', 'Displayname', sprintf("PRN %.f", PRNs(k)));
    end
    xlabel("Time [hr]"); ylabel("Visibility")
linkaxes(ax, 'x');


