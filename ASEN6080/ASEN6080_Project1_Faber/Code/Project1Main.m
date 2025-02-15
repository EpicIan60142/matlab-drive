%% ASEN 6080 Project 1 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
addpath(genpath("..\..\Utilities"))

earthConst = getPlanetConst();
scConst = getSCConst();
stations = makeStations();

%% Load and parse provided data
data = readmatrix("..\Data\project.txt");

stations = readData(stations, data);

%% Plot measurements
titleText = "Range and Range-Rate Measurements";
xLabel = "Time [sec]"; yLabel = ["Range [m]","Range-Rate [m/s]"];
plotMeasurements(stations, titleText, xLabel, yLabel);

%% Filter setup
stationX0 = [];
for k = 1:length(stations)
    stationX0 = [stationX0; stations(k).X0];
end

tspan = data(:,1); % Timespan of the data

X0 = [scConst.X0_cart; earthConst.mu; earthConst.J2; scConst.Cd; stationX0];
P0 = diag([
            1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e20, 1e6, 1e6, ...
            1e-10, 1e-10, 1e-10, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6
          ]);

opt = odeset('AbsTol',1e-12,'RelTol',1e-12);

%% Program setup
choice = input(sprintf(...
                        "Choose an option to run by typing in the number to the left of it:\n" + ...
                        "\t 1. Question 2/3 - All measurements, original P0\n" + ...
                        "\t 2. Question 5 - Range vs. range rate data strengths\n" + ...
                        "\t 3. Question 6 - No fixed stations\n" + ...
                        "\t 4. Question 6 - Fix station 337 instead of 101\n" + ...
                        "\t 5. Question 6 - Fix station 394 instead of 101\n" + ...
                        ... "\t 6. Animate the problem based on the results of questions 2/3!" + ...
                        "\t Press 'Ctrl-C' to exit\n" + ...
                        "\tChoice:\t" ...
                       ) ...
               );
switch choice
    case 1
        runQuestion = 1;
    case 2
        runQuestion = 2;
    case 3
        runQuestion = 3;
    case 4
        runQuestion = 4;
    case 5
        runQuestion = 5;
    case 6
        runQuestion = 6;
    otherwise
        fprintf("Invalid input, provide the number of the option you'd like to run!")
end

%% Run code depending on the desired situation
if runQuestion == 1
    fprintf("\n\nRunning Questions 2/3!\n\n")

    % Run Batch - all measurements, original P0
    measInclude = true(1,2);
    batchRun = runBatch(X0, zeros(size(X0)), P0, earthConst, scConst, stations, tspan, [], opt, 10, measInclude);
    
    fprintf("\nPress 'Enter' to run the LKF.\n")
    pause

    % Run LKF - original P0
    LKFRun = runLKF(X0, zeros(size(X0)), P0, earthConst, scConst, stations, 5);

    fprintf("\n\nDone, have a great day!\n\n")
elseif runQuestion == 2
    fprintf("\n\nRunning Question 5!\n\n")

    % Run Batch - Range only, original P0
    fprintf("\n\tRange data only:\n")
    measInclude = [true, false];
    batchRunRange = runBatch(X0, zeros(size(X0)), P0, earthConst, scConst, stations, tspan, [], opt, 10, measInclude);

    fprintf("\nPress 'Enter' to continue to range rate data only.\n")
    pause

    % Run Batch - Range rate only, original P0
    fprintf("\n\tRange rate data only:\n")
    measInclude = [false, true];
    batchRunRangeRate = runBatch(X0, zeros(size(X0)), P0, earthConst, scConst, stations, tspan, [], opt, 10, measInclude);

    fprintf("\n\nDone, have a great day!\n\n")
elseif runQuestion == 3
    fprintf("\n\nRunning Question 6 - no fixed stations!\n\n")

    % Remove fixed station
    P0 = diag([
                1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e20, 1e6, 1e6, ...
                1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6
              ]);

    % Run Batch - all measurements, new P0
    measInclude = true(1,2);
    batchRun = runBatch(X0, zeros(size(X0)), P0, earthConst, scConst, stations, tspan, [], opt, 10, measInclude);
    
    fprintf("\nPress 'Enter' to run the LKF.\n")
    pause

    % Run LKF - new P0
    LKFRun = runLKF(X0, zeros(size(X0)), P0, earthConst, scConst, stations, 5);

    fprintf("\n\nDone, have a great day!\n\n")
elseif runQuestion == 4
    % Set station 337 to be fixed
    P0 = diag([
                1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e20, 1e6, 1e6, ...
                1e6, 1e6, 1e6, 1e-10, 1e-10, 1e-10, 1e6, 1e6, 1e6
              ]);

    % Run Batch - all measurements, new P0
    measInclude = true(1,2);
    batchRun = runBatch(X0, zeros(size(X0)), P0, earthConst, scConst, stations, tspan, [], opt, 10, measInclude);
    
    fprintf("\nPress 'Enter' to run the LKF.\n")
    pause

    % Run LKF - new P0
    LKFRun = runLKF(X0, zeros(size(X0)), P0, earthConst, scConst, stations, 5);

    fprintf("\n\nDone, have a great day!\n\n")
elseif runQuestion == 5
    % Set station 394 to be fixed
    P0 = diag([
                1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e20, 1e6, 1e6, ...
                1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e-10, 1e-10, 1e-10
              ]);

    % Run Batch - all measurements, new P0
    measInclude = true(1,2);
    batchRun = runBatch(X0, zeros(size(X0)), P0, earthConst, scConst, stations, tspan, [], opt, 10, measInclude);
    
    fprintf("\nPress 'Enter' to run the LKF.\n")
    pause

    % Run LKF - new P0
    LKFRun = runLKF(X0, zeros(size(X0)), P0, earthConst, scConst, stations, 5);

    fprintf("\n\nDone, have a great day!\n\n")
end
