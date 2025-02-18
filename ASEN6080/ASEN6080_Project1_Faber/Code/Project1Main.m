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
                        "\t 6. Extra - Potter Algorithm Comparison to LKF\n" + ...
                        "\t 7. Extra - Reasonable P0 Estimate\n" + ... 
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
    case 7
        runQuestion = 7;
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

    choice = input(sprintf("Animate solution? [y/n]\n\tChoice:\t"), "s");
    switch choice
        case 'y'
            titleText = "ASEN 6080 Project 1 Problems 2/3 - Batch";
            movieTitle = "../Videos/ASEN6080_Project1_P2-3_Batch";
            saveMovie = false;
            animateProblem(batchRun.X_batchFilt(2:end,:)', stations, earthConst, titleText, movieTitle, saveMovie)

            titleText = "ASEN 6080 Project 1 Problems 2/3 - LKF";
            movieTitle = "../Videos/ASEN6080_Project1_P2-3_LKF";
            saveMovie = false;
            animateProblem(LKFRun.X_LKF, stations, earthConst, titleText, movieTitle, saveMovie)

            fprintf("\n\nDone, have a great day!\n\n")
        case 'n'
            fprintf("\n\nDone, have a great day!\n\n")
        otherwise
            fprintf("\n\nDone, have a great day!\n\n")
    end
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

    choice = input(sprintf("Animate solution? [y/n]\n\tChoice:\t"), "s");
    switch choice
        case 'y'
            titleText = "ASEN 6080 Project 1 Problem 5 - Range Data Only";
            movieTitle = "../Videos/ASEN6080_Project1_P5_Range";
            saveMovie = false;
            animateProblem(batchRunRange.X_batchFilt(2:end,:)', stations, earthConst, titleText, movieTitle, saveMovie)

            titleText = "ASEN 6080 Project 1 Problem 5 - Range Rate Data Only";
            movieTitle = "../Videos/ASEN6080_Project1_P5_RangeRate";
            saveMovie = false;
            animateProblem(batchRunRangeRate.X_batchFilt(2:end,:)', stations, earthConst, titleText, movieTitle, saveMovie)

            fprintf("\n\nDone, have a great day!\n\n")
        case 'n'
            fprintf("\n\nDone, have a great day!\n\n")
        otherwise
            fprintf("\n\nDone, have a great day!\n\n")
    end
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
    LKFRun = runLKF(X0, zeros(size(X0)), P0, earthConst, scConst, stations, 6);

    choice = input(sprintf("Animate solution? [y/n]\n\tChoice:\t"), "s");
    switch choice
        case 'y'
            titleText = "ASEN 6080 Project 1 Problem 6 - No fixed stations, Batch";
            movieTitle = "../Videos/ASEN6080_Project1_P6_NoFix_Batch";
            saveMovie = false;
            animateProblem(batchRun.X_batchFilt(2:end,:)', stations, earthConst, titleText, movieTitle, saveMovie)

            titleText = "ASEN 6080 Project 1 Problems 6 - No fixed stations, LKF";
            movieTitle = "../Videos/ASEN6080_Project1_P6_NoFix_LKF";
            saveMovie = false;
            animateProblem(LKFRun.X_LKF, stations, earthConst, titleText, movieTitle, saveMovie)

            fprintf("\n\nDone, have a great day!\n\n")
        case 'n'
            fprintf("\n\nDone, have a great day!\n\n")
        otherwise
            fprintf("\n\nDone, have a great day!\n\n")
    end
elseif runQuestion == 4
    fprintf("\n\nRunning Question 6 - Station 337 fixed!\n\n")

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
    LKFRun = runLKF(X0, zeros(size(X0)), P0, earthConst, scConst, stations, 7);

    choice = input(sprintf("Animate solution? [y/n]\n\tChoice:\t"), "s");
    switch choice
        case 'y'
            titleText = "ASEN 6080 Project 1 Problem 6 - Station 337 fixed, Batch";
            movieTitle = "../Videos/ASEN6080_Project1_P6_Fixed337_Batch";
            saveMovie = false;
            animateProblem(batchRun.X_batchFilt(2:end,:)', stations, earthConst, titleText, movieTitle, saveMovie)

            titleText = "ASEN 6080 Project 1 Problems 6 - Station 337 fixed, LKF";
            movieTitle = "../Videos/ASEN6080_Project1_P6_Fixed337_LKF";
            saveMovie = false;
            animateProblem(LKFRun.X_LKF, stations, earthConst, titleText, movieTitle, saveMovie)

            fprintf("\n\nDone, have a great day!\n\n")
        case 'n'
            fprintf("\n\nDone, have a great day!\n\n")
        otherwise
            fprintf("\n\nDone, have a great day!\n\n")
    end
elseif runQuestion == 5
    fprintf("\n\nRunning Question 6 - Station 394 fixed!\n\n")

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
    LKFRun = runLKF(X0, zeros(size(X0)), P0, earthConst, scConst, stations, 8);

    choice = input(sprintf("Animate solution? [y/n]\n\tChoice:\t"), "s");
    switch choice
        case 'y'
            titleText = "ASEN 6080 Project 1 Problem 6 - Station 394 fixed, Batch";
            movieTitle = "../Videos/ASEN6080_Project1_P6_Fixed394_Batch";
            saveMovie = false;
            animateProblem(batchRun.X_batchFilt(2:end,:)', stations, earthConst, titleText, movieTitle, saveMovie)

            titleText = "ASEN 6080 Project 1 Problems 6 - Station 394 fixed, LKF";
            movieTitle = "../Videos/ASEN6080_Project1_P6_Fixed394_LKF";
            saveMovie = false;
            animateProblem(LKFRun.X_LKF, stations, earthConst, titleText, movieTitle, saveMovie)

            fprintf("\n\nDone, have a great day!\n\n")
        case 'n'
            fprintf("\n\nDone, have a great day!\n\n")
        otherwise
            fprintf("\n\nDone, have a great day!\n\n")
    end
elseif runQuestion == 6
    fprintf("\n\nRunning Extra code - LKF vs. Potter!\n\n")

    % Run LKF
    LKFRun = runLKF(X0, zeros(size(X0)), P0, earthConst, scConst, stations, 5);

    fprintf("\nPress 'Enter' to run the Potter Algorithm.\n")
    pause

    % Run Potter
    PotterRun = runPotter(X0, zeros(size(X0)), P0, earthConst, scConst, stations, 5);

    choice = input(sprintf("Animate solution? [y/n]\n\tChoice:\t"), "s");
    switch choice
        case 'y'
            titleText = "ASEN 6080 Project 1 Extra - LKF";
            movieTitle = "../Videos/ASEN6080_Project1_Extra_LKF";
            saveMovie = false;
            animateProblem(LKFRun.X_LKF, stations, earthConst, titleText, movieTitle, saveMovie)

            titleText = "ASEN 6080 Project 1 Extra - Potter";
            movieTitle = "../Videos/ASEN6080_Project1_Extra_Potter";
            saveMovie = false;
            animateProblem(PotterRun.X_Potter(:,2:end), stations, earthConst, titleText, movieTitle, saveMovie)

            fprintf("\n\nDone, have a great day!\n\n")
        case 'n'
            fprintf("\n\nDone, have a great day!\n\n")
        otherwise
            fprintf("\n\nDone, have a great day!\n\n")
    end
elseif runQuestion == 7
    fprintf("\n\nRunning Extra code - Reasonable P0 estimates!\n\n")

    % Make P0 on the order of meters and mm/s
    P0 = diag([1 1 1 1e-3 1e-3 1e-3 10 1e-4 1e-2 ...
               1 1 1 1 1 1 1 1 1]);

    % Run Batch
    measInclude = true(1,2);
    batchRun = runBatch(X0, zeros(size(X0)), P0, earthConst, scConst, stations, tspan, [], opt, 10, measInclude);
    
    fprintf("\nPress 'Enter' to run the LKF.\n")
    pause

    % Run LKF
    LKFRun = runLKF(X0, zeros(size(X0)), P0, earthConst, scConst, stations, 5);

    choice = input(sprintf("Animate solution? [y/n]\n\tChoice:\t"), "s");
    switch choice
        case 'y'
            titleText = "ASEN 6080 Project 1 Extra - Batch, reasonable P0";
            movieTitle = "../Videos/ASEN6080_Project1_Extra_Batch_P0";
            saveMovie = false;
            animateProblem(batchRun.X_batchFilt(2:end,:)', stations, earthConst, titleText, movieTitle, saveMovie)

            titleText = "ASEN 6080 Project 1 Extra - LKF, reasonable P0";
            movieTitle = "../Videos/ASEN6080_Project1_Extra_LKF_P0";
            saveMovie = false;
            animateProblem(LKFRun.X_LKF, stations, earthConst, titleText, movieTitle, saveMovie)

            fprintf("\n\nDone, have a great day!\n\n")
        case 'n'
            fprintf("\n\nDone, have a great day!\n\n")
        otherwise
            fprintf("\n\nDone, have a great day!\n\n")
    end
end
