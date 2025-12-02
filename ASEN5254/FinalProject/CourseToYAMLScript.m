%% ASEN 5254 Race Course to YAML script
% By: Ian Faber, 12/02/2025

%% Housekeeping
clc; clear; close all;

addpath("../../ASEN6020/Final Project/OfficialRuns/");

%% Setup
% Choose dataset
load("Run_05-May-2025_215651_Scenario2.mat");

% Make YAML file object
file = fopen("RaceCourse.yml",'w');

%% Populate YAML file
    % Course header
fprintf(file, "Map:\n");

    % World dimensions
fprintf(file, "\tDimensions:\n");

fig = plotCourse(startRing,rings,endRing,cubesats,1,"Race Course","Radial [m]","Along-Track [m]","Cross-Track [m]",true);

lims = zeros(3,2);
lims(1,:) = 1.1*fig.Children(3).XLim;
lims(2,:) = 1.1*fig.Children(3).YLim;
lims(3,:) = 1.1*fig.Children(3).ZLim;

for k = 1:3
    for kk = 1:2
        fprintf(file,"\t\t- %.3f\n",lims(k,kk));
    end
end

    % Rings
rings = [startRing; rings; endRing];
fprintf(file, "\tRings:\n");
for k = 1:length(rings)
    fprintf(file, "\t\t- ring%.0f\n",k-1);
        % Ring center
    for kk = 1:3
        fprintf(file, "\t\t\t- %.3f\n", rings(k).center(kk));
    end
        % Ring normal vector
    for kk = 1:3
        fprintf(file, "\t\t\t- %.3f\n", rings(k).normal(kk));
    end
        % Ring semi-major and semi-minor axes
    for kk = 1:2
        fprintf(file, "\t\t\t- %.3f\n", rings(k).S(kk,kk));
    end
        % DCM to Inertial from Ring frame
    for kk = 1:3
        for idx = 1:3
            fprintf(file, "\t\t\t- %.3f\n", rings(k).NR(kk, idx));
        end
    end
end

    % Cubesats header
fprintf(file, "Agents:\n");

    % Cubesats
for k = 1:length(cubesats)
    fprintf(file, "\tagent%.0f\n", k-1);
        % Cubesat dynamics
    fprintf(file, "\t\tModel: CubesatCWH\n");
        % Start state
    fprintf(file, "\t\tStart:\n");
    for kk = 1:6
        fprintf(file, "\t\t\t- %.3f\n", cubesats(k).X0(kk));
    end
        % Max control
    fprintf(file, "\t\tuMax: %.3f\n", cubesats(k).uMax);
end

fclose(file);

