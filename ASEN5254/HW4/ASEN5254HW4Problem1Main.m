%% ASEN 5254 HW 4 Problem 1 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
theta = linspace(0, pi/4, 12);

obstacleX = [0, 1, 0, 0];
obstacleY = [0, 2, 2, 0];
obstacle = [obstacleX; obstacleY];
robotX = obstacleX;
robotY = obstacleY;
invRobotX = [-1 0 0, -1];
invRobotY = [-2 -2 0, -2];

%% Make c-space obstacle for every theta
cSpaceObs = [];
for k = 1:length(theta)
    invRobot = [cos(theta(k)) -sin(theta(k)); sin(theta(k)) cos(theta(k))]*[invRobotX; invRobotY];

    [invRobotX, invRobotY] = poly2ccw(invRobot(1,:), invRobot(2,:));

    invRobot = [invRobotX; invRobotY];

    ii = 1; jj = 1; n = length(obstacleX); m = length(invRobotX);
    cSpaceObsX = []; cSpaceObsY = [];
    while ii ~= n + 1 && jj ~= m + 1
        %fprintf("theta = %.3f deg, ii = %.0f, jj = %.0f\n", rad2deg(theta(k)), ii, jj)

        cSpaceObsX = [cSpaceObsX, obstacleX(ii) + invRobotX(jj)];
        cSpaceObsY = [cSpaceObsY, obstacleY(ii) + invRobotY(jj)];
        
        if ii == length(obstacleX)
            obsVec = obstacle(:,1) - obstacle(:,ii);
        else
            obsVec = obstacle(:,ii+1) - obstacle(:,ii);
        end
        obsAngle = wrapTo2Pi(atan2(obsVec(2), obsVec(1)));

        if jj == length(invRobotX)
            invVec = invRobot(:,1) - invRobot(:,jj);
        else
            invVec = invRobot(:,jj+1) - invRobot(:,jj);
        end
        invAngle = wrapTo2Pi(atan2(invVec(2), invVec(1)));

        if(obsAngle < invAngle)
            ii = ii + 1;
        elseif(obsAngle > invAngle)
            jj = jj + 1;
        else
            ii = ii + 1;
            jj = jj + 1;
        end
    end
    cSpaceObs = [cSpaceObs; {[cSpaceObsX; cSpaceObsY]}];
end

figure;
hold on; grid on;
title("C-space object")
for k = 1:length(theta)
    plot3(cSpaceObs{k}(1,:), cSpaceObs{k}(2,:), rad2deg(theta(k))*ones(size(cSpaceObs{k}(1,:))), 'b-')
    xlabel("X [m]"); ylabel("Y [m]"); zlabel("\theta [deg]")

    view([30 35])






end