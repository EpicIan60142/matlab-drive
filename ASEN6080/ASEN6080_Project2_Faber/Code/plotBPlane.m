function fig = plotBPlane(BdotR, BdotT, X_crossing, P_BPlane, STR2ECI, pConst, boundLevel, titleText, xLabel, yLabel, zLabel, figNum)
% Plots the location of the B vector in the B plane for a problem as well
% as the true crossing point and 3 sigma uncertainty
%   Inputs:
%       - BdotR: Component of B vector in the Rhat direction
%       - BdotT: Component of B vector in the That direction
%       - X_crossing: State at the B plane crossing
%       - P_BPlane: Covariance matrix at the B plane crossing
%       - STR2ECI: DCM to convert STR coordinates to ECI coordinates
%       - pConst: Planetary constants structure as defined in
%                 getPlanetConst.m
%       - boundLevel: What level of uncertainty to block, i.e. 1, 2, or 3
%                     sigma. To plot 3 sigma, pass in the number 3.
%       - titleText: Title to be displayed at the top of the plot as a 
%                    string
%       - xLabel: Label of the x-axis as a string
%       - yLabel: Label of the y-axis as a string
%       - zLabel: Label of the z-axis as a string
%       - figNum: Number of the figure, specify to plot multiple ellipses
%                 across multiple plotBPlane calls
%   Outputs:
%       - fig: Ellipse plot figure handle
%
%   By: Ian Faber, 04/15/2025
%

    % Plotting setup
I = imread('2k_earth_daymap.jpg'); % Image from https://www.solarsystemscope.com/textures/
theta0 = 0;
[earthX_s, earthY_s, earthZ_s] = sphere(20);
earthX_s = pConst.Ri*earthX_s;
earthY_s = pConst.Ri*earthY_s;
earthX = earthX_s*cos(theta0) - earthY_s*sin(theta0);
earthY = earthX_s*sin(theta0) + earthY_s*cos(theta0);
earthZ = pConst.Ri*earthZ_s;

    % Define helper variables
theta = linspace(0, 2*pi, 100);

    % Find eigenvectors and values of P_BPlane in R and T directions
[eigVec, Lambda] = eig(P_BPlane(2:3, 2:3));

    % Find axis scales
Lambda = chol(Lambda);
scales = [];
for k = 1:size(Lambda,1)
    scales = [scales; Lambda(k,k)];
end

    % Calculate basic circle
circle.t = cos(theta);
circle.r = sin(theta);

    % Scale circle
circleScaled.t = boundLevel*scales(1)*circle.t;
circleScaled.r = boundLevel*scales(2)*circle.r;

    % Rotate ellipse
rot = eigVec*[circleScaled.t; circleScaled.r];
ellipseRot.t = rot(1,:);
ellipseRot.r = rot(2,:);
ellipseRot.s = zeros(size(ellipseRot.t));

    % Convert ellipse to ECI coordinates and center around the Bplane
    % crossing
coords = [ellipseRot.s; ellipseRot.t; ellipseRot.r];
coords = STR2ECI*coords;
coords = coords + X_crossing(1:3)';
ellipse.x = coords(1,:);
ellipse.y = coords(2,:);
ellipse.z = coords(3,:);

    % Convert B vector to ECI coordinates
B = STR2ECI'*[0; BdotT; BdotR];

    % Make labels
BVecLabel = sprintf("B Plane Target: \nBdotR = %.4e km,\nBdotT = %.4e km", BdotR, BdotT);
crossingLabel = sprintf("Estimated B Plane Crossing");
uncertLabel = sprintf("+/- %.0f\\sigma B plane crossing uncertainty", boundLevel);

    % Plot the Bplane!
fig = figure(figNum);
hold on; grid on; axis equal
title(titleText, 'FontSize', 12)
    % Plot Earth
earth = surf(earthX, earthY, earthZ, 'FaceAlpha', 0.5);
set(earth,'FaceColor','texturemap','cdata',I,'edgecolor','none');
    % Plot B plane target
BVec = quiver3(0,0,0,B(1),B(2),B(3),1,'Filled','b-','LineWidth',2);
    % Plot B plane crossing and uncertainty
crossing = plot3(X_crossing(1), X_crossing(2), X_crossing(3),'kx', 'MarkerSize', 10);
uncert = plot3(ellipse.x, ellipse.y, ellipse.z, 'r-');
    % Labels and legend
xlabel(xLabel); ylabel(yLabel); zlabel(zLabel); view([30 35]);
legend([BVec, crossing, uncert], [BVecLabel, crossingLabel, uncertLabel], 'location', 'northwest');

end