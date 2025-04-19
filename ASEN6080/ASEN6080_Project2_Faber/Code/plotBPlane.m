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

    % Convert B vector to ECI coordinates
BVec = STR2ECI'*[0; BdotT; BdotR];

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
    % target
coords = [ellipseRot.s; ellipseRot.t; ellipseRot.r];
coords = STR2ECI*coords;
coords = coords + BVec;%X_crossing(1:3)';
ellipse.x = coords(1,:);
ellipse.y = coords(2,:);
ellipse.z = coords(3,:);

    % Make Bplane object
N = 10;
s = zeros(1,N);
t = linspace(0,2*BdotT,N);
r = linspace(0,2*BdotR,N);

% points = [s;t;r];
% points = STR2ECI*points;
% A = STR2ECI(1,1); B = STR2ECI(2,1); C = STR2ECI(3,1); D = 0;
% [X_Bplane, Y_Bplane] = meshgrid(points(1,:),points(2,:));
% Z_Bplane = -1./C*(A*X_Bplane + B*Y_Bplane + D);

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
Bvec = quiver3(0,0,0,BVec(1),BVec(2),BVec(3),1,'Filled','b-','LineWidth',2);
    % Plot STR axes
Svec = quiver3(0,0,0,STR2ECI(1,1),STR2ECI(2,1),STR2ECI(3,1), pConst.Ri, 'filled', 'b--', 'LineWidth', 2);
Tvec = quiver3(0,0,0,STR2ECI(1,2),STR2ECI(2,2),STR2ECI(3,2), pConst.Ri, 'filled', 'r--', 'LineWidth', 2);
Rvec = quiver3(0,0,0,STR2ECI(1,3),STR2ECI(2,3),STR2ECI(3,3), pConst.Ri, 'filled', 'k--', 'LineWidth', 2);
    % Plot B plane crossing and uncertainty
crossing = plot3(X_crossing(1), X_crossing(2), X_crossing(3),'kx', 'MarkerSize', 10);
uncert = plot3(ellipse.x, ellipse.y, ellipse.z, 'r-');
    % Plot B plane
Bplane = surf(X_Bplane, Y_Bplane, Z_Bplane, 'EdgeColor', 'k', 'FaceColor', 'g', 'FaceAlpha', 0.4);
    % Labels and legend
xlabel(xLabel); ylabel(yLabel); zlabel(zLabel); view([30 35]);
legend([Bvec, crossing, uncert, Svec, Tvec, Rvec, Bplane], [BVecLabel, crossingLabel, uncertLabel, "S Vector", "T Vector", "R Vector", "Bplane"], 'location', 'northwest');

end