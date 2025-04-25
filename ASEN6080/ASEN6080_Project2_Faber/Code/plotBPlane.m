function fig = plotBPlane(BdotR, BdotT, X_crossing, P_BPlane, STR2ECI, pConst, boundLevel, titleText, xLabel, yLabel, zLabel, ellipseLabel, ellipseColor, BVecLabel, figNum, newFig)
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
BVec = STR2ECI*[0; BdotT; BdotR];

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
% s = zeros(1,N);
mag = 1.25*max(BdotT, BdotR);
t = linspace(0,mag,N);
r = linspace(0,mag,N);

A = 0; B = 0; C = 1; D = 0; % Pull out normal vector
[T_Bplane, R_Bplane] = meshgrid(t, r);
S_Bplane = -1./C*(A*T_Bplane + B*R_Bplane + D);

coords = zeros(3,10,10);
coords(1,:,:) = S_Bplane;
coords(2,:,:) = T_Bplane;
coords(3,:,:) = R_Bplane;

coords = tensorprod(STR2ECI',coords,1,1);

    % Make labels
% BVecLabel = sprintf("B Plane Target after %.3f days of data: \nBdotR = %.4e km,\nBdotT = %.4e km", BdotR, BdotT);
crossingLabel = sprintf("Trajectory propagated to LTOF");
% uncertLabel = sprintf("+/- %.0f\\sigma B plane crossing uncertainty", boundLevel);

    % Plot the Bplane!
fig = figure(figNum); fig.WindowState = 'maximized';
if newFig
    tl = tiledlayout(1,3);
    title(tl, titleText, 'FontSize', 12)
    nexttile(1)
        hold on; grid on; axis equal
        title("B plane in ECI frame")
            % Plot Earth
        earth = surf(earthX, earthY, earthZ, 'FaceAlpha', 0.5);
        set(earth,'FaceColor','texturemap','cdata',I,'edgecolor','none');
        
            % Plot STR axes
        Svec = quiver3(0,0,0,STR2ECI(1,1),STR2ECI(2,1),STR2ECI(3,1), pConst.Ri, 'filled', 'b--', 'LineWidth', 2);
        Tvec = quiver3(0,0,0,STR2ECI(1,2),STR2ECI(2,2),STR2ECI(3,2), pConst.Ri, 'filled', 'r--', 'LineWidth', 2);
        Rvec = quiver3(0,0,0,STR2ECI(1,3),STR2ECI(2,3),STR2ECI(3,3), pConst.Ri, 'filled', 'k--', 'LineWidth', 2);
            % Plot B plane crossing and uncertainty
        crossing = plot3(X_crossing(1), X_crossing(2), X_crossing(3),'kx', 'MarkerSize', 10);
        % plot3(ellipse.x, ellipse.y, ellipse.z, '-', 'Color', ellipseColor, 'DisplayName', ellipseLabel);
            % Plot B plane target
        % quiver3(0,0,0,BVec(1),BVec(2),BVec(3),1,'Filled','-','Color', ellipseColor, 'LineWidth',2,'DisplayName',BVecLabel);
            % Plot B plane
        % Bplane = surf(X_Bplane, Y_Bplane, Z_Bplane, 'EdgeColor', 'none', 'FaceColor', 'g', 'FaceAlpha', 0.4);
        Bplane = surf(squeeze(coords(1,:,:)),squeeze(coords(2,:,:)),squeeze(coords(3,:,:)), 'EdgeColor', 'none', 'FaceColor', 'g', 'FaceAlpha', 0.4);
            % Labels and viewing angle
        xlabel(xLabel); ylabel(yLabel); zlabel(zLabel); view([-15 25]);
    nexttile(2)
        hold on; grid on;
        title("Bplane in STR frame")
            % Bplane target
        plot(BdotT, BdotR, 'x', 'Color', ellipseColor);
            % Uncertainty
        plot(ellipseRot.t + BdotT, ellipseRot.r + BdotR, '-', 'Color', ellipseColor);
        xlabel("T [km]"); ylabel("R [km]");
    
    lgnd = legend([crossing, Svec, Tvec, Rvec, Bplane], [crossingLabel, "S Vector", "T Vector", "R Vector", "B Plane"], 'location', 'layout');
    lgnd.Layout.Tile = 3;

    nexttile(1)
        hold on;
            % Plot B plane target in ECI
        quiver3(0,0,0,BVec(1),BVec(2),BVec(3),1,'Filled','-','Color', ellipseColor, 'LineWidth',2, 'DisplayName', BVecLabel);
            % Plot uncertainty in ECI
        plot3(ellipse.x, ellipse.y, ellipse.z, '-', 'Color', ellipseColor, 'DisplayName', ellipseLabel);    
else
        % Allow for dynamic B plane target and uncertainty plotting
    nexttile(1)
        hold on;
            % Plot B plane target in ECI
        quiver3(0,0,0,BVec(1),BVec(2),BVec(3),1,'Filled','-','Color', ellipseColor, 'LineWidth',2, 'DisplayName', BVecLabel);
            % Plot uncertainty in ECI
        plot3(ellipse.x, ellipse.y, ellipse.z, '-', 'Color', ellipseColor, 'DisplayName', ellipseLabel);   
    nexttile(2)
        hold on;
            % Plot B plane target in STR
        plot(BdotT, BdotR, 'x', 'Color', ellipseColor);
            % Plot uncertainty in STR
        plot(ellipseRot.t + BdotT, ellipseRot.r + BdotR, '-', 'Color', ellipseColor);
end

drawnow;

end