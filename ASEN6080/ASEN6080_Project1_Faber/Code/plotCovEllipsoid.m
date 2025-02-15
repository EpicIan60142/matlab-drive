function fig = plotCovEllipsoid(C, mu, titleText, xLabel, yLabel, zLabel, cBarText)
% Function that plots the covariance ellipsoid for a given 3x3 covariance
% matrix
%   Inputs:
%       - C: 3x3 covariance matrix
%       - mu: Mean of the data for the covariance matrix
%       - titleText: Title to be displayed at the top of the plot as a 
%                    string
%       - xLabel: Label of the x-axis as a string
%       - yLabel: Label of the y-axis as a string
%       - zLabel: Label of the z-axis as a string
%       - cBarText: Label of the colorbar of the ellispoid
%   Outputs:
%       - fig: Ellispoid plot figure handle
%
%   By: Ian Faber, 02/24/2025
%

    % Define helper variables
theta = linspace(0,2*pi,100); % Latitude/Azimuth
phi = linspace(0,pi,100); % Longitude/Elevation

    % Find eigenvectors and eigenvalues of C
[eigVec, Lambda] = eig(C);

    % Take square root of eigenvalues
Lambda = chol(Lambda);
axes = [];
for k = 1:size(Lambda,1)
    axes = [axes; Lambda(k,k)];
end

    % Calculate ellipse coordinates
[X, Y, Z] = sphere(100);
Xscaled = axes(1)*X;
Yscaled = axes(2)*Y;
Zscaled = axes(3)*Z;

coords = [];
coords(:,:,1) = Xscaled;
coords(:,:,2) = Yscaled;
coords(:,:,3) = Zscaled;

dist = zeros(size(Xscaled));
coordsRot = zeros(size(coords));
for k = 1:size(coords,1)
    for kk = 1:size(coords,2)
        coordsRot(k,kk,:) = eigVec*squeeze(coords(k,kk,:)) + mu;
        dist(k,kk) = norm(squeeze(coordsRot(k,kk,:))-mu);
    end
end
% coordsRot = tensorprod(eigVec',coords,1,3) + mu;
% coords = eigVec*[axes(1)*sin(phi).*cos(theta); axes(2)*sin(phi).*sin(theta); axes(3)*cos(phi)];
coords = coords + reshape(mu,1,1,3);

    % Plot ellipse
fontSize = 10;
scale = 3.0;
fig = figure;
hold on; grid on; axis equal
title(titleText, 'FontSize', fontSize)
% oneSig = surf(coordsRot(:,:,1), coordsRot(:,:,2), coordsRot(:,:,3), 'FaceColor', 'b', 'FaceAlpha', 1, 'LineStyle', 'none');
% twoSig = surf(coordsRot(:,:,1), coordsRot(:,:,2), coordsRot(:,:,3), 'FaceColor', 'r', 'FaceAlpha', 0.75, 'LineStyle', 'none');
% threeSig = surf(coordsRot(:,:,1), coordsRot(:,:,2), coordsRot(:,:,3), 'FaceColor', 'y', 'FaceAlpha', 0.5, 'LineStyle', 'none');
surf(coordsRot(:,:,1), coordsRot(:,:,2), coordsRot(:,:,3), dist, 'FaceAlpha',0.5, 'LineStyle','none')
a = quiver3(mu(1),mu(2),mu(3),eigVec(1,1),eigVec(2,1),eigVec(3,1),scale*Lambda(1,1),'filled','b');
b = quiver3(mu(1),mu(2),mu(3),eigVec(1,2),eigVec(2,2),eigVec(3,2),scale*Lambda(2,2),'filled','r');
c = quiver3(mu(1),mu(2),mu(3),eigVec(1,3),eigVec(2,3),eigVec(3,3),scale*Lambda(3,3),'filled','k');
xlabel(xLabel); ylabel(yLabel); zlabel(zLabel);
legend([oneSig, twoSig, threeSig, a,b,c], ["1\sigma", "2\sigma", "3\sigma","v_1", "v_2", "v_3"], 'location', 'eastoutside')
view([30 35]); colormap("turbo"); c = colorbar('southoutside');
c.Label.String = cBarText; drawnow;
% drawnow;
end