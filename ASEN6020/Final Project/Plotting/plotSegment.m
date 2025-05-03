function fig = plotSegment(cubesat, ring, t, X, figNum, titleText, xLabel, yLabel, zLabel)
% Function that plots a single segment of the cubesat race course
%   Inputs:
%       - cubesat: Cubesat structure for the cubesat of interest
%       - ring: The last ring of this segment
%       - t: Time vector for this segment
%       - X: State vector for this segment
%       - figNum: Figure number for this segment
%       - 
%   Outputs:
%       - fig: Figure handle for this segment

fig = figure(figNum);
title(titleText);
hold on; grid on; axis equal
% if figNum - 1 == 1
%     plotRing(startRing, 'g-');
% else
%     plotRing(ring.params.lastRing, 'g-');
% end
plotRing(ring.params.lastRing, 'g-');
plot3(X(1,1), X(1,2), X(1,3), 'k', 'Marker', cubesat.marker);
plotRing(ring, 'r-');
plot3(X(:,1), X(:,2), X(:,3), '--', 'Color', cubesat.color);
xlabel(xLabel); ylabel(yLabel); zlabel(zLabel);
view([30 35]); drawnow;

end