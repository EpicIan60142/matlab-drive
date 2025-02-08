clear; clc; close all;

w0 = 2001;
L = 27.25;
maxP = 50;

score = zeros(1,maxP);
cost = zeros(1,maxP);
error = zeros(1,maxP);

for p = 1:maxP
    discretized = discretize_load(p, L, w0);
    [score(p), cost(p), error(p), ~, ~] = bestTree(p, discretized, w0, L);    
end

chosenP = 15;

[~, ~, ~, Mdist, Mpoint] = bestTree(chosenP, discretize_load(chosenP, L, w0), w0, L);

figure
subplot(3,1,1)
hold on;
title("Wiffle Tree cost vs. p point loads");
plot(cost)
xline(chosenP);
xlabel("Number of point loads");
ylabel("Cost ($)");
legend("Cost","Chosen number of point loads: " + num2str(chosenP),'Location','best')
hold off;

subplot(3,1,2)
hold on;
title(" Wiffle Tree error vs. p point loads");
plot(error)
xline(chosenP);
xlabel("Number of point loads");
ylabel("Error (N*m)");
legend("Error","Chosen number of point loads: " + num2str(chosenP),'Location','best')
hold off;

subplot(3,1,3)
hold on;
title("Wiffle Tree score vs. p point loads");
plot(score)
xline(chosenP);
xlabel("Number of point loads");
ylabel("Wiffle Tree score");
legend("Score","Chosen number of point loads: " + num2str(chosenP),'Location','best')
hold off;

figure
subplot(2,1,1)
hold on;
title("Analytical beam moment and point beam moment, unzoomed");
plot(Mdist);
plot(Mpoint);
xlabel("Percent of total length (%)");
ylabel("Moment (N*m)");
legend("Analytical moment", "Point moment", 'Location', 'best');
hold off;

zoom = 10;

subplot(2,1,2)
hold on;
title("Analytical beam moment and point beam moment, zoomed in");
plot(Mdist(length(Mdist)/2-zoom:length(Mdist)/2+zoom));
plot(Mpoint(length(Mpoint)/2-zoom:length(Mpoint)/2+zoom));
xlabel("Percent of total length (%)");
ylabel("Moment (N*m)");
legend("Analytical moment", "Point moment", 'Location', 'best');
hold off;
