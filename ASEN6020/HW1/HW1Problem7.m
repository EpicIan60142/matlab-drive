%% ASEN 6020 HW 1 Problem 7 Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
    % Physical params
mu = 398600.4415; % km^3/s^2
r = 7000:3000:30000; % km
% r = 10000;

    % Solution space - ratio of l/r
ratio = [linspace(0.1,1,100), linspace(2,10,99)]; % l/r - make 10 ratios below 1 and 10 ratios at or above 1

    % Helper functions
V_lc = @(rad) sqrt(mu./rad); % km/s
l = @(rat, r) rat*r; % km

    % Cost functions
J1i = @(Vinf, r) (sqrt(Vinf.^2 + 2*V_lc(r).^2))/V_lc(r) - 1;
J2ia = @(Vinf, el, r) sqrt((2*el)./(el+1)) - 1 + sqrt(Vinf.^2 + 2*V_lc(l(el,r)).^2)/V_lc(r) - sqrt(2./(el.*(el+1))); % el is script l
J2ip = @(Vinf, el, r) 1 - sqrt((2.*el)./(el+1)) + sqrt(Vinf.^2 + 2*V_lc(l(el,r)).^2)/V_lc(r) - sqrt(2./(el.*(el+1)));

%% Solve for optimal costs

X = []; Y = [];
cost_1imp = []; cost_2impA = []; cost_2impP = [];
for k = 1:length(r)
        % Solution space - Vinf
    Vinf = 0:0.1:2*V_lc(r(k));

        % Create grids
    [VINF, RATIO] = meshgrid(Vinf, ratio);
    [VINF_2iA, RATIO_2iA] = meshgrid(Vinf, ratio(ratio >= 1)); % Apoapsis maneuver only valid when l > r -> ratio > 1
    [VINF_2iP, RATIO_2iP] = meshgrid(Vinf, ratio(ratio <= 1)); % Periapsis maneuver only valid when l < r -> ratio < 1
    X = [X; {VINF, VINF_2iA, VINF_2iP}];
    Y = [Y; {RATIO, RATIO_2iA, RATIO_2iP}];

        % Calculate costs
    cost_1imp = [cost_1imp; {J1i(VINF, r(k))}];
    cost_2impA = [cost_2impA; {J2ia(VINF_2iA, RATIO_2iA, r(k))}];
    cost_2impP = [cost_2impP; {J2ip(VINF_2iP, RATIO_2iP, r(k))}];
end
%% Plot costs at each r
for k = 1:length(r)
    figure
    hold on; grid on;
    titleText = sprintf("Escape Maneuver Cost Comparison, r = %.1f km", r(k));
    title(titleText);
    surf(X{k,1}, Y{k,1}, cost_1imp{k}, 'EdgeColor', 'none', 'FaceColor', 'y', 'FaceAlpha', 0.9);
    surf(X{k,2}, Y{k,2}, cost_2impA{k}, 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.4);
    surf(X{k,3}, Y{k,3}, cost_2impP{k}, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.75);
    patch(sqrt(2)*V_lc(r(k))*ones(4,1), [0;0;ratio(end);ratio(end)], [0; 3; 3; 0], 'k', 'FaceAlpha', 0.25)
    xlabel("V_{\infty} [km/s]"); ylabel("l/r"); zlabel("\DeltaV Cost [km/s]")
    legend("1-Impulse", "2-Impulse, Apoapsis", "2-Impulse, Periapsis", "V_{\infty} = \\sqrt(2)V_{lc}(r)",'Location', 'northwest')
    view([-45 35])
end


