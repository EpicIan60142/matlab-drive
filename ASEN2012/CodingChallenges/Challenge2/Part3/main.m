%% Ian Faber - SID: 108577813
% ASEN 2012: Coding Challenge 2
% Last modified: 9/14/21

%% Housekeeping
clc % Clear command window
close all % Close out of any open figures

%% Given Information
g0 = 9.81; % [CONSTANT] Acceleration due to gravity (m/s^2)

i_sp_bar = 459; % Best estimate for specific impulse (s)
sigma_i_sp = 11; % Standard deviation for specific impulse (s)

m_s_bar = 13050; % Best estimate for structure mass (kg)
sigma_m_s = 60; % Standard deviation for structure mass (kg)

m_p_bar = 71800; % Best estimate for propellant mass (kg)
sigma_m_p = 300; % Standard deviation for propellant mass (kg)

%% Monte Carlo Simulation

% Number of Monte Carlo samples. This should be sufficiently high such that your simulation
% starts to approach the analytical solution.
n_monte_carlo = 10000; 

% Here, you want to create arrays of normally distributed values for each of these parameters,
% with the correct mean and error values. This can be done either with vectors or with for-
% loops, though the latter is much more lines and runs significantly slower with large n

Isp = i_sp_bar + sigma_i_sp*randn(n_monte_carlo,1); % Randomly distributed Isp values (s)
m_s = m_s_bar + sigma_m_s*randn(n_monte_carlo,1); % Randomly distributed structure mass values (kg)
m_p = m_p_bar + sigma_m_p*randn(n_monte_carlo,1); % Randomly distributed propellant mass values (kg)

% Use the arrays above to determine the randomly distributed delta-v values, along with
% the associated mean and error
delta_v = g0*Isp.*log((m_s+m_p)./m_s); % Determine values for delta v using the random variables
delta_v_bar = mean(delta_v); % Determine the mean of that data
sigma_delta_v = std(delta_v); % As well as the standard deviation


%% Output and Format Histogram

% Plot the histogram of the data. With enough points, this should approach a bell-curve
figure
hold on

CUgold = '#CFB87C';
CUsilver = '#A2A4A3';
delta_v_bar_graph = round(delta_v_bar,0);
sigma_delta_v_graph = round(sigma_delta_v,0);

[title,subtitle] = title('Monte Carlo Simulation of the Tsiolkovsky Rocket Equation',['N = ', num2str(n_monte_carlo),' samples']);
xlabel('\Deltav')
ylabel('Frequency')
Histo = histogram(delta_v,'FaceColor',CUgold,'edgeAlpha',0.5);
v_line = xline(delta_v_bar,'LineWidth',2,'color','r');
set(gca,'color',CUsilver)
legend([v_line,Histo],{['\Deltav = ',num2str(delta_v_bar_graph),'\pm',num2str(sigma_delta_v_graph),' m/s'],'Calculated \Deltav'},'Location','best')