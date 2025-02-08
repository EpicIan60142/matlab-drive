%% Ian Faber - SID: 108577813
% ASEN 2012: Coding Challenge 2
% Last modified: 9/9/21

%% Housekeeping
 % Clear command window
 % Close out of any open figures
clc; close all;

%% Given Information
syms g0 i_sp_bar sigma_i_sp m_s_bar sigma_m_s m_p_bar sigma_m_p positive
syms delta_v_bar sigma_delta_v real



%g0 = 9.81; % [CONSTANT] Acceleration due to gravity (m/s^2)

%i_sp_bar = 459; % Best estimate for specific impulse (s)
sigma_i_sp = 11; % Standard deviation for specific impulse (s)

%m_s_bar = 13050; % Best estimate for structure mass (kg)
sigma_m_s = 60; % Standard deviation for structure mass (kg)

%m_p_bar = 71800; % Best estimate for propellant mass (kg)
sigma_m_p = 300; % Standard deviation for propellant mass (kg)


% Note: Any quantity can be expressed as its best estimate +/- its absolute
% uncertainty, and the absolute uncertainty can be expressed as its best estimate * its fractional uncertainty.
    
% Note: You may create your own functions (e.g. for propagating error for addition, multiplication, etc.)
% or you may propagate all error from the main script.

%% Applying the ideal rocket equation
delta_v_bar = g0*i_sp_bar*log((m_s_bar+m_p_bar)/m_s_bar) % Best estimate of Delta_v (change in speed [m/s])

divIsp = diff(delta_v_bar, i_sp_bar)
divMs = diff(delta_v_bar, m_s_bar)
divMp = diff(delta_v_bar, m_p_bar)

sigma_delta_v = sqrt((divIsp*sigma_i_sp)^2+(divMs*sigma_m_s)^2+(divMp*sigma_m_p)^2) % Absolute uncertainty of Delta_v (change in speed [m/s])

del_v = vpa(subs(delta_v_bar,[g0,i_sp_bar,m_s_bar,m_p_bar],[9.81,459,13050,71800]),4)
sig_v = vpa(subs(sigma_delta_v,[g0,i_sp_bar,m_s_bar,m_p_bar],[9.81,459,13050,71800]),4)


