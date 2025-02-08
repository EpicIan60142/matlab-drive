%% Ian Faber - SID: 108577813
% ASEN 2012: Coding Challenge 2
% Last modified: 9/9/21

%% Housekeeping
 % Clear command window
 % Close out of any open figures
clc; close all;

%% Given Information
g0 = 9.81; % [CONSTANT] Acceleration due to gravity (m/s^2)

i_sp_bar = 459; % Best estimate for specific impulse (s)
sigma_i_sp = 11; % Standard deviation for specific impulse (s)

m_s_bar = 13050; % Best estimate for structure mass (kg)
sigma_m_s = 60; % Standard deviation for structure mass (kg)

m_p_bar = 71800; % Best estimate for propellant mass (kg)
sigma_m_p = 300; % Standard deviation for propellant mass (kg)


% Note: Any quantity can be expressed as its best estimate +/- its absolute
% uncertainty, and the absolute uncertainty can be expressed as its best estimate * its fractional uncertainty.
    
% Note: You may create your own functions (e.g. for propagating error for addition, multiplication, etc.)
% or you may propagate all error from the main script.

[f, sigf] = addSubError([m_s_bar, m_p_bar], [sigma_m_s, sigma_m_p], true)
[g, sigg] = multDivError([f, m_s_bar], [sigf, sigma_m_s], false)
[h, sigh] = natLogError(g,sigg)
[v,sigv] = multDivError([h, i_sp_bar], [sigh, sigma_i_sp], true)

%% Applying the ideal rocket equation
delta_v_bar = g0*v % Best estimate of Delta_v (change in speed [m/s])
sigma_delta_v = g0*sigv % Absolute uncertainty of Delta_v (change in speed [m/s])

function [z, sigz] = addSubError(X, sigX, adding)
    if(adding)
        z = sum(X);
    else
        z = X(1)-sum(X(2:end));
    end

    sigz = sqrt(sum(sigX.^2));
end

function [z, sigz] = multDivError(X, sigX, multiplying)
    z = X(1);
    if(multiplying)
       for k = 2:length(X)
          z = z*X(k);
       end
    else
       for k = 2:length(X)
          z = z/X(k);
       end 
    end
    
    sigz = z*sqrt(sum((sigX./X).^2));
end

function [z, sigz] = natLogError(x, sigx)
    z = log(x);
    sigz = sigx/x;
end