function sc = getSCConst()
% Function that returns relevant spacecraft constants for OD problems. 
% Update as necessary!
%   Inputs:
%       - None
%   Outputs:
%       - sc: Structure of spacecraft paremeters at the start of an OD 
%             problem. Made of the following fields:
%           - A: Cross-sectional area [km^2]
%           - m: Mass [kg]
%           - Cd: Drag coefficient [n.d.]
%           - X0_cart: Initial cartesian state in km and km/s:
%                      [X0; Y0; Z0; Xdot0; Ydot0; Zdot0]
%
%   By: Ian Faber, 02/10/2025
%

    % Physical Parameters
sc.A = 3/(1000^2); % m^2 -> km^2
sc.m = 970;
sc.Cd = 2;
sc.X0_cart = [
                757700.0;
                5222607.0;
                4851500.0;
                2213.21;
                4678.34;
                -5371.30;
             ]/1000; % m -> km

end