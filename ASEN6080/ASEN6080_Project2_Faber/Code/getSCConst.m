function sc = getSCConst(part)
% Function that returns relevant spacecraft constants for OD problems. 
% Update as necessary!
%   Inputs:
%       - part: Which part of project 3 the constants are for
%   Outputs:
%       - sc: Structure of spacecraft paremeters at the start of an OD 
%             problem. Made of the following fields:
%           - AtoM: Spacecraft area to mass ratio [m^2/kg]
%           - C_R: Expected coefficient of reflectivity [n.d.]
%           - X0_cart: Initial cartesian state in km and km/s:
%                      [X0; Y0; Z0; Xdot0; Ydot0; Zdot0]
%           - Phi: Solar pressure flux at the location of the spacecraft's
%                  orbit for the problem (i.e. at 1 AU) [W/km^2]
%
%   By: Ian Faber, 04/03/2025
%

    % Physical Parameters
sc.AtoM = 0.01/(1000^2); % Area to mass ratio, m^2/kg -> km^2/kg
if part == 2
    sc.C_R = 1.2;
    sc.X0_cart = [
                    -274096790.0;
                    -92859240.0;
                    -40199490.0;
                    32.67;
                    -8.94;
                    -3.88;
                 ]; % km, km/s
else
    sc.C_R = 1.0;
    sc.X0_cart = [
                    -274096770.76544;
                    -92859266.4499061;
                    -40199493.6677441;
                    32.6704564599943;
                    -8.93838913761049;
                    -3.87881914050316;
                 ]; % km, km/s
end
sc.Phi = 1357; % W/m^2 -> kg/s^3, no need for unit conversion...

end