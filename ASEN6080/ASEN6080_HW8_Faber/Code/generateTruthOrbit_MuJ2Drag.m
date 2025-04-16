function data = generateTruthOrbit_MuJ2Drag(pConst, oConst, scConst, x0Perturb, filename, tspan)
% Function that generates truth orbit and measurement data based on a set
% of planetary constants, orbital parameters, and potential initial orbital
% perturbations including only dynamics from Mu and J2.
%   Inputs:
%       - pConst: Planetary constants structure as specified by
%                 getPlanetConst.m
%       - oConst: Orbital elements structure as specified by
%                 getOrbitConst.m
%       - scConst: Spacecraft cosntants structure as specified by
%                  getSCConst.m
%       - x0Perturb: Initial orbit perturbation organized as follows:
%                    [dX; dY; dZ; dXdot; dYdot; dZdot]
%       - filename: Name of the file to save truth data to. 
%                   Ex: "HW2Problem1Data.mat"
%       - tspan: Time period to generate data over as a vector
%   Outputs:
%       - data: Data structure with the following fields:
%               - X0: Initial cartesian state based on the orbital elements
%                     in oConst, organized as [X0; Y0; Z0; Xdot0; Ydot0;
%                     Zdot0]
%               - x0Perturb: Initial state perturbation as specified by
%                            x0Perturb
%               - X_ref: Truth orbit trajectory, organized as a matrix of
%                        (numOrbits*T)/dt 1x6 state row vectors:
%                        [ 
%                          [X, Y, Z, Xdot, Ydot, Zdot]; 
%                          [X, Y, Z, Xdot, Ydot, Zdot]; ...
%                        ]
%               - t_ref: Time vector corresponding to each row of X_ref
%               - tspan: The time span used by ode45 to generate
%                        measurements and the orbit trajectory
%               - opt: Options used by ode45 to generate the measurements
%                      and orbit trajectory
%
%   Note: It is recommended to just load the created mat file after running
%   this function, which will provide all of the fields defined in the data
%   output as their own variables!
%
%   By: Ian Faber, 02/01/2025
%

orbital = oConst;

    % Convert orbital elements to cartesian for the initial orbit state
X0 = convOrbitalToCartesian([pConst.mu; orbital.a; orbital.e; orbital.i; orbital.RAAN; orbital.argPeri; orbital.truAnom]);

opt = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% XPhi0 = [X0+x0Perturb; reshape(eye(length(X0)),length(X0)^2,1)];

[t_ref, X_ref] = ode45(@(t,X)orbitEOM_MuJ2Drag(t,X,pConst, scConst), tspan, X0+x0Perturb, opt);
% [t_ref, XPhi_ref] = ode45(@(t,XPhi)STMEOM_J2(t,XPhi,planetConst.mu,planetConst.J2,planetConst.Ri),tspan,XPhi0,opt);

% X_ref = XPhi_ref(:,1:length(X0));

    % Save data to file
if ~isempty(filename)
    save(filename, "X0", "x0Perturb", "X_ref", "t_ref", "tspan", "opt", '-mat')
end
    % Set output structure
data = struct("X0", X0, "x0Perturb", x0Perturb, ...
              "X_ref", X_ref, "t_ref", t_ref, "tspan", tspan, "opt", opt);

end