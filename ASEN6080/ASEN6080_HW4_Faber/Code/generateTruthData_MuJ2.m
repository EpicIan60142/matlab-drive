function data = generateTruthData_MuJ2(pConst, oConst, x0Perturb, stations, filename, numOrbits, measNoise, dt)
% Function that generates truth orbit and measurement data based on a set
% of planetary constants, orbital parameters, and potential initial orbital
% perturbations including only dynamics from Mu and J2.
%   Inputs:
%       - pConst: Planetary constants structure as specified by
%                 getPlanetConst.m
%       - oConst: Orbital elements structure as specified by
%                 getOrbitConst.m
%       - x0Perturb: Initial orbit perturbation organized as follows:
%                    [dX; dY; dZ; dXdot; dYdot; dZdot]
%       - stations: Stations structure as defined by makeStations.m
%       - filename: Name of the file to save truth data to. 
%                   Ex: "HW2Problem1Data.mat"
%       - numOrbits: Number of orbits to generate data for
%       - measNoise: Boolean indicating whether measurement noise is
%                    included or not
%       - dt: Timestep separation for generating measurements (optional).
%             If not provided, defaults to dt = T/1000, where T is the
%             period calculated from oConst
%   Outputs:
%       - data: Data structure with the following fields:
%               - X0: Initial cartesian state based on the orbital elements
%                     in oConst, organized as [X0; Y0; Z0; Xdot0; Ydot0;
%                     Zdot0]
%               - x0Perturb: Initial state perturbation as specified by
%                            x0Perturb
%               - stations: Stations structure as defined by
%                           makeStations.m, now with propagated station 
%                           states, measurements, and measurement times
%                           according to the planetary parameters in pConst
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

planetConst = pConst;
orbital = oConst;

    % Convert orbital elements to cartesian for the initial orbit state
X0 = convOrbitalToCartesian([planetConst.mu; orbital.a; orbital.e; orbital.i; orbital.RAAN; orbital.argPeri; orbital.truAnom]);

    % Propagate reference orbit
T = 2*pi*sqrt(orbital.a^3/planetConst.mu);

if ~exist("dt", "var")
    dt = T/1000;
end

tspan = 0:dt:(numOrbits*T);
opt = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% XPhi0 = [X0+x0Perturb; reshape(eye(length(X0)),length(X0)^2,1)];

[t_ref, X_ref] = ode45(@(t,X)orbitEOM_MuJ2(t,X,planetConst.mu,planetConst.J2,planetConst.Ri), tspan, X0+x0Perturb, opt);
% [t_ref, XPhi_ref] = ode45(@(t,XPhi)STMEOM_J2(t,XPhi,planetConst.mu,planetConst.J2,planetConst.Ri),tspan,XPhi0,opt);

% X_ref = XPhi_ref(:,1:length(X0));

    % Propagate station states and generate clean measurements
for k = 1:length(t_ref)
    dTheta = planetConst.wPlanet*t_ref(k);
    for kk = 1:length(stations)
        r = rotZ(dTheta)*stations(kk).X0;
        v = cross([0;0;planetConst.wPlanet], r);

        stations(kk).Xs = [stations(kk).Xs; [r', v', t_ref(k)]];
        
        y = generateRngRngRate(X_ref(k,:), stations(kk).Xs(k,:), stations(kk).elMask);

        if ~isnan(y)
            stations(kk).rho = [stations(kk).rho; y(1)];
            stations(kk).rhoDot = [stations(kk).rhoDot; y(2)];
            stations(kk).elAngle = [stations(kk).elAngle; y(3)];
            stations(kk).tMeas = [stations(kk).tMeas; t_ref(k)];
            if measNoise
                R = diag([stations(kk).sigRho^2, stations(kk).sigRhoDot^2]);
                stations(kk).R = [stations(kk).R; {R}];
            else
                R = zeros(2,2);
                stations(kk).R = [stations(kk).R; {R}];
            end
        end
    end
end

    % Add noise to measurements based on measurement uncertainties
if measNoise
    for k = 1:length(stations)
        noise = randn(length(stations(k).rho),2);
        stations(k).rho = stations(k).rho + stations(k).sigRho*noise(:,1);
        stations(k).rhoDot = stations(k).rhoDot + stations(k).sigRhoDot*noise(:,2);
    end
end

    % Save data to file
if ~isempty(filename)
    save(filename, "X0", "x0Perturb", "stations", "X_ref", "t_ref", "tspan", "opt", '-mat')
end
    % Set output structure
data = struct("X0", X0, "x0Perturb", x0Perturb, "stations", stations, ...
              "X_ref", X_ref, "t_ref", t_ref, "tspan", tspan, "opt", opt);

end