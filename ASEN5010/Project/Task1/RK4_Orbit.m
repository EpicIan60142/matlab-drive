function out = RK4_Orbit(x0, t0, dt, tf)
% Function that implements the Runga-Kutta 4 algorithm to integrate 
% circular orbital motion based on a set of initial conditions
%   Inputs:
%       - x0: Initial state vector, with w written in orbit coordinates
%               [radius; EA_0; w_0]
%       - t0: Time that integration will start, in seconds
%       - dt: Time step for integration, in seconds
%       - tf: Time that integration will stop, in seconds
%
%   Outputs:
%       - out: Integration output matrix, each column is a vector with the
%              same number of elements n as there were timesteps
%               [t (nx1), r (nx3), rDot (nx3), EA (nx3), w (nx3)]
%
    radius = x0(1);
    EA_0 = x0(2:4);
    w_0 = x0(5:7);

    X = [EA_0; w_0];
    t = t0;
    
    [~, r, rDot] = calculateOrbit(X, radius);

    out = zeros(length(t0:dt:tf)-1, 13);
    out(1,:) = [t0, r', rDot', X']; % t, r(1:3), rDot(1:3), EA(1:3), w(1:3)
    k = 1;

    while t < tf
        k1 = dt*calculateOrbit(X,radius);
        k2 = dt*calculateOrbit(X+k1/2,radius);
        k3 = dt*calculateOrbit(X+k2/2,radius);
        k4 = dt*calculateOrbit(X+k3,radius);

        X = X + (1/6)*(k1 + 2*k2 + 2*k3 + k4);

        t = t + dt;
        k = k + 1;

        [~, r, rDot] = calculateOrbit(X, radius);

        out(k, :) = [t, r', rDot', X'];

    end

end