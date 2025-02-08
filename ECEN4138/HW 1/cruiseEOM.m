function dX = cruiseEOM(t,X,const)
% EOM function for simulating a cruise control system with ode45
%   Inputs:
%       t: time [sec]
%       X: state vector
%           v
%       const: vector of constants for simulation
%           [m; b; K; vr; disturbTime]
%
%   Outputs:
%       dX: rate of change vector
%           [ vx; vy; ax; ay ]
%
%   By: Ian Faber, 09/05/2023
%

    v = X;

    m = const(1);
    b = const(2);
    K = const(3);
    vr = const(4);
    disturbTime = const(5);

    if t < disturbTime
        vr = 0;
    end

    a = (1/m)*(K*vr - (b+K)*v);

    dX = a;

end