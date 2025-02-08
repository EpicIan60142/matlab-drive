function dX = suspensionEOM(t,X,const)
% EOM function for simulating a simplified car suspension with ode45
%   Inputs:
%       t: time [sec]
%       X: state vector
%           [ x; y; vx; vy ]
%       const: vector of constants for simulation
%           [m1; m2; Kw; Ks; b; r; disturbTime]
%
%   Outputs:
%       dX: rate of change vector
%           [ vx; vy; ax; ay ]
%
%   By: Ian Faber, 09/05/2023
%

    x = X(1);
    y = X(2);
    vx = X(3);
    vy = X(4);
    
    m1 = const(1);
    m2 = const(2);
    Kw = const(3);
    Ks = const(4);
    b = const(5);
    r = const(6);
    disturbTime = const(7);
    
    if t < disturbTime % Don't apply the disturbance r until specified in the main script
        r = 0;
    end
    
    ax = (Kw/m1)*(r-x) + (Ks/m1)*(y-x) + (b/m1)*(vy-vx);
    ay = (Ks/m2)*(x-y) + (b/m2)*(vx-vy);
    
    dX = [vx; vy; ax; ay];

end