function dX = pendulumEOM(t,X,g,config)
% EOM function for simulating a simple pendulum with ode45
%   Inputs:
%       t: time [sec]
%       X: state vector
%           [ theta; omega ]
%       config: Type of EOM to simulate
%           0 = Linear, 1 = Nonlinear, defaults to linear if not 0 or 1
%
%   Outputs:
%       dX: rate of change vector
%           [ omega; alpha ]
%
%   By: Ian Faber, 09/05/2023
%

    % Extract state variables
    theta = X(1);
    omega = X(2);
    
    % Choose equation set to simulate
    switch config
        case 0 % Linear
            alpha = -g*theta;
        case 1 % Nonlinear
            alpha = -g*sin(theta);
        otherwise % Dumb user
            alpha = -g*theta;
    end
    
    dX = [omega; alpha];

end