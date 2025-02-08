function [e, c_L, c_Di] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, N)
%   PLLT - Prandtl Lifting Line Theory Fundamental Equation Solver
%       Solves the fundamental equation of Prandtl Lifting Line Theory for
%       finite wings with thick airfoils
%       Inputs:
%           b: Span (ft)
%           a0_t: Cross-sectional lift slope at the tips (1/rad)
%           a0_r: Cross-sectional lift slope at the root (1/rad)
%           c_t: Chord at the tips (ft)
%           c_r: Chord at the root (ft)
%           aero_t: Zero-lift angle of attack at the tips (deg)
%           aero_r: Zero-lift angle of attack at the root (deg)
%           geo_t: Geometric angle of attack at the tips (deg)
%           geo_r: Geometric angle of attack at the root (deg)
%           N: Number of odd terms to include in the series expansion for
%              circulation
%       Outputs:
%           e: Span efficiency factor, AKA (1 + delta)^-1
%           c_L: Coefficient of lift
%           c_Di: Coefficient of induced drag
%
%       Author: Ian Faber, 03/14/2023

    % Define spanwise variations
    a = @(theta) a0_r + (a0_t-a0_r)*cos(theta); % 1/rad
    c = @(theta) c_r + (c_t - c_r)*cos(theta); % ft
    aero = @(theta) deg2rad(aero_r + (aero_t - aero_r)*cos(theta)); % deg -> rad
    geo = @(theta) deg2rad(geo_r + (geo_t - geo_r)*cos(theta)); % deg -> rad
    
    wing = linspace(0,pi/2,1000);
    S = 2*trapz(wing, c(wing).*(b/2).*sin(wing));
    
    % Define testing locations
    i = 1:N;
    theta = (pi/(2*N))*i';
    
    B = geo(theta) - aero(theta);
    
    column = @(theta, n) ((4*b)./(a(theta).*c(theta))).*sin(n*theta) + (n*sin(n*theta))./sin(theta);
    
    A = zeros(N,N);
    
    for k = 1:N
        A(:,k) = column(theta,2*k-1);
    end
    
    x = linsolve(A,B); % Solve for An coefficients
    
    delta = sum((2*i(2:end)'-1).*((x(2:end)./x(1)).^2));
    
    AR = b^2/S;
    
    e = (1+delta)^-1;
    c_L = x(1)*pi*AR;
    c_Di = ((c_L^2)/(pi*AR))*(1+delta);

end