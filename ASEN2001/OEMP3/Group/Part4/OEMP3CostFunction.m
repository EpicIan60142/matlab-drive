%% OEMP 3 Group 2 Part 4
% Group 11, 10:40 Lab
%
% Code by: Nathan Evans

clear; clc; close all;

%% Constants
w0 = 2001/12; % lb/in
L = 27.25*12; % in

maxDim = 1*12; % in

%density and sigma
rhovec = [.098, .283, .304, .284, .16];
sigmaYieldvec = [35, 70, 35, 115, 120] .* 1000;
price = [8.03, 8.07, 52.78, 29.63, 115.36]; 
%% Hollow Square
for i = 1:5
    
    rho = rhovec(i);
    sigmaYield = sigmaYieldvec(i);
    moment = @(A, x) (w0/2)*(((-x^3)/(3*L)) + x^2 - L*x + (L^2)/2) - (rho* A * (L-x)^2)/2;

    inertiaCircle = @(r) (pi*r^4)/4;
    inertiaRectangle = @(b, h) (1/12)*b*h^3;

    bendingStress = @(M, y, I) (M*y)/I;

    factorOfSafety = @(sigmaYield, sigmaApplied) sigmaYield/sigmaApplied;

    maxX = @(A) 0;%-L*(2*A*rho - w0)/w0;

    factors.shape = [];
    factors.area = [];
    factors.factorOfSafety = [];
for k = 1:4*maxDim % 1/4 in to 12 in length
    for l = 1:4*maxDim % 1/4 in to 12 in height
        for m = 1:4*(maxDim/2) % 1/4 in to 6 in wall thickness
            b = k/4;
            h = l/4;
            t = m/4;
            
            if(t >= b/2 || t >= h/2)
                break;
            end
            
            b1 = t;
            b2 = b - (2*t);
            h1 = h;
            h2 = t;
            
            area1 = b1*h1;
            area2 = b2*h2;
            area = 2*area1 + 2*area2;
            
            inertia = 2*(inertiaRectangle(b1, h1) + area1*(b-(t/2))^2 + inertiaRectangle(b2, h2) + area2*(h-(t/2))^2);
            
            x = maxX(area);
            
            sigmaBend = bendingStress(moment(area, x), h/2, inertia);
            factorHollow = factorOfSafety(sigmaYield, sigmaBend);
            
            if(factorHollow >= 1.5 && factorHollow <= 1.53)
                description = sprintf("hollow square, b = %.3f, h = %.3f, t = %0.3f", b, h, t);
                factors.shape = [factors.shape; description];
                factors.area = [factors.area; area];
                factors.factorOfSafety = [factors.factorOfSafety; factorHollow];
            end
        end
    end
end
hollow = factors;
hollowBeamResult = struct2table(hollow);
hollowBest = find(hollow.area == min(hollow.area));
hollowFinal = hollowBeamResult(hollowBest,:)
volume = hollowFinal.area * L;
weight = volume * rho;
cost = weight * price(i)
end



