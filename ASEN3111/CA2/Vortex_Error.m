function Vortex_Error(c, alpha, V_inf, NMax, bounds, nPoints)
%   Vortex_Error - Function that plots the velocity and pressure error as a
%       function of the number of discrete vortices, N
%       Inputs:
%           c: chord length (m)
%           alpha: angle of attack (deg)
%           V_inf: Freestream velocity (m/s)
%           NMax: maximum number of vortices to test
%       Outputs:
%           Plots of velocity and pressure error as a function of N
%       
%       Author: Ian Faber, 2/25/2023
    
    % Define domain
    xVec = linspace(bounds(1,1), bounds(1,2), nPoints);
    yVec = linspace(bounds(2,1), bounds(2,2), nPoints);
    
    [x,y] = meshgrid(xVec, yVec);

    alpha = deg2rad(alpha);


    % Define anonymous function handles
    gamma = @(x) 2*alpha*V_inf*sqrt((1-(x/c))./(x/c));
    theta = @(x,x1,y,y1) mod(atan2(y-y1,x-x1),2*pi);
%     radius = @(x,x1,y,y1) sqrt((x-x1).^2 + (y-y1).^2);

    % Define the number of vortices for the "true" value to compare against
    NTrue = NMax*2;

    % Define divisions (true solution)
    dxTrue = c/NTrue;

    vortXTrue = linspace(dxTrue/2, c-(dxTrue/2), NTrue);
    vortXTrue = reshape(vortXTrue,1,1,[]);

    % Define elementary potentials (true solution)
    Phi_UF = V_inf*(y*sin(alpha) + x*cos(alpha));
    
    thetaVortTrue = theta(x,vortXTrue,y,0);
    Phi_VortTrue = -((gamma(vortXTrue)*dxTrue)/(2*pi)).*thetaVortTrue;
    Phi_VortTrue = sum(Phi_VortTrue,3);

    % Build full potential (true solution)
    Phi_True = Phi_UF + Phi_VortTrue;

    % Define V and Cp (true solution)
    [u, v] = gradient(Phi_True, mean(diff(xVec)), mean(diff(yVec)));
    v(v>100) = 100;
    v(v<-100) = -100;
    V_True = sqrt(u.^2 + v.^2);
    Cp_True = 1 - (V_True/V_inf).^2;
    Cp_True(Cp_True < -3) = -3;

    % "Truth" test values
    trueValV = mean(Phi_True, 'all');
    trueValCp = mean(Cp_True, 'all');

    % Preallocate error test
    Phi = zeros(size([x,y]));
    Cp = zeros(size([x,y]));

    testValV = zeros(NMax,1);
    testValCp = zeros(NMax,1);
    
    % Run error test
    for N = 1:NMax

        dx = c/N;
    
        vortX = linspace(dx/2, c-(dx/2), N);
        vortX = reshape(vortX,1,1,[]);

        thetaVort = theta(x,vortX,y,0);
        Phi_Vort = -((gamma(vortX)*dx)/(2*pi)).*thetaVort;
        Phi_Vort = sum(Phi_Vort,3);
        
        % Build full potential
        Phi = Phi_UF + Phi_Vort;

        [u, v] = gradient(Phi, mean(diff(xVec)), mean(diff(yVec)));
        v(v>100) = 100;
        v(v<-100) = -100;
        V = sqrt(u.^2 + v.^2);
        Cp = 1 - (V/V_inf).^2;
        Cp(Cp < -3) = -3;

        testValV(N) = mean(Phi,'all');
        testValCp(N) = mean(Cp,'all');
    end

    % Define error vector
    VError = abs((testValV - trueValV)/trueValV);
    CpError = abs((testValCp - trueValCp)/trueValCp);

    % Find indices for 3% error
    for k = 1:length(testValV)
        if VError(k) < 0.03
            nV = k;
            break;
        end
    end

    for k = 1:length(testValCp)
        if CpError(k) < 0.03
            nCp = k;
            break;
        end
    end

    % Plot error and display indices
    figure
    hold on
    title("Velocity Error")
    plot(VError);
    yline(0.03, 'k--')
    xlabel("Number of vortices")
    ylabel("Velocity error")
    legend("Error", "3% Error Threshold")

    figure
    hold on
    title("Cp Error")
    plot(CpError);
    yline(0.03, 'k--')
    xlabel("Number of vortices")
    ylabel("Cp error")
    legend("Error", "3% Error Threshold")

    fprintf("Number of vortices for 3%% velocity error: %.0f\n", nV);
    fprintf("Number of vortices for 3%% Cp error: %.0f\n", nCp);
end