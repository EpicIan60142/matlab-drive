function Plot_Airfoil_Flow(c, alpha, V_inf, p_inf, rho_inf, N, nLevels, bounds, nPoints)
%   Plot_Airfoil_Flow - Airfoil streamline plotting function
%       Plots the streamlines around a symmetric thin airfoil using a
%       vortex sheet with N discrete vortices
%       Inputs:
%           c: chord length (m)
%           alpha: angle of attack (deg)
%           V_inf: Freestream velocity (m/s)
%           p_inf: Freestream pressure (Pa)
%           rho_inf: Freestream density (kg/m^3)
%           N: Number of discrete vortices
%           nLevels: Number of levels for [streamlines, equipotential lines, 
%                    and pressure contours]
%           bounds: Max and Min x and y values [xmin, xmax; ymin, ymax]
%       Outputs:
%           Plot of streamlines around symmetric thin airfoil
%       
%       Author: Ian Faber, 2/23/2023

    % Define domain
    xVec = linspace(bounds(1,1), bounds(1,2), nPoints);
    yVec = linspace(bounds(2,1), bounds(2,2), nPoints);
    
    [x,y] = meshgrid(xVec, yVec);
    
    alpha = deg2rad(alpha);
    
    % Define anonymous function handles
    gamma = @(x) 2*alpha*V_inf*sqrt((1-(x/c))./(x/c));
    theta = @(x,x1,y,y1) mod(atan2(y-y1,x-x1),2*pi);
    radius = @(x,x1,y,y1) sqrt((x-x1).^2 + (y-y1).^2);

    % Define divisions
    dx = c/N;
    
    vortX = linspace(dx/2, c-(dx/2), N);
    vortX = reshape(vortX,1,1,[]);
    
    % Define elementary flows
    Psi_UF = V_inf*(y*cos(alpha) - x*sin(alpha));
    
    Psi_Vort = ((gamma(vortX)*dx)/(2*pi)).*log(radius(x, vortX, y, 0));
    Psi_Vort = sum(Psi_Vort,3);
    
    % Build full flow
    Psi = Psi_UF + Psi_Vort;
    
    % Plot flow
    levmin = min(min(Psi));
    levmax = max(max(Psi));
    levels = linspace(levmin,levmax,nLevels(1))';
    
    figure(1)
    hold on
    axis equal
    titleText = sprintf("Streamlines at c = %.0f m, \\alpha = %.0f^o, and V_{\\infty} = %.0f m/s", c, rad2deg(alpha), V_inf);
    title(titleText)
    contourf(x, y, Psi, levels, 'LineWidth', 1.5)
    colormap('cool')
    colorbar
    
    % Plot Body
    plot(linspace(0,c,nPoints), zeros(length(xVec)),'w-', 'LineWidth', 4)
    
    % Define elementary potentials
    Phi_UF = V_inf*(y*sin(alpha) + x*cos(alpha));
    
    thetaVort = theta(x,vortX,y,0);
    Phi_Vort = -((gamma(vortX)*dx)/(2*pi)).*thetaVort;
    Phi_Vort = sum(Phi_Vort,3);
    
    % Build full potential
    Phi = Phi_UF + Phi_Vort;
    
    % Plot potential
    levmin2 = min(min(Phi));
    levmax2 = max(max(Phi));
    levels2 = linspace(levmin2,levmax2,nLevels(2))';
    
    figure(2)
    hold on
    axis equal
    titleText = sprintf("Equipotential lines at c = %.0f m, \\alpha = %.0f^o, and V_{\\infty} = %.0f m/s", c, rad2deg(alpha), V_inf);
    title(titleText)
    contourf(x,y,Phi,levels2,'LineWidth',1.5)
    colormap('cool')
    colorbar
    
    % Plot body
    plot(linspace(0,c,nPoints), zeros(length(xVec)),'w-', 'LineWidth', 4)
    
    % Define Cp
    [u, v] = gradient(Phi, mean(diff(xVec)), mean(diff(yVec)));
    % Limit vertical speeds to 100 m/s due to discontinuities in the field
    v(v>100) = 100;
    v(v<-100) = -100;
    V = sqrt(u.^2 + v.^2);
    Cp = 1 - (V/V_inf).^2;
    Cp(Cp < -3) = -3; % Limit Cp to -3 for better contour coloring

    % Calculate pressure
    p = p_inf + (0.5*rho_inf*V_inf^2)*Cp;
     
    % Plot pressure
    levmin3 = min(min(p));
    levmax3 = max(max(p));
    levels3 = linspace(levmin3,levmax3,nLevels(3))';
    
    figure(3)
    hold on
    axis equal
    titleText = sprintf("Pressure contours at c = %.0f m, \\alpha = %.0f^o, and V_{\\infty} = %.0f m/s", c, rad2deg(alpha), V_inf);
    title(titleText)
    contour(x,y,p,levels3,'LineWidth',1.5)
    set(gca,'Color','k')
    colormap('cool')
    colorbar
    
    % Plot body
    plot(linspace(0,c,nPoints), zeros(length(xVec)),'w-', 'LineWidth', 4)

end



