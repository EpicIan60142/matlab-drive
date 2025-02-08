function Sensitivity_Analysis(c, alpha, V_inf, p_inf, rho_inf, N, nLevels, bounds, nPoints)
%   Sensitivity_Analysis - Function that tests different values of c,
%       alpha, and V_inf to see how the streamlines and equipotential lines 
%       are affected
%       Inputs:
%           c: vector of chord lengths to check (1xn)
%           alpha: vector of angles of attack to check (1xn)
%           V_inf: vector of freestream velocities to check (1xn)
%           p_inf: Freestream pressure (Pa)
%           rho_inf: Freestream density (kg/m^3)
%           N: Number of discrete vortices
%           nLevels: Number of levels for [streamlines, equipotential lines, 
%                    and pressure contours]
%           bounds: Max and Min x and y values [xmin, xmax; ymin, ymax]
%       Outputs:
%           Plot of streamlines around symmetric thin airfoil
%
%       Author: Ian Faber, 2/28/2023

    % Define domain
    xVec = linspace(bounds(1,1), bounds(1,2), nPoints);
    yVec = linspace(bounds(2,1), bounds(2,2), nPoints);
    
    [x,y] = meshgrid(xVec, yVec);

    alpha = deg2rad(alpha);

    % Define anonymous function handles
    gamma = @(x,ai,vi,ci) 2*alpha(ai)*V_inf(vi)*sqrt((1-(x/c(ci)))./(x/c(ci)));
    theta = @(x,x1,y,y1) mod(atan2(y-y1,x-x1),2*pi);
    radius = @(x,x1,y,y1) sqrt((x-x1).^2 + (y-y1).^2);

    % Initial/Base value indices
    cIdx = 5;
    aIdx = 3;
    VIdx = 6;
    
    % Test chord length and create a tiled layouts for each contour plot
    figure
    SLc = tiledlayout(3,2);
    title(SLc,"Streamline Variation with c")
    xlabel(SLc, "x [m]")
    ylabel(SLc, "y [m]")
    
    figure
    EPc = tiledlayout(3,2);
    title(EPc,"Equipotential Line Variation with c")
    xlabel(EPc, "x [m]")
    ylabel(EPc, "y [m]")
    
    figure
    PCc = tiledlayout(3,2);
    title(PCc,"Pressure Contour Variation with c")
    xlabel(PCc, "x [m]")
    ylabel(PCc, "y [m]")

    for k = 1:length(c)
        
        dx = c(k)/N;
    
        vortX = linspace(dx/2, c(k)-(dx/2), N);
        vortX = reshape(vortX,1,1,[]);

        % Define elementary flows
        Psi_UF = V_inf(VIdx)*(y*cos(alpha(aIdx)) - x*sin(alpha(aIdx)));
        
        Psi_Vort = ((gamma(vortX,aIdx,VIdx,k)*dx)/(2*pi)).*log(radius(x, vortX, y, 0));
        Psi_Vort = sum(Psi_Vort,3);
        
        % Build full flow
        Psi = Psi_UF + Psi_Vort;
        
        % Plot flow
        levmin = min(min(Psi));
        levmax = max(max(Psi));
        levels = linspace(levmin,levmax,nLevels(1))';
        
        nexttile(SLc)
        hold on
        titleText = sprintf("c = %.0f m", c(k));
        title(titleText)
        contour(x, y, Psi, levels, 'LineWidth', 1.5)
        colormap('cool')
        colorbar
        set(gca,'Color','k')
        
        % Plot Body
        plot(linspace(0,c(k),nPoints), zeros(length(xVec)),'w-', 'LineWidth', 3)

        xlim(bounds(1,:))
        ylim(bounds(2,:))
        
        % Define elementary potentials
        Phi_UF = V_inf(VIdx)*(y*sin(alpha(aIdx)) + x*cos(alpha(aIdx)));
        
        thetaVort = theta(x,vortX,y,0);
        Phi_Vort = -((gamma(vortX,aIdx,VIdx,k)*dx)/(2*pi)).*thetaVort;
        Phi_Vort = sum(Phi_Vort,3);
        
        % Build full potential
        Phi = Phi_UF + Phi_Vort;
        
        % Plot potential
        levmin2 = min(min(Phi));
        levmax2 = max(max(Phi));
        levels2 = linspace(levmin2,levmax2,nLevels(2))';
        
        nexttile(EPc)
        hold on
        titleText = sprintf("c = %.0f m", c(k));
        title(titleText)
        contour(x,y,Phi,levels2,'LineWidth',1.5)
        colormap('cool')
        colorbar
        set(gca,'Color','k')
        
        % Plot body
        plot(linspace(0,c(k),nPoints), zeros(length(xVec)),'w-', 'LineWidth', 3)

        xlim(bounds(1,:))
        ylim(bounds(2,:))
        
        % Define Cp
        [u, v] = gradient(Phi, mean(diff(xVec)), mean(diff(yVec)));
        v(v>100) = 100;
        v(v<-100) = -100;
        V = sqrt(u.^2 + v.^2);
        Cp = 1 - (V/V_inf(VIdx)).^2;
        Cp(Cp < -3) = -3;
        
        % Calculate pressure
        p = p_inf + (0.5*rho_inf*V_inf(VIdx)^2)*Cp;
        
        % Plot pressure
        levmin3 = min(min(p));
        levmax3 = max(max(p));
        levels3 = linspace(levmin3,levmax3,nLevels(3))';
        
        nexttile(PCc)
        hold on
        titleText = sprintf("c = %.0f m", c(k));
        title(titleText)
        contour(x,y,p,levels3,'LineWidth',1.5)
        colormap('cool')
        colorbar
        set(gca,'Color','k')
        
        % Plot body
        plot(linspace(0,c(k),nPoints), zeros(length(xVec)),'w-', 'LineWidth', 3)

        xlim(bounds(1,:))
        ylim(bounds(2,:))
    end

    % Test alpha
    figure
    SLa = tiledlayout(3,2);
    title(SLa,"Streamline Variation with \alpha")
    xlabel(SLa, "x [m]")
    ylabel(SLa, "y [m]")
    
    figure
    EPa = tiledlayout(3,2);
    title(EPa,"Equipotential Line Variation with \alpha")
    xlabel(EPa, "x [m]")
    ylabel(EPa, "y [m]")
    
    figure
    PCa = tiledlayout(3,2);
    title(PCa,"Pressure Contour Variation with \alpha")
    xlabel(PCa, "x [m]")
    ylabel(PCa, "y [m]")

    for k = 1:length(alpha)
        dx = c(cIdx)/N;
    
        vortX = linspace(dx/2, c(cIdx)-(dx/2), N);
        vortX = reshape(vortX,1,1,[]);

        % Define elementary flows
        Psi_UF = V_inf(VIdx)*(y*cos(alpha(k)) - x*sin(alpha(k)));
        
        Psi_Vort = ((gamma(vortX,k,VIdx,cIdx)*dx)/(2*pi)).*log(radius(x, vortX, y, 0));
        Psi_Vort = sum(Psi_Vort,3);
        
        % Build full flow
        Psi = Psi_UF + Psi_Vort;
        
        % Plot flow
        levmin = min(min(Psi));
        levmax = max(max(Psi));
        levels = linspace(levmin,levmax,nLevels(1))';
        
        nexttile(SLa)
        hold on
        titleText = sprintf("\\alpha = %.0f\\circ", rad2deg(alpha(k)));
        title(titleText)
        contour(x, y, Psi, levels, 'LineWidth', 1.5)
        colormap('cool')
        colorbar
        set(gca,'Color','k')
        
        % Plot Body
        plot(linspace(0,c(cIdx),nPoints), zeros(length(xVec)),'w-', 'LineWidth', 3)

        xlim(bounds(1,:))
        ylim(bounds(2,:))
        
        % Define elementary potentials
        Phi_UF = V_inf(VIdx)*(y*sin(alpha(k)) + x*cos(alpha(k)));
        
        thetaVort = theta(x,vortX,y,0);
        Phi_Vort = -((gamma(vortX,k,VIdx,cIdx)*dx)/(2*pi)).*thetaVort;
        Phi_Vort = sum(Phi_Vort,3);
        
        % Build full potential
        Phi = Phi_UF + Phi_Vort;
        
        % Plot potential
        levmin2 = min(min(Phi));
        levmax2 = max(max(Phi));
        levels2 = linspace(levmin2,levmax2,nLevels(2))';
        
        nexttile(EPa)
        hold on
        titleText = sprintf("\\alpha = %.0f\\circ", rad2deg(alpha(k)));
        title(titleText)
        contour(x,y,Phi,levels2,'LineWidth',1.5)
        colormap('cool')
        colorbar
        set(gca,'Color','k')
        
        % Plot body
        plot(linspace(0,c(cIdx),nPoints), zeros(length(xVec)),'w-', 'LineWidth', 3)

        xlim(bounds(1,:))
        ylim(bounds(2,:))
        
        % Define Cp
        [u, v] = gradient(Phi, mean(diff(xVec)), mean(diff(yVec)));
        v(v>100) = 100;
        v(v<-100) = -100;
        V = sqrt(u.^2 + v.^2);
        Cp = 1 - (V/V_inf(VIdx)).^2;
        Cp(Cp < -3) = -3;
        
        % Calculate pressure
        p = p_inf + (0.5*rho_inf*V_inf(VIdx)^2)*Cp;
        
        % Plot pressure
        levmin3 = min(min(p));
        levmax3 = max(max(p));
        levels3 = linspace(levmin3,levmax3,nLevels(3))';
        
        nexttile(PCa)
        hold on
        titleText = sprintf("\\alpha = %.0f\\circ", rad2deg(alpha(k)));
        title(titleText)
        contour(x,y,p,levels3,'LineWidth',1.5)
        colormap('cool')
        colorbar
        set(gca,'Color','k')
        
        % Plot body
        plot(linspace(0,c(cIdx),nPoints), zeros(length(xVec)),'w-', 'LineWidth', 3)

        xlim(bounds(1,:))
        ylim(bounds(2,:))
    end

    % Test V_inf
    figure
    SLv = tiledlayout(3,2);
    title(SLv,"Streamline Variation with V_{\infty}")
    xlabel(SLa, "x [m]")
    ylabel(SLa, "y [m]")
    
    figure
    EPv = tiledlayout(3,2);
    title(EPv,"Equipotential Line Variation with V_{\infty}")
    xlabel(EPa, "x [m]")
    ylabel(EPa, "y [m]")
    
    figure
    PCv = tiledlayout(3,2);
    title(PCv,"Pressure Contour Variation with V_{\infty}")
    xlabel(PCa, "x [m]")
    ylabel(PCa, "y [m]")

    for k = 1:length(V_inf)
        dx = c(cIdx)/N;
    
        vortX = linspace(dx/2, c(cIdx)-(dx/2), N);
        vortX = reshape(vortX,1,1,[]);

        % Define elementary flows
        Psi_UF = V_inf(k)*(y*cos(alpha(aIdx)) - x*sin(alpha(aIdx)));
        
        Psi_Vort = ((gamma(vortX,aIdx,k,cIdx)*dx)/(2*pi)).*log(radius(x, vortX, y, 0));
        Psi_Vort = sum(Psi_Vort,3);
        
        % Build full flow
        Psi = Psi_UF + Psi_Vort;
        
        % Plot flow
        levmin = min(min(Psi));
        levmax = max(max(Psi));
        levels = linspace(levmin,levmax,nLevels(1))';
        
        nexttile(SLv)
        hold on
        titleText = sprintf("V_{\\infty} = %.0f m/s", V_inf(k));
        title(titleText)
        contour(x, y, Psi, levels, 'LineWidth', 1.5)
        colormap('cool')
        colorbar
        set(gca,'Color','k')
        
        % Plot Body
        plot(linspace(0,c(cIdx),nPoints), zeros(length(xVec)),'w-', 'LineWidth', 3)

        xlim(bounds(1,:))
        ylim(bounds(2,:))
        
        % Define elementary potentials
        Phi_UF = V_inf(k)*(y*sin(alpha(aIdx)) + x*cos(alpha(aIdx)));
        
        thetaVort = theta(x,vortX,y,0);
        Phi_Vort = -((gamma(vortX,aIdx,k,cIdx)*dx)/(2*pi)).*thetaVort;
        Phi_Vort = sum(Phi_Vort,3);
        
        % Build full potential
        Phi = Phi_UF + Phi_Vort;
        
        % Plot potential
        levmin2 = min(min(Phi));
        levmax2 = max(max(Phi));
        levels2 = linspace(levmin2,levmax2,nLevels(2))';
        
        nexttile(EPv)
        hold on
        titleText = sprintf("V_{\\infty} = %.0f m/s", V_inf(k));
        title(titleText)
        contour(x,y,Phi,levels2,'LineWidth',1.5)
        colormap('cool')
        colorbar
        set(gca,'Color','k')
        
        % Plot body
        plot(linspace(0,c(cIdx),nPoints), zeros(length(xVec)),'w-', 'LineWidth', 3)

        xlim(bounds(1,:))
        ylim(bounds(2,:))
        
        % Define Cp
        [u, v] = gradient(Phi, mean(diff(xVec)), mean(diff(yVec)));
        v(v>100) = 100;
        v(v<-100) = -100;
        V = sqrt(u.^2 + v.^2);
        Cp = 1 - (V/V_inf(k)).^2;
        Cp(Cp < -3) = -3;
        
        % Calculate pressure
        p = p_inf + (0.5*rho_inf*V_inf(k)^2)*Cp;
        
        % Plot pressure
        levmin3 = min(min(p));
        levmax3 = max(max(p));
        levels3 = linspace(levmin3,levmax3,nLevels(3))';
        
        nexttile(PCv)
        hold on
        titleText = sprintf("V_{\\infty} = %.0f m/s", V_inf(k));
        title(titleText)
        contour(x,y,p,levels3,'LineWidth',1.5)
        colormap('cool')
        colorbar
        set(gca,'Color','k')
        
        % Plot body
        plot(linspace(0,c(cIdx),nPoints), zeros(length(xVec)),'w-', 'LineWidth', 3)

        xlim(bounds(1,:))
        ylim(bounds(2,:))
    end

end