%% Problem 2

% Bullet 2

    % All units are in meters and Newtons and values are taken from ANSYS
    % data.

    % Reading in the data
    ReactionForces = load("Data\ReactionForces1.txt");
    BarForces = load("Data\BarForces1.txt");
    
    % Ratio needed for incorrect calculations in ANSYS
    ratio = 2.5/(3.959e-5);
    
    % The order for these are "node","force value"
    ReactF = [1, 55.562; 52, 55.638; 17, 55.638; 68, 55.562];
    BarF = [128, -440.02; 131, -440.02; 178, -388.15; 180, -388.22];
    
    % The order for this is "node","displacement"
    NodeDisp = [26,ratio*(-0.29207e-007);43,ratio*(-0.29263e-007)];
    
    % Info for Figures
    BarFfront = linspace(0,1,10)*351.233 - 400.34; % force in bar on bottom front of truss
    dispMiddle = linspace(0,1,10)*(ratio*(-0.29207e-007)); % Displacement of nodes in the middle of truss  
    extLoad = linspace(0,1,10)*(222.4); % external load value (total)
    ReactF1 = linspace(0,1,10)*55.562; % reaction force at node 1 (pin)
    ReactF52 = linspace(0,1,10)*55.638; % reaction force at node 52 (pin)
    ReactFRoller = linspace(0,1,10)*(55.638+55.562);
    
    % Plots
        figure(1)
        hold on
        grid on
        plot(extLoad,ReactF1);
        xlabel("External Load(N)")
        ylabel("Reaction Forces(N)")
        title("Reaction Forces vs External Load (Joint 1)")
        
        figure(2)
        hold on
        grid on
        plot(extLoad,ReactF52);
        xlabel("External Load(N)")
        ylabel("Reaction Forces(N)")
        title("Reaction Forces vs External Load (Joint 52)")
        
        figure(3)
        hold on
        grid on
        plot(extLoad,ReactFRoller);
        xlabel("External Load(N)")
        ylabel("Reaction Forces(N)")
        title("Reaction Forces vs External Load (Roller)")
        
        figure(4)
        hold on
        grid on
        plot(extLoad,BarFfront);
        xlabel("External Load (N)")
        ylabel("Bar Force (N)")
        title("Bar Forces vs External Load")
        
        figure(5)
        hold on 
        grid on
        plot(extLoad,dispMiddle)
        xlabel("External Load (N)")
        ylabel("Displacement (m)")
        title("Displacement of Truss Center vs External Load")
    
    
    