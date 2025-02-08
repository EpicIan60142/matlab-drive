%% ASEN 1320 Final Project - MATLAB
clc;clear;close all;

Sat1Init = [1986.2;6388.3;-1237.2;-4.93;0.402;-5.831]; %Set the initial state for the ISS
Sat2Init = [6480.8;1108.2;-2145.5;-0.29;7.0712;2.747]; %Set the initial state for the Hubble Space Telescope
simulationInterval = 0:60:86400; %Set the simulation duration as 1 day in terms of 60 second intervals

EOMFun = @OrbitEOM; %Create a function handle for the equations of motion
setup = odeset('RelTol',1e-12,'AbsTol',1e-12); %Create a setup variable for ode45
[Sat1Time, Sat1State] = ode45(EOMFun, simulationInterval, Sat1Init, setup); %Calculate the trajectory for ISS
[Sat2Time, Sat2State] = ode45(EOMFun, simulationInterval, Sat2Init, setup); %Calculate the trajectory for Hubble

Sat1Position = Sat1State(:,1:3); %Extract the coordinates of ISS from the calculated trajectory
Sat2Position = Sat2State(:,1:3); %Extract the coordinates of Hubble from the calculated trajectory
GSPosition = load('GSPosition.csv','-ascii');
Sat1Visibility = load('Sat1Visibility.csv','-ascii');
Sat2Visibility = load('Sat2Visibility.csv','-ascii');

csvwrite('Sat1Position.csv',Sat1Position); %Save the ISS's position data into a csv file
csvwrite('Sat2Position.csv',Sat2Position); %Save the Hubble's position data into a csv file

angle = 0;

figh = figure;

Sat1VisibilityPlot(1441,3) = [0];
Sat2VisibilityPlot(1441,3) = [0];

for k = 1:length(simulationInterval)
    if(Sat1Visibility(k) == 1)
        Sat1VisibilityPlot(k,:) = Sat1Position(k,:);
    else
        Sat1VisibilityPlot(k,:) = " ";
    end
end

for k = 1:length(simulationInterval)
    if(Sat2Visibility(k) == 1)
        Sat2VisibilityPlot(k,:) = Sat2Position(k,:);
    else
        Sat2VisibilityPlot(k,:) = " ";
    end
end

for k = 1:length(simulationInterval)
    clf
    
    hold on
    grid on

    ISSTrajectory = plot3(Sat1Position(:,1),Sat1Position(:,2),Sat1Position(:,3),'r'); %Plot the ISS's orbit in 3D
    HubbleTrajectory = plot3(Sat2Position(:,1),Sat2Position(:,2),Sat2Position(:,3),'b'); %Plot the Hubble's orbit in 3D
    GroundStationTrajectory = plot3(GSPosition(:,1),GSPosition(:,2),GSPosition(:,3),'g');
    
    Sat1VisibilityGraphic = plot3(Sat1VisibilityPlot(:,1),Sat1VisibilityPlot(:,2),Sat1VisibilityPlot(:,3),'k.','MarkerSize',25);
    Sat2VisibilityGraphic = plot3(Sat2VisibilityPlot(:,1),Sat2VisibilityPlot(:,2),Sat2VisibilityPlot(:,3),'k.','MarkerSize',25);
    
    t_k = simulationInterval(k);
    
    %Satellite 1 Animation
    x1_k = Sat1Position(k,1);
    y1_k = Sat1Position(k,2);
    z1_k = Sat1Position(k,3);
    if(Sat1Visibility(k) > 0)
        Sat1 = plot3(x1_k,y1_k,z1_k,'m.','LineWidth',10,'MarkerSize',25);
    else
        Sat1 = plot3(x1_k,y1_k,z1_k, 'r.', 'LineWidth', 3, 'MarkerSize', 15);
    end
    
    %Satellite 2 Animation
    x2_k = Sat2Position(k,1);
    y2_k = Sat2Position(k,2);
    z2_k = Sat2Position(k,3);
    if(Sat2Visibility(k) > 0)
        Sat2 = plot3(x2_k,y2_k,z2_k,'m.','LineWidth',10,'MarkerSize',25);
    else
        Sat2 = plot3(x2_k,y2_k,z2_k, 'b.', 'LineWidth', 3, 'MarkerSize', 15);
    end
    
    %Ground Station Animation
    x3_k = GSPosition(k,1);
    y3_k = GSPosition(k,2);
    z3_k = GSPosition(k,3);
    GroundStation = plot3(x3_k,y3_k,z3_k, 'g.', 'LineWidth', 3, 'MarkerSize', 15);
    
    %Earth Animation
    angle = angle + 60*(360/86400)*(pi/180);
    [earthX, earthY, earthZ] = sphere;
    earthXrot = earthX*cos(angle)-earthY*sin(angle);
    earthYrot = earthX*sin(angle)+earthY*cos(angle);
    earthZ = earthZ;
    earth = surf(6370*earthXrot, 6370*earthYrot, 6370*earthZ);
    I = imread('EarthMap.jpeg');
    %I = flipImage(I);
    set(earth,'FaceColor','texturemap','cdata',I,'edgecolor','none');
    
    subset = [ISSTrajectory, HubbleTrajectory, GroundStationTrajectory,Sat1,Sat2,GroundStation];
    
    %Create a title, label the axes, and create a legend for the orbit visualization
    title(['Orbits of the ISS and Hubble Space Telescope at t = ', num2str(t_k),' seconds']);
    xlabel('X Axis (km)');
    ylabel('Y Axis (km)');
    zlabel('Z Axis (km)');
    view([30 35]);
    legend(subset,'ISS Trajectory','Hubble Telescope Trajectory','Ground Station Path','ISS','Hubble Telescope','Ground Station','Location','southOutside','Orientation','vertical');
    
    drawnow
    %movieVector(k) = getframe(figh, [0 0 560 420]);
end



% movie = VideoWriter('FinalProjectFullMovie','MPEG-4');
% movie.FrameRate = 60;
% 
% %Open the VideoWriter object, write the movie, and close the file
% open(movie);
% writeVideo(movie, movieVector);
% close(movie);
