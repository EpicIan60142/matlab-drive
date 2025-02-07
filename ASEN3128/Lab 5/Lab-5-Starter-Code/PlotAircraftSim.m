% Problem 2c Modified Plotting function

function figs = PlotAircraftSim(TOUT, aircraft_state, control_surfaces, background_wind_array, col)


%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
 
subplot(311);
h1= plot(TOUT, aircraft_state(:,1),col);hold on;
title('Position v Time');   
ylabel('X [m]')    

subplot(312);
plot(TOUT, aircraft_state(:,2),col); hold on;
ylabel('Y [m]')    
 
figs(1) = subplot(313);
plot(TOUT, aircraft_state(:,3),col);hold on;
ylabel('Z [m]')    
xlabel('time [sec]');

%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
 
subplot(311);
plot(TOUT, (180/pi)*aircraft_state(:,4),col);hold on;
title('Euler Angles v Time');   
ylabel('Roll [deg]')    

subplot(312);
plot(TOUT, (180/pi)*aircraft_state(:,5),col);hold on;
 ylabel('Pitch [deg]')    
 
figs(2) = subplot(313);
plot(TOUT, (180/pi)*aircraft_state(:,6),col);hold on;
ylabel('Yaw [deg]')    
xlabel('time [sec]');

%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
 
subplot(311);
plot(TOUT, aircraft_state(:,7),col);hold on;
title('Velocity v Time');   
ylabel('uE [m/s]')    

subplot(312);
plot(TOUT, aircraft_state(:,8),col);hold on;
 ylabel('vE [m/s]')    
 
figs(3) = subplot(313);
plot(TOUT, aircraft_state(:,9),col);hold on;
ylabel('wE [m/s]')    
xlabel('time [sec]');

%%%%%%%%%%%%%%%%%%%%%%%%
figure(4);
 
subplot(311);
plot(TOUT, (180/pi)*aircraft_state(:,10),col);hold on;
title('Angular Velocity v Time');   
ylabel('p [deg/s]')    

subplot(312);
plot(TOUT, (180/pi)*aircraft_state(:,11),col);hold on;
 ylabel('q [deg/s]')    
 
figs(4) = subplot(313);
plot(TOUT, (180/pi)*aircraft_state(:,12),col);hold on;
ylabel('r [deg/s]')    
xlabel('time [sec]');

%%%%%%%%%%%%%%%%%%%%%%%%
figure(5);
figs(5) = plot3(aircraft_state(:,1),aircraft_state(:,2),-aircraft_state(:,3),col);hold on;
title("3D Plane Trajectory")
xlabel("X [m]")
ylabel("Y [m]")
zlabel("Z [m]")


%%%%%%%%%%%%%%%%%%%%%%%%
if (~isempty(control_surfaces))
    figure(6);
    
    subplot(411);
    plot(TOUT, control_surfaces(:,1),col);hold on;
    title('Control Surfaces v Time');   
    ylabel('elevator [rad]')    

    subplot(412);
    plot(TOUT, control_surfaces(:,2),col);hold on;
    ylabel('aileron [rad]')      
 
    subplot(413);
    plot(TOUT, control_surfaces(:,3),col);hold on;
    ylabel('rudder [rad]')       
    
    figs(6) = subplot(414);
    plot(TOUT, control_surfaces(:,4),col);hold on;
    ylabel('throttle [frac]')     
    xlabel('time [sec]');

end

if ~isempty(background_wind_array)
    figure(7);
    
    V_B = [];
    wind_angles = [];
    for k = 1:length(TOUT)
        vel = aircraft_state(k, 7:9)';%TransformFromInertialToBody(a, aircraft_state(k,4:6));
        angles = WindAnglesFromVelocityBody(vel);
        V_B = [V_B, vel];
        wind_angles = [wind_angles, angles];
    end
% 
%     wind_angles = WindAnglesFromVelocityBody(V_B);
    
    figs(7) = subplot(311);
    plot(TOUT, wind_angles(1,:), col); hold on;
    title("Wind Angles v Time")
    ylabel("Airspeed [m/s]")

    subplot(312);
    plot(TOUT, wind_angles(2,:), col)
    ylabel("\beta [rad]")

    subplot(313);
    plot(TOUT, wind_angles(3,:), col)
    ylabel("\alpha [rad]")
    xlabel("time [sec]")

end

