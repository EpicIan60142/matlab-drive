% Use plot3

clc; clear; close all;

g = 9.81;
h0 = 125;

s = 0:pi/200:5*pi;

x = ((2*s+10)).*cos(s);
y = ((2*s+10)).*sin(s);
z = (-3*s);

coords = [x; y; z];
coords2 = rotate(coords, 0, 0, 0);

figure
hold on;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
view([30 35]);

%plot3(coords(1,:), coords(2,:), coords(3,:), 'LineWidth',5);
plot = plot3(coords2(1,:), coords2(2,:), coords2(3,:), 'LineWidth', 5);

StartHeight = input('Cone high point')
GForces = [length(coords)];
Gindex = 1;

for i = 1:length(coords)
    Radius = sqrt(((x(Gindex)^2) + (y(Gindex)^2)));
    
    GForces(Gindex) = (2 * 9.81 * (125 - (StartHeight + z(Gindex))))/ (Radius * 9.81);
    
    Gindex = Gindex + 1;
end




function newCoords = rotate(oldCoords, yaw, pitch, roll)
    % Yaw rotates about z axis
    % Pitch rotates about y axis, 
    % Roll rotates about x axis
    % Angles taken in degrees
    
    pitch = pi*(pitch/180);
    yaw = pi*(yaw/180);
    roll = pi*(roll/180);

    rotateMatrix = [
                    cos(yaw)*cos(pitch), cos(yaw)*sin(pitch)*sin(roll) - sin(yaw)*cos(roll), cos(yaw)*sin(pitch)*cos(roll) + sin(yaw)*sin(roll);
                    sin(yaw)*cos(pitch), sin(yaw)*sin(pitch)*sin(roll) + cos(yaw)*cos(roll), sin(yaw)*sin(pitch)*cos(roll) - cos(yaw)*sin(roll);
                    -sin(pitch)        , cos(pitch)*sin(roll)                              , cos(pitch)*cos(roll)
                   ];
    
    newCoords = rotateMatrix*oldCoords;
end
