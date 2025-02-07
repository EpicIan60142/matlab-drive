% Use plot3

clc; clear; close all;

s = 0:pi/200:4*pi;

x = s;
y = zeros(1,length(s));
z = cos(2*(s-(acos(-1/4)+2*pi)/2))+cos(4*(s-(acos(-1/4)+2*pi)/2));

coords = [x; y; z];
coords2 = rotate(coords, 0, 90, 0);

figure
hold on;

xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
view([30 35]);

plot3(coords(1,:), coords(2,:), coords(3,:), 'LineWidth',5);

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