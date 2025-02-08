% Use plot3

clc; clear; close all;

s = 0:pi/200:2*pi;

x = (5*cos(s));
y = (5*sin(s));
z = (s/3);

coords = [x; y; z];
coords2 = rotate(coords, 0, 90, 0);

figure
hold on;
xlim([-20, 20]);
ylim([-20, 20]);
zlim([-20, 20]);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
view([30 35]);

plot3(coords(1,:), coords(2,:), coords(3,:), 'LineWidth',5);
plot3(coords2(1,:) + 10, coords2(2,:) + 3, coords2(3,:) + 5, 'LineWidth', 5);

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
