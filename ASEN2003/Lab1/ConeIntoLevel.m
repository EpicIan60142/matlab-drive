clear;
clc;
close all;

s = (2*pi) - atan(1/3):pi/200:2*pi;

x = 0*s;
y = (40).*sin(s + pi);
z = (40).*cos(s + pi) + 40;
coords = [x; y; z];

figure
grid on;
hold on;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
xlim([-45,45]);
ylim([-45,45]);
zlim([-45,45]);
view([30 35]);

% StartHeight = 77.8761;
% GRadials = [length(coords)];
% Gindex = 1;
% 
% for i = 1:length(coords)
%     Radius = 40;
%     
%     GRadials(Gindex) = (2 * 9.81 * (125 - (StartHeight + z(Gindex))))/ (Radius * 9.81);
%     Gindex = Gindex + 1;
% end

plot = plot3(coords(1,:), coords(2,:), coords(3,:), 'LineWidth', 5);