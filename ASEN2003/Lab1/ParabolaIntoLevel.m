clear
clc

clear 
clc

s = pi+atan(1/32):pi/200:3*pi/2;

x = 0*s;
y = (40).*cos(s);
z = (40).*sin(s) + 40;
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

% StartHeight = 77.8761 - 2.1;
% GRadials = [length(coords)];
% Gindex = 1;
% 
% for i = 1:length(coords)
%     Radius = 40;
%     
%     GRadials(Gindex) = (2 * 9.81 * (125 - (StartHeight + z(Gindex))))/ (Radius * 9.81); % The 37.9473319220205 account for an offset due to the used graphing method
%     
%     Gindex = Gindex + 1;
% end

plot = plot3(coords(1,:), coords(2,:), coords(3,:), 'LineWidth', 5);