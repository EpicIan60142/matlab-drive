%% Starting Parameters
h0 = 125;
pos_in = [0, 0, h0];
n = 100; % number of elements for each track feature

%%%%%%%%%%%%%%%%%%% PATH CREATION %%%%%%%%%%%%%%%%%%%%
% Call each track element in succession. Output velocity becomes input velocity of
% next element. Last point in output path becomes pos_in for next element.
% Append path, distance, and g_s outputs to variables of the same name for full coaster. 
%% Element 1: Downward Ramp %%
d = 30;
theta = -pi/4;
[path, distance, vel_out, g_s] = ramp(pos_in, h0, d, theta, n);

%% Transition 1 %%
%transition params
vel_in = vel_out;
pos_in = path(end,:);
vel_exit_unitvec = [1 0 0];
r = 10;
% call transition here

%% Element 2: Banked Turn %%
% turn params
bank_angle = pi/4;
turn_dir = 1;
vel_in = vel_out;
pos_in = path(end,:);
r = 10;
% call banked turn here

%% Transistion 2 %%
vel_in = vel_out;
vel_in(3) = round(vel_in(3)); %transition angles can be sensitive to small errors in vel so we round
pos_in = path(end,:);
vel_exit_unitvec = [-1 0 1];
r = 10;
% call transition here

%% Element 3: Zero-G Parabola %%
vel_in = vel_out;
pos_in = path(end,:);
a = -9.81; % zero-g
d = 50;
% call parabola here

%% Transition 3 %%
vel_in = vel_out;
pos_in = path(end,:);
vel_exit_unitvec = [-1 0 0];
r = 20;
% call transition here

%% Element 4: Banked Turn %%
bank_angle = pi/4;
turn_dir = 1;
vel_in = vel_out;
pos_in = path(end,:);
r = 40;
% call banked turn here

%% Element 5: Circular Loop %%
prop = 1;
r = 17;
vel_in = vel_out;
pos_in = path(end,:);
% call circular loop here

%% Element 6: Braking Section %%
vel_in = vel_out;
pos_in = path(end,:);
d = 50;
% call braking section here

%%%%%%%%%%%%%%%% PATH VISUALIZATION %%%%%%%%%%%%%%%
%% 3D plot track
figure
hold on
plot3(path(:,1),path(:,2),path(:,3),'Color','[0.8500 0.3250 0.0980]', 'linewidth', 2)
plot3(path(1,1),path(1,2),path(1,3), 'go', 'linewidth', 2);
plot3(path(end,1),path(end,2),path(end,3), 'ro', 'linewidth', 2);
grid on
title('Path of Rollercoaster')
xlabel('x-direction(m)')
ylabel('y-direction(m)')
zlabel('z-direction(m)')
legend('Path','Start','Finish')
axis equal
view([-.4 -1 .45])

%% Plot Gs(vertical, lateral and horizontal)
% horizontal
figure
hold on
grid on
plot(distance, g_s(:,1),'Color','[0.8500 0.3250 0.0980]', 'linewidth', 1.5)
plot(distance, ones(1,length(distance))*5, 'r--', 'linewidth', 1.5)
plot(distance, -ones(1,length(distance))*4, 'r--', 'linewidth', 1.5)
legend('Gs','Allowable Gs Bound')
title('Horizontal Gs vs Distance Traveled along Track')
xlabel('Distance(m)')
ylabel('Gs')
xlim([0 distance(end)])
ylim([-5 6])

% lateral 
figure
hold on
grid on
plot(distance, g_s(:,2),'Color','[0.8500 0.3250 0.0980]', 'linewidth', 1.5)
plot(distance, ones(1,length(distance))*3, 'r--', 'linewidth', 1.5)
plot(distance, ones(1,length(distance))*(-3), 'r--', 'linewidth', 1.5)
legend('Gs','Allowable Gs Bound')
title('Lateral Gs vs Distance Traveled along Track')
xlabel('Distance(m)')
ylabel('Gs')
xlim([0 distance(end)])
ylim([-4 4])

% vertical
figure
hold on
grid on
plot(distance, g_s(:,3),'Color','[0.8500 0.3250 0.0980]', 'linewidth', 1.5)
plot(distance, ones(1,length(distance))*6, 'r--', 'linewidth', 1.5)
plot(distance, ones(1,length(distance))*(-1), 'r--', 'linewidth', 1.5)
legend('Gs','Allowable Gs Bound')
title('Vertical Gs vs Distance Traveled along Track')
xlabel('Distance(m)')
ylabel('Gs')
xlim([0 distance(end)])
ylim([-2 7])

%% FUNCTIONS - Paste your functions in here

%%%%%LOOP%%%%%

%%%%%%%%%%%%%%

%%%%PARABOLA%% 

%%%%%%%%%%%%%%

%%TRANSITION%%

%%%%%%%%%%%%%%

%%BANKED TURN%

%%%%%%%%%%%%%%

%%%%BRAKING%%%%

%%%%%%%%%%%%%%