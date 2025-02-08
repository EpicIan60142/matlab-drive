%% Starting Parameters
h0 = 125;
pos_in = [0, 0, h0];
n = 100; % number of elements for each track feature

%%%%%%%%%%%%%%%%%%% PATH CREATION %%%%%%%%%%%%%%%%%%%%
%% Element 1: Downward Ramp - IMPLEMENTED FOR YOU %%
d = 30;
theta = -pi/4;
[path, distance, vel_out, g_s] = ramp(pos_in, h0, d, theta, n);

%% Transition 1 - IMPLEMENTED FOR YOU %%
vel_in = vel_out;
pos_in = path(end,:);
vel_exit_unitvec = [1 0 0];
r = 10;
[path(end+1:end+n,:), distance_el, vel_out, g_s(end+1:end+n,:)] = ...
    transition(vel_in, pos_in, vel_exit_unitvec, r, 1,h0, n);
distance(end+1:end+n) = distance_el + distance(end);

%% Element 2: Banked Turn - YOU IMPLEMENT %%
[] = banked_turn(vel_in, pos_in, n, r, turn_dir, bank_angle, h0, NaN, NaN, NaN);

%% Transistion 2 - IMPLEMENTED FOR YOU %%
vel_in = vel_out;
pos_in = path(end,:);
vel_exit_unitvec = [-1 0 1];
r = 10;
[path(end+1:end+n,:), distance_el, vel_out, g_s(end+1:end+n,:)] = ...
    transition(vel_in, pos_in, vel_exit_unitvec, r, 1, h0, n);
distance(end+1:end+n) = distance_el + distance(end);

%% Element 3: Zero-G Parabola - YOU IMPLEMENT%%
[] = parabola(vel_in, pos_in, n, a, d, h0, NaN,NaN,NaN);

%% Transition 3 - IMPLEMENTED FOR YOU %%
vel_in = vel_out;
pos_in = path(end,:);
vel_exit_unitvec = [-1 0 0];
r = 20;
[path(end+1:end+n,:), distance_el, vel_out, g_s(end+1:end+n,:)] = ...
    transition(vel_in, pos_in, vel_exit_unitvec, r, -1, h0, n);
distance(end+1:end+n) = distance_el + distance(end);

%% Element 4: Banked Turn - YOU IMPLEMENT %%
[] = banked_turn(vel_in, pos_in, n, r, turn_dir, bank_angle, h0, NaN, NaN, NaN);

%% Element 5: Circular Loop - YOU IMPLEMENT%%
[] = loop(vel_in, pos_in, n, prop, r, h0, NaN,NaN,NaN);

%% Element 6: Braking Section - YOU IMPLEMENT%%
[] = braking(vel_in, pos_in, n, d, h0, NaN, NaN);

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
