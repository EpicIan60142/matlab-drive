%% Senior Projects Packet Reception Rate Plotting
% - Ian Faber, with modified code from Ryan Caputo

%% Housekeeping
clc; clear; close all

%% Setup
maxPackets_WiFi = 22; % Max number of WiFi packets transmitted
maxPackets_BT = 91; % Max number of Bluetooth packets transmitted (~4 BT packets per 1 WiFi packet)

%% Extract data - WiFi
fprintf("Processing WiFi Packets!\n\n")

% Below code heavily modified from Ryan Caputo's heat map algorithm

% Wednesday test
folder = "range_and_gain_test_data_april10"; % Note, angles are reversed (90 is 270 and vice versa)
addpath(folder)

file_names = dir(folder);
file_names = file_names(3:end);

idx = 21; % last wifi packet

count = 0;
summary = []; % Summary of each coordinate

for i = 1:idx % Only look at wifi for now
    data = readstruct(file_names(i).name);
    
    range = str2num(cell2mat(extractBetween(file_names(i).name, "_range_", "_altitude_")));
    altitude = str2num(cell2mat(extractBetween(file_names(i).name, "_altitude_","_angle_")));
    try
        angle = str2num(cell2mat(extractBetween(file_names(i).name, "_angle_","_")));
    catch
        angle = str2num(cell2mat(extractBetween(file_names(i).name, "_angle_","."))); % There was an "attempt_2"
    end

    if isempty(angle) % Check for test cases
        angle = 9999;
    end

    % loop through the number of RID packets per data frame
    for j = 1:length(data)        
        % extract latitude and longitude as "strings" (double quotations!)
        lat_s = data(j).x_source.layers.opendroneid.opendroneid_message_pack.opendroneid_message_location.OpenDroneID_loc_lat; % deg
        lon_s = data(j).x_source.layers.opendroneid.opendroneid_message_pack.opendroneid_message_location.OpenDroneID_loc_lon; % deg
        alt_s = data(j).x_source.layers.opendroneid.opendroneid_message_pack.opendroneid_message_location.OpenDroneID_loc_geoAlt; % ft

        %lat_s = data(j).rid.opendroneid_message_system.OpenDroneID_system_lat; % deg
        %lon_s = data(j).rid.opendroneid_message_system.OpenDroneID_system_lon; % deg

        if lat_s == '0'
            % ...then GPS has not locked
            fprintf('%s\n', 'GPS not locked, skipping RID packet...');
            continue
        end

        count = count + 1; % increment the count variable for indexing

        % convert from string to char
        lat_s = char(lat_s);
        lon_s = char(lon_s);
        alt_s = char(alt_s);

        % add a decimal point, and then convert to double
        lat_all(count) = str2double([lat_s(1:2), '.', lat_s(3:end)]);
        lon_all(count) = str2double([lon_s(1:4), '.', lon_s(5:end)]);
        % alt_all(count) = str2double(alt_s)/3.283; % convert to m
        alt_all(count) = altitude;

        strength(count) = str2double(data(j).x_source.layers.radiotap.radiotap_dbm_antsignal); % dBm (watts)

        rate(count) = length(data)/maxPackets_WiFi; % Calculate reception rate 

        % Tracker for diagnostics
        summarySlice = [range, altitude, angle, lat_all(count), lon_all(count), alt_all(count), rate(count)];
        summary = [summary; summarySlice];
    end
end

% Friday test
folder = "range_and_gain_test_data_april12"; % Angles back to normal
addpath(folder)

file_names = dir(folder);
file_names = file_names(3:end);

idx = 48; % last wifi packet
test2idx = count + 1;

for i = 1:idx % Only look at wifi for now
    data = readstruct(file_names(i).name);
    
    range = str2num(cell2mat(extractBetween(file_names(i).name, "_range_", "_altitude_")));
    altitude = str2num(cell2mat(extractBetween(file_names(i).name, "_altitude_","_angle_")));
    try
        angle = str2num(cell2mat(extractBetween(file_names(i).name, "_angle_","_")));
    catch
        angle = str2num(cell2mat(extractBetween(file_names(i).name, "_angle_","."))); % There was no "attempt_2"
    end

    if isempty(angle) % Check for test cases
        angle = 9999;
    end

    % loop through the number of RID packets per data frame
    for j = 1:length(data)        
        % extract latitude and longitude as "strings" (double quotations!)
        lat_s = data(j).x_source.layers.opendroneid.opendroneid_message_pack.opendroneid_message_location.OpenDroneID_loc_lat; % deg
        lon_s = data(j).x_source.layers.opendroneid.opendroneid_message_pack.opendroneid_message_location.OpenDroneID_loc_lon; % deg
        alt_s = data(j).x_source.layers.opendroneid.opendroneid_message_pack.opendroneid_message_location.OpenDroneID_loc_geoAlt; % ft

        %lat_s = data(j).rid.opendroneid_message_system.OpenDroneID_system_lat; % deg
        %lon_s = data(j).rid.opendroneid_message_system.OpenDroneID_system_lon; % deg

        if lat_s == '0'
            % ...then GPS has not locked
            fprintf('%s\n', 'GPS not locked, skipping RID packet...');
            continue
        end

        count = count + 1; % increment the count variable for indexing

        % convert from string to char
        lat_s = char(lat_s);
        lon_s = char(lon_s);
        alt_s = char(alt_s);

        % add a decimal point, and then convert to double
        lat_all(count) = str2double([lat_s(1:2), '.', lat_s(3:end)]);
        lon_all(count) = str2double([lon_s(1:4), '.', lon_s(5:end)]);
        % alt_all(count) = str2double(alt_s)/3.283; % Convert to m
        alt_all(count) = altitude;

        strength(count) = str2double(data(j).x_source.layers.radiotap.radiotap_dbm_antsignal); % dBm (watts)

        rate(count) = 100*(length(data)/maxPackets_WiFi); % Calculate reception rate

        % Tracker for diagnostics
        summarySlice = [range, altitude, angle, lat_all(count), lon_all(count), alt_all(count), rate(count)];
        summary = [summary; summarySlice];

    end
end

%% Find [x, y] Coordinates - WiFi

zRot = [3,1,1]; % Dummy 3-1-1 EA rotation for z-axis rotating
EAfix = EA2DCM(deg2rad([120,0,0]),zRot); % 120 degree clockwise rotation to fix orientation between Wednesday and Friday

x = zeros(length(lat_all), 1);
y = zeros(length(lat_all), 1);
z = zeros(length(lat_all), 1);
for i = 1:test2idx - 1
    [x(i), y(i), z(i)] = latlon_to_xyz(lat_all(i), lon_all(i), alt_all(i), lat_all(1), lon_all(1), alt_all(1));
    
    Xfix = EAfix*[x(i); y(i); z(i)]; % twist the node to match Friday

    if summary(i,3) ~= 9999
        DCM = EA2DCM(deg2rad([-summary(i,3),0,0]),zRot); % Rotate node, 90 and 270 are reversed
    else
        DCM = eye(3);
    end

    X = DCM*Xfix;

    x(i) = X(1);
    y(i) = X(2);
    z(i) = X(3);
end

for i = test2idx:length(lat_all) % Change node location on Friday
    [x(i), y(i), z(i)] = latlon_to_xyz(lat_all(i), lon_all(i), alt_all(i), lat_all(test2idx), lon_all(test2idx), alt_all(test2idx));

    if summary(i,3) ~= 9999
        DCM = EA2DCM(deg2rad([summary(i,3),0,0]),zRot); % Rotate node
    else
        DCM = eye(3);
    end

    X = [x(i); y(i); z(i)];

    X = DCM*X;

    x(i) = X(1);
    y(i) = X(2);
    z(i) = X(3);
end

%% Plot points in 3D space - WiFi
figure
hold on
grid on
title("WiFi RID Packet Reception Rate Test", 'FontSize', 16)
scatter3(x,y,z,20,rate,'filled')
colormap("jet")
cb = colorbar;
ylabel(cb, "Packet Reception Rate [%]", 'FontSize', 16, 'Rotation', 270)
xlabel("Range [m]",'FontSize',16)
ylabel("Crossrange [m]",'FontSize',16)
zlabel("Altitude [m]",'FontSize',16)

view([30, 35])

%% Extract data - Bluetooth
fprintf("\nProcessing Bluetooth Packets!\n\n")

clear lat_all lon_all alt_all strength rate % Clear variables of interest

% Below code heavily modified from Ryan Caputo's heat map algorithm

% Wednesday test
folder = "range_and_gain_test_data_april10"; % Note, angles are reversed (90 is 270 and vice versa)
addpath(folder)

file_names = dir(folder);
file_names = file_names(3:end);

idx = 21; % last wifi packet

count = 0;
summary = []; % Summary of each coordinate

for i = idx+1:length(file_names) % Look at bluetooth only
    data = readstruct(file_names(i).name);
    
    range = str2num(cell2mat(extractBetween(file_names(i).name, "_range_", "_altitude_")));
    altitude = str2num(cell2mat(extractBetween(file_names(i).name, "_altitude_","_angle_")));
    try
        angle = str2num(cell2mat(extractBetween(file_names(i).name, "_angle_","_")));
    catch
        angle = str2num(cell2mat(extractBetween(file_names(i).name, "_angle_","."))); % There was an "attempt_2"
    end

    if isempty(angle) % Check for test cases
        angle = 9999;
    end

    % loop through the number of RID packets per data frame
    for j = 1:length(data)
        try
            % extract latitude and longitude as "strings" (double quotations!)
            lat_s = data(j).x_source.layers.opendroneid.opendroneid_message_pack.opendroneid_message_location.OpenDroneID_loc_lat; % deg
            lon_s = data(j).x_source.layers.opendroneid.opendroneid_message_pack.opendroneid_message_location.OpenDroneID_loc_lon; % deg
            alt_s = data(j).x_source.layers.opendroneid.opendroneid_message_pack.opendroneid_message_location.OpenDroneID_loc_geoAlt; % ft
    
            %lat_s = data(j).rid.opendroneid_message_system.OpenDroneID_system_lat; % deg
            %lon_s = data(j).rid.opendroneid_message_system.OpenDroneID_system_lon; % deg
    
            if lat_s == '0'
                % ...then GPS has not locked
                fprintf('%s\n', 'GPS not locked, skipping RID packet...');
                continue
            end
    
            count = count + 1; % increment the count variable for indexing
    
            % convert from string to char
            lat_s = char(lat_s);
            lon_s = char(lon_s);
            alt_s = char(alt_s);
    
            % add a decimal point, and then convert to double
            lat_all(count) = str2double([lat_s(1:2), '.', lat_s(3:end)]);
            lon_all(count) = str2double([lon_s(1:4), '.', lon_s(5:end)]);
            % alt_all(count) = str2double(alt_s)/3.283; % convert to m
            alt_all(count) = altitude;
    
            strength(count) = str2double(data(j).x_source.layers.nordic_ble.nordic_ble_rssi); % dBm (watts)
    
            rate(count) = 100*(length(data)/maxPackets_BT); % Calculate reception rate

            % Tracker for diagnostics
            summarySlice = [range, altitude, angle, lat_all(count), lon_all(count), alt_all(count), rate(count)];
            summary = [summary; summarySlice];
        catch
            fprintf("No Message Pack, no data... \n")
        end
        
    end
end

% Friday test
folder = "range_and_gain_test_data_april12"; % Angles back to normal
addpath(folder)

file_names = dir(folder);
file_names = file_names(3:end);

idx = 48; % last wifi packet
test2idx = count + 1;

for i = idx+1:length(file_names) % Look at Bluetooth only
    data = readstruct(file_names(i).name);
    
    range = str2num(cell2mat(extractBetween(file_names(i).name, "_range_", "_altitude_")));
    altitude = str2num(cell2mat(extractBetween(file_names(i).name, "_altitude_","_angle_")));
    try
        angle = str2num(cell2mat(extractBetween(file_names(i).name, "_angle_","_")));
    catch
        angle = str2num(cell2mat(extractBetween(file_names(i).name, "_angle_","."))); % There was no "attempt_2"
    end

    if isempty(angle) % Check for test cases
        angle = 9999;
    end

    % loop through the number of RID packets per data frame
    for j = 1:length(data)
        try
            % extract latitude and longitude as "strings" (double quotations!)
            lat_s = data(j).x_source.layers.opendroneid.opendroneid_message_pack.opendroneid_message_location.OpenDroneID_loc_lat; % deg
            lon_s = data(j).x_source.layers.opendroneid.opendroneid_message_pack.opendroneid_message_location.OpenDroneID_loc_lon; % deg
            alt_s = data(j).x_source.layers.opendroneid.opendroneid_message_pack.opendroneid_message_location.OpenDroneID_loc_geoAlt; % ft
    
            %lat_s = data(j).rid.opendroneid_message_system.OpenDroneID_system_lat; % deg
            %lon_s = data(j).rid.opendroneid_message_system.OpenDroneID_system_lon; % deg
    
            if lat_s == '0'
                % ...then GPS has not locked
                fprintf('%s\n', 'GPS not locked, skipping RID packet...');
                continue
            end
    
            count = count + 1; % increment the count variable for indexing
    
            % convert from string to char
            lat_s = char(lat_s);
            lon_s = char(lon_s);
            alt_s = char(alt_s);
    
            % add a decimal point, and then convert to double
            lat_all(count) = str2double([lat_s(1:2), '.', lat_s(3:end)]);
            lon_all(count) = str2double([lon_s(1:4), '.', lon_s(5:end)]);
            % alt_all(count) = str2double(alt_s)/3.283; % Convert to m
            alt_all(count) = altitude;
    
            strength(count) = str2double(data(j).x_source.layers.nordic_ble.nordic_ble_rssi); % dBm (watts)
    
            rate(count) = length(data)/maxPackets_BT; % Calculate reception rate

            % Tracker for diagnostics
            summarySlice = [range, altitude, angle, lat_all(count), lon_all(count), alt_all(count), rate(count)];
            summary = [summary; summarySlice];
        catch
            fprintf("No Message Pack, no data... \n")
        end
    end
end

%% Find [x, y] Coordinates - Bluetooth

zRot = [3,1,1]; % Dummy 3-1-1 EA rotation for z-axis rotating
EAfix = EA2DCM(deg2rad([120,0,0]),zRot); % 120 degree clockwise rotation to fix orientation between Wednesday and Friday

x = zeros(length(lat_all), 1);
y = zeros(length(lat_all), 1);
z = zeros(length(lat_all), 1);
for i = 1:test2idx - 1
    [x(i), y(i), z(i)] = latlon_to_xyz(lat_all(i), lon_all(i), alt_all(i), lat_all(1), lon_all(1), alt_all(1));
    
    Xfix = EAfix*[x(i); y(i); z(i)]; % twist the node to match Friday

    if summary(i,3) ~= 9999
        DCM = EA2DCM(deg2rad([-summary(i,3),0,0]),zRot); % Rotate node, 90 and 270 are reversed
    else
        DCM = eye(3);
    end

    X = DCM*Xfix;

    x(i) = X(1);
    y(i) = X(2);
    z(i) = X(3);
end

for i = test2idx:length(lat_all) % Change node location on Friday
    [x(i), y(i), z(i)] = latlon_to_xyz(lat_all(i), lon_all(i), alt_all(i), lat_all(test2idx), lon_all(test2idx), alt_all(test2idx));

    if summary(i,3) ~= 9999
        DCM = EA2DCM(deg2rad([summary(i,3),0,0]),zRot); % Rotate node
    else
        DCM = eye(3);
    end

    X = [x(i); y(i); z(i)];

    X = DCM*X;

    x(i) = X(1);
    y(i) = X(2);
    z(i) = X(3);
end

%% Plot points in 3D space - Bluetooth
figure
hold on
grid on
title("Bluetooth RID Packet Reception Test", 'FontSize', 16)
scatter3(x,y,z,20,rate,'filled')
colormap("jet")
cb = colorbar;
ylabel(cb, "Packet Reception Rate [%]", 'FontSize', 16, 'Rotation', 270)
xlabel("Range [m]",'FontSize',16)
ylabel("Crossrange [m]",'FontSize',16)
zlabel("Altitude [m]",'FontSize',16)
xlim([-150,150]);
ylim([-150,150]);

view([30, 35])


%% Function to go from Lat/Lon/Alt to XYZ Coordinates (Relative to Node Reference Point)

function [x, y, z, r, theta, r2] = latlon_to_xyz(lat, lon, alt, reflat, reflon, refalt)
    R_earth = 6371e3; % radius of earth in meter
    % convert to radians
    lat = lat * pi/180;
    lon = lon * pi/180;
    reflat = reflat * pi/180;
    reflon = reflon * pi/180;

    r = R_earth * acos(sin(lat) * sin(reflat) + cos(lat) * cos(reflat) * cos(lon - reflon));

    x = (lon - reflon) * cos((lat + reflat) / 2);
    y = lat - reflat;
    r2 = R_earth * sqrt(x^2 + y^2);

    theta = atan2(sin(lon-reflon) * cos(lat-reflat), cos(reflat) * sin(lat) - sin(reflat) * cos(lat) * cos(lon - reflon));

    %x = r * cos(theta);
    %y = r * sin(theta);
    x = r * sin(theta); % north is defined as zero degrees!
    y = r * cos(theta); % north is defined as zero degrees!
    z = alt - refalt;
end





