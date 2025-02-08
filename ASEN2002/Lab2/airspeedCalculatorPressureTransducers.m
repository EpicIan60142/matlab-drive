clear; clc; close all;

%% file loading
numFiles = 7;

for i = 1:numFiles %find name of each group data file
    % grab the data name and put it in a string
    files{i} = dir(['ASEN 2002_1030_',num2str(i),'*.csv']);
    % maybe name for each section like "files1030"
end

%for j=1:5
%    x=files{j};
%    num_elements(j) = sum(arrayfun(@(x) ~isempty(x.name),x)); %number of datafiles for each port
%end

clear i j;
R = 287;
A_1 = 9.5;
A_2 = 1;

for i=1:numFiles
        % load data using the data name into an array
        data = load(files{i}.name);    % need to call info{1,#} to get data out
        
        % find index where voltages change and use them to make smaller sections
        % of data
        % create a vector of voltages
        voltages = data(:,13);
        
        % create a vector of indexes where voltage changes
        v_Change_Index = find(diff(voltages) > 1);
        
        % add start and end to the indexes -- these are the different groups of
        % voltages
        v_Change_Index = [0;v_Change_Index;length(voltages)];
        
        % preallocate arraya for speed
        sections = zeros(length(v_Change_Index)-1,13);
        
        % loop through and put data into multiple sections
        for k = 1:length(v_Change_Index)-1
            % collapse the data into a row for each section by taking mean of the
            % section
            sections(k,:)= mean(data((v_Change_Index(k)+1:v_Change_Index(k+1)),:));
        end
        
        %% need to put measurements in structs
        
        % Column 1: Atmospheric temperature [K]
        % Column 2: Atmospheric pressure [Pa]
        % Column 3: Atmospheric density (kg/m^3)
        % Column 4: Airspeed (m/s)
        % Column 5: Transducer 1 Differential pressure [Pa]
        % Column 6: Transducer 2 Differential pressure [Pa]
        % Column 7: Angle of Attack [deg]
        % Column 8: Sting Normal Force [N]
        % Column 9: Sting Axial Force [N]
        % Column 10: Sting Pitching moment [N*m]
        % Column 11: ELD Probe X-axis [mm]
        % Column 12: ELD Probe Y-axis [mm]
        % Column 13: Applied voltage [V]
        
        T_atm = sections(:,1);
        P_atm = sections(:,2);
        rho_atm = sections(:,3);
        airspeedTest = sections(:,4);
        venturiPressure = sections(:,5);
        pitotPressure = sections(:,6);
        angleAttack = sections(:,7);
        normalForce = sections(:,8);
        axialForce = sections(:,9);
        pitchingMoment = sections(:,10);
        ELDX = sections(:,11);
        ELDY = sections(:,12);
        voltage = sections(:,13);
        
        measurements(i).T_atm = T_atm;
        measurements(i).P_atm = P_atm;
        measurements(i).venturiPressure = venturiPressure;
        measurements(i).pitotPressure = pitotPressure;
        measurements(i).voltage = voltage;
        measurements(i).airSpeedVenturi = sqrt(2*venturiPressure.*((R*T_atm)./P_atm));
        measurements(i).airSpeedPitot = sqrt((2*pitotPressure*R.*T_atm)./(P_atm*(1-(A_2/A_1)^2)));
        
    
end

%% Pressure transducers
measurements;
measurements.T_atm;
measurements.P_atm;
measurements.venturiPressure;
measurements.pitotPressure;
measurements.airSpeedVenturi;
measurements.airSpeedPitot;
measurements.voltage;

allMeasurements = struct2table(measurements)

final.voltage = [];
final.airSpeedVenturi = [];
final.airSpeedPitot = [];

clear i j;

for i = 1:length(allMeasurements.voltage)
    final.voltage = [final.voltage; allMeasurements.voltage{i}];
    final.airSpeedVenturi = [final.airSpeedVenturi; allMeasurements.airSpeedVenturi{i}];
    final.airSpeedPitot = [final.airSpeedPitot; allMeasurements.airSpeedPitot{i}];
end

clear i j;

final.voltage;
final.airSpeedVenturi;
final.airSpeedPitot;

for i = 1:length(final.voltage)
    for j = (i+1):length(final.voltage)
        if(final.voltage(i) > final.voltage(j))
            
            temp = final.voltage(i);
            final.voltage(i) = final.voltage(j);
            final.voltage(j) = temp;
            
            temp = final.airSpeedVenturi(i);
            final.airSpeedVenturi(i) = final.airSpeedVenturi(j);
            final.airSpeedVenturi(j) = temp;
            
            temp = final.airSpeedPitot(i);
            final.airSpeedPitot(i) = final.airSpeedPitot(j);
            final.airSpeedPitot(j) = temp;
            
        end
    end
end

final.voltage;
final.airSpeedVenturi;
final.airSpeedPitot;

A = [final.voltage, ones(length(final.voltage),1)];
d = final.airSpeedVenturi;
x_hat_1 = (A'*A)^-1*A'*d
LeastSquares1 = @(x)x_hat_1(1)*x + x_hat_1(2);

B = [final.voltage, ones(length(final.voltage),1)];
e = final.airSpeedPitot;
x_hat_2 = (B'*B)^-1*B'*e
LeastSquares2 = @(x)x_hat_2(1)*x + x_hat_2(2);

figure

subplot(2,1,1)
hold on
title("Voltage vs. Venturi Tube-Pressure Tranducer Airspeed")
plot(final.voltage, final.airSpeedVenturi, 'b.')
plot(final.voltage, LeastSquares1(final.voltage),'r')
xlabel("Voltage (V)")
ylabel("Venturi Airspeed (m/s)")
regLabel1 = sprintf("v = %.3f*V + %.3f", x_hat_1(1), x_hat_1(2));
legend("Raw data", regLabel1,'Location','best')
hold off

subplot(2,1,2)
hold on
title("Voltage vs. Pitot-Static Probe-Pressure Transducer Airspeed")
plot(final.voltage, final.airSpeedPitot, 'b.')
plot(final.voltage, LeastSquares2(final.voltage),'r')
xlabel("Voltage (V)")
ylabel("Airspeed (m/s)")
regLabel2 = sprintf("v = %.3f*V + %.3f", x_hat_2(1), x_hat_2(2));
legend("Raw data", regLabel2,'Location','best')
hold off

