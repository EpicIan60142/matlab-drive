%% ASEN 2002 Lab 2
% Wind Tunnel Airspeed Calculator
% Code by Group 7, Subteam A: Ian Faber, Andrew Kabos, Tsuening Lee

%% Housekeeping and setup
clear; clc; close all;

clear i j;
R = 287; % Gas constant for air, J/(kg*K)
A_1 = 12*12*0.0254^2; % Inlet area, 12x12 in -> m^2
A_2 = A_1 / 9.5; % Test section area, A_2/A_1 = 9.5, m^2
SG = 0.826; % Specific gravity of manometer fluid
sigma_PT = 68.95; % 0.01(1-0)*6894.76, 1% of max range, psi -> Pa
sigma_Patm = 3450; % 0.015*(230-20), 1.5% of max range, Pa
sigma_Man = 10.28; % (1000*SG)*9.8*0.05*0.0254, +/- 0.1 in gradations, Pa


%% Pressure Transducer file loading
numFiles = 7; % Load 7 files
for i = 1:numFiles % Find the name of each group data file
   % Grab the name and put it in a string
   files{i} = dir(['ASEN 2002_1030_',num2str(i),'*.csv']);
end

files

%% Pressure Transducer file processing
for i=1:numFiles % Loop once over each file
       % Load data using the data name into an array
       data = load(files{i}.name);
      
       % Find index where voltages change and use them to make smaller sections
       % of data
       
       % Create a vector of voltages
       voltages = data(:,13);
      
       % Create a vector of indexes where voltage changes signigicantly
       v_Change_Index = find(diff(voltages) > 1);
      
       % Add start and end to the indexes -- these are the different groups of
       % voltages
       v_Change_Index = [0;v_Change_Index;length(voltages)];
      
       % Preallocate arrays for speed
       sections = zeros(length(v_Change_Index)-1,13);
      
       % Loop through and put data into multiple sections
       for k = 1:length(v_Change_Index)-1
           % Collapse the data into a row for each section by taking mean of the
           % section
           sections(k,:)= mean(data((v_Change_Index(k)+1:v_Change_Index(k+1)),:));
       end
      
       % Put measurements in structs
      
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
       
       % Calculate airspeeds using the appropriate equation
       measurements(i).airSpeedVenturi = sqrt(2*venturiPressure.*((R*T_atm)./P_atm));
       measurements(i).airSpeedPitot = sqrt((2*pitotPressure*R.*T_atm)./(P_atm*(1-(A_2/A_1)^2)));
      
  
end
 
%% Pressure Transducer processing
% Unsuppress for debugging
measurements;
measurements.T_atm;
measurements.P_atm;
measurements.venturiPressure;
measurements.pitotPressure;
measurements.airSpeedVenturi;
measurements.airSpeedPitot;
measurements.voltage;

% Convert measurements struct into a table
allMeasurements = struct2table(measurements);
 
% Preallocate final struct for pressure transducers
final.T_atm = [];
final.P_atm = [];
final.voltage = [];
final.airSpeedVenturi = [];
final.airSpeedPitot = [];

% Populate the final struct with information from the measurements struct
for i = 1:length(allMeasurements.voltage)
   final.T_atm = [final.T_atm; allMeasurements.T_atm{i}];
   final.P_atm = [final.P_atm; allMeasurements.P_atm{i}];
   final.voltage = [final.voltage; allMeasurements.voltage{i}];
   final.airSpeedVenturi = [final.airSpeedVenturi; allMeasurements.airSpeedVenturi{i}];
   final.airSpeedPitot = [final.airSpeedPitot; allMeasurements.airSpeedPitot{i}];
end

% Sort the final struct in ascending order

% Unsuppress for debugging
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

% Unsuppress for debugging
final;

% Perform a linear Least Squares regression on the voltage and airspeed 
% information in the final struct

A = [final.voltage, ones(length(final.voltage),1)];
d = final.airSpeedVenturi;
x_hat_1 = (A'*A)^-1*A'*d
LeastSquares1 = @(x)x_hat_1(1)*x + x_hat_1(2);
sanityVenturi = A\d;
sigmaVVenturiPT = Venturi_Tube_Uncertainty_Function(final.airSpeedVenturi, final.T_atm, final.P_atm, R, A_1, A_2, sigma_PT, sigma_Patm);

B = [final.voltage, ones(length(final.voltage),1)];
e = final.airSpeedPitot;
x_hat_2 = (B'*B)^-1*B'*e
LeastSquares2 = @(x)x_hat_2(1)*x + x_hat_2(2);
sanityPitot = B\e;
sigmaVPitotPT = Pitot_Static_Uncertainty_Function(final.airSpeedPitot, final.T_atm, final.P_atm, R, sigma_PT, sigma_Patm);
 
%% Manometers
Table= readtable('water_manometer_readings.csv');
venturiTubeMan= [];
pitotStaticProbeMan = [];
for i = 1:length(Table.WasTheManometerConnectedToTheVenturiTubeOrPitotStaticProbe_)
    if Table.WasTheManometerConnectedToTheVenturiTubeOrPitotStaticProbe_(i) == "Venturi tube"    
       %T(i,(4:end))
        venturiTubeMan = [venturiTubeMan ; Table(i,(4:end))];
    
%     end
% end
% 
% for  i = 1:length(Table.WasTheManometerConnectedToTheVenturiTubeOrPitotStaticProbe_{i})
%     if Table.WasTheManometerConnectedToTheVenturiTubeOrPitotStaticProbe_(i) == "Pitot Static Probe"    
    else
        pitotStaticProbeMan = [pitotStaticProbeMan ; Table(i,(4:end))];
    end
end
measureMan.venturiVoltage = [venturiTubeMan.Voltage_1; venturiTubeMan.Voltage_2;venturiTubeMan.Voltage_3;venturiTubeMan.Voltage_4;venturiTubeMan.Voltage_5];
measureMan.venturiHeights = [venturiTubeMan.ReadingForVoltage_1_heightChange_inH2O__; venturiTubeMan.ReadingForVoltage_2_heightChange_inH2O__;venturiTubeMan.ReadingForVoltage_3_heightChange_inH2O__;venturiTubeMan.ReadingForVoltage_4_heightChange_inH2O__;venturiTubeMan.ReadingForVoltage_5_heightChange_inH2O__];
 
 
measureMan.pitotVoltage = [pitotStaticProbeMan.Voltage_1; pitotStaticProbeMan.Voltage_2;pitotStaticProbeMan.Voltage_3;pitotStaticProbeMan.Voltage_4;pitotStaticProbeMan.Voltage_5];
measureMan.pitotHeights = [pitotStaticProbeMan.ReadingForVoltage_1_heightChange_inH2O__; pitotStaticProbeMan.ReadingForVoltage_2_heightChange_inH2O__;pitotStaticProbeMan.ReadingForVoltage_3_heightChange_inH2O__;pitotStaticProbeMan.ReadingForVoltage_4_heightChange_inH2O__;pitotStaticProbeMan.ReadingForVoltage_5_heightChange_inH2O__];
 
g = 9.81;
 
measureMan.venturiDiffpress = 1000*SG*measureMan.venturiHeights * 0.0254 * g;
measureMan.pitotDiffpress = 1000*SG*measureMan.pitotHeights * 0.0254 * g;
T_atmAvg = mean(T_atm);
P_atmAvg = mean(P_atm);
 
measureMan.venturiAirspeed = sqrt(2*measureMan.venturiDiffpress*((R*T_atmAvg)/P_atmAvg) );
 
measureMan.pitotAirspeed = sqrt((2*measureMan.pitotDiffpress*R*T_atmAvg)./(P_atmAvg*(1-(A_2/A_1)^2)));       
 
Mat = [measureMan.venturiVoltage measureMan.venturiAirspeed];
[UA,~,idx] = unique(Mat(:,1));
newMat3 = [UA,accumarray(idx,Mat(:,2),[],@mean)];
 
measureMan.venturiVoltage = newMat3(:,1);
measureMan.venturiAirspeed = newMat3(:,2);
 
Mat = [measureMan.pitotVoltage measureMan.pitotAirspeed];
[UA,~,idx] = unique(Mat(:,1));
newMat4 = [UA,accumarray(idx,Mat(:,2),[],@mean)];
 
measureMan.pitotVoltage = newMat4(:,1);
measureMan.pitotAirspeed = newMat4(:,2);
 
A1 = [measureMan.venturiVoltage, ones(length(measureMan.venturiVoltage),1)];
d1 = measureMan.venturiAirspeed;
xhatMan1 = A1\d1;
LeastSquaresMan1 = @(x) xhatMan1(1)*x+xhatMan1(2) ;
sigmaVVenturiMan = Venturi_Tube_Uncertainty_Function(measureMan.venturiAirspeed, T_atmAvg, P_atmAvg, R, A_1, A_2, sigma_Man, sigma_Patm);
 
A2 = [measureMan.pitotVoltage, ones(length(measureMan.pitotVoltage),1)];
d2 = measureMan.pitotAirspeed;
xhatMan2 = A2\d2;
LeastSquaresMan2 = @(x) xhatMan2(1)*x+xhatMan2(2);
sigmaVPitotMan = Pitot_Static_Uncertainty_Function(measureMan.pitotAirspeed, T_atmAvg, P_atmAvg, R, sigma_Man, sigma_Patm);
 
x_hat_avg = (x_hat_1 + x_hat_2 + xhatMan1 + xhatMan2)/4;
LeastSquaresAvg = @(x) x_hat_avg(1)*x + x_hat_avg(2);

%% Plots

% 4 plots, with uncertainty
figure
subplot(2,2,1)
hold on
title("Voltage vs. Venturi Tube - Pressure Transducer Airspeed")
plot(final.voltage, final.airSpeedVenturi, 'b.')
plot(final.voltage, LeastSquares1(final.voltage),'r')
errorbar(final.voltage,LeastSquares1(final.voltage),sigmaVVenturiPT,'LineStyle','none');
xlabel("Voltage (V)")
ylabel("Venturi Airspeed (m/s)")
regLabel1 = sprintf("v = %.3f*V + %.3f", x_hat_1(1), x_hat_1(2));
legend("Raw data", regLabel1,'Location','best')
hold off
 
subplot(2,2,2)
hold on
title("Voltage vs. Pitot-Static Probe - Pressure Transducer Airspeed")
plot(final.voltage, final.airSpeedPitot, 'b.')
plot(final.voltage, LeastSquares2(final.voltage),'r')
errorbar(final.voltage,LeastSquares2(final.voltage),sigmaVPitotPT,'LineStyle','none');
xlabel("Voltage (V)")
ylabel("Airspeed (m/s)")
regLabel2 = sprintf("v = %.3f*V + %.3f", x_hat_2(1), x_hat_2(2));
legend("Raw data", regLabel2,'Location','best')
hold off
 
subplot(2,2,3)
hold on
title("Voltage vs. Venturi Tube - Water Manometer Airspeed")
plot(measureMan.venturiVoltage, measureMan.venturiAirspeed, 'b.')
plot(measureMan.venturiVoltage, LeastSquaresMan1(measureMan.venturiVoltage),'r')
errorbar(measureMan.venturiVoltage,LeastSquaresMan1(measureMan.venturiVoltage),sigmaVVenturiMan,'LineStyle','none');
xlabel("Voltage (V)")
ylabel("Venturi Airspeed (m/s)")
regLabel3 = sprintf("v = %.3f*V + %.3f", xhatMan1(1), xhatMan1(2));
legend("Raw data", regLabel3,'Location','best')
hold off
 
subplot(2,2,4)
hold on
title("Voltage vs. Pitot-Static Probe - Water Manometer Airspeed")
plot(measureMan.pitotVoltage, measureMan.pitotAirspeed, 'b.')
plot(measureMan.pitotVoltage, LeastSquaresMan2(measureMan.pitotVoltage),'r')
errorbar(measureMan.pitotVoltage,LeastSquaresMan2(measureMan.pitotVoltage),sigmaVPitotMan,'LineStyle','none');
xlabel("Voltage (V)")
ylabel("Airspeed (m/s)")
regLabel4 = sprintf("v = %.3f*V + %.3f", xhatMan2(1), xhatMan2(2));
legend("Raw data", regLabel4,'Location','best')
hold off

% 2 plots, comparing uncertainty
figure
subplot(2,1,1)
hold on

title("Voltage vs. Venturi Tube Airspeeds")
rawPT1 = plot(final.voltage, final.airSpeedVenturi, 'b.');
linePT1 = plot(final.voltage, LeastSquares1(final.voltage),'r');
uncertaintyPT1 = errorbar(final.voltage,LeastSquares1(final.voltage),sigmaVVenturiPT,'LineStyle','none');

rawWM1 = plot(measureMan.venturiVoltage, measureMan.venturiAirspeed, 'm.');
lineWM1 = plot(measureMan.venturiVoltage, LeastSquaresMan1(measureMan.venturiVoltage),'g');
uncertaintyWM1 = errorbar(measureMan.venturiVoltage,LeastSquaresMan1(measureMan.venturiVoltage),sigmaVVenturiMan,'LineStyle','none');

subset1 = [linePT1, uncertaintyPT1, lineWM1, uncertaintyWM1];

xlabel("Voltage (V)")
ylabel("Venturi Airspeed (m/s)")
legend(subset1, regLabel1, "PT Uncertainty", regLabel3, "WM Uncertainty", 'Location','best')

subplot(2,1,2)
hold on

title("Voltage vs. Pitot-Static Probe Airspeeds")
rawPT2 = plot(final.voltage, final.airSpeedPitot, 'b.');
linePT2 = plot(final.voltage, LeastSquares2(final.voltage),'r');
uncertaintyPT2 = errorbar(final.voltage,LeastSquares2(final.voltage),sigmaVPitotPT,'LineStyle','none');

rawWM2 = plot(measureMan.pitotVoltage, measureMan.pitotAirspeed, 'm.');
lineWM2 = plot(measureMan.pitotVoltage, LeastSquaresMan2(measureMan.pitotVoltage),'g');
uncertaintyWM2 = errorbar(measureMan.pitotVoltage,LeastSquaresMan2(measureMan.pitotVoltage),sigmaVPitotMan,'LineStyle','none');

subset2 = [linePT2, uncertaintyPT2, lineWM2, uncertaintyWM2];

xlabel("Voltage (V)")
ylabel("Pitot-Static Airspeed (m/s)")
legend(subset2, regLabel2, "PT Uncertainty", regLabel4, "WM Uncertainty", 'Location','best')

hold off

% 1 plot, comparing airspeed models
figure
hold on

linePTFin1 = plot(final.voltage, LeastSquares1(final.voltage),'r');

lineWMFin1 = plot(measureMan.venturiVoltage, LeastSquaresMan1(measureMan.venturiVoltage),'b');

linePTFin2 = plot(final.voltage, LeastSquares2(final.voltage),'g');

lineWMFin2 = plot(measureMan.pitotVoltage, LeastSquaresMan2(measureMan.pitotVoltage),'m');

lineAvg = plot(final.voltage, LeastSquaresAvg(final.voltage),'k--');

subset3 = [linePTFin1, lineWMFin1, linePTFin2, lineWMFin2, lineAvg];

title("Voltage vs. Airspeeds")
regLabel5 = sprintf("v = %.3f*V + %.3f", x_hat_avg(1), x_hat_avg(2));
xlabel("Voltage (V)")
ylabel("Airspeed (m/s)")
legend(subset3, "Venturi + PT", "Venturi + WM", "Pitot-Static + PT", "Pitot-Static + WM", regLabel5);

hold off;








