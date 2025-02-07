clear all; clc; close all;

%what file do you want to run (1-3)
indices = [7,6,5];
for k = 1:3
    fileNum = k;
    
    %load data
    fileValues = ["6" "9" "12"];
    fileName = append('Data\Sec12_Group6_temp',fileValues(fileNum));
    data = load(fileName);
    
    %find start and end of data using pulse voltage (1,0)
    startInd = find(data(:,8),1,'first');
    endInd = find(data(:,8),1,'last');
    
    %removing start and end of data
    data = data(startInd:endInd,:);
    % data(:,1) = data(:,1) - data(1,1);
    
    %find rpm using tachorpm which measures difference between pulses
    dt = mean(diff(data(:,1)));
    hz = 1/dt;
    rpm = mean(tachorpm(data(:,8),hz));
    
    % find plate temperatures avg of top and bottom couple per plate
    TopPlateTemp = mean([mean(data(:,3)) mean(data(:,4))]);
    BotPlateTemp = mean([mean(data(:,5)) mean(data(:,6))]);
    
    diamPowerPiston = 19.5 / 1000; %mm
    PPcrossArea = pi*(diamPowerPiston/2)^2; %m^2
    PPtravel = 11/1000; %mm
    PPvolume = PPcrossArea*PPtravel; %m^3
    
    chamberHeight = 21 / 1000; %mm
    chamberDiam = 150 / 1000; %mm
    chamberVolume = chamberHeight*pi*(chamberDiam/2)^2; %m^3
    
    displacerDiam = 144/1000; %mm
    displacerHeight = 11/1000; %mm
    displacerVolume = pi*(displacerDiam/2)^2 * displacerHeight; %m^3
    
    TotalVolumeMax = chamberVolume - displacerVolume + PPvolume; %m^3
    TotalVolumeMin = chamberVolume - displacerVolume; %m^3
    
    %load data
    fileValues = ["6" "9" "12"];
    fileName = append('Data\solidWorksVolume_',fileValues(fileNum),'.xlsx');
    volData = xlsread(fileName);
    
    volData(:,3) = volData(:,3) - min(volData(:,3));
    
    %mm to meters
    PPdispVol = (volData(:,3)/1000).*PPcrossArea;
    
    zeroP = find(diff(sign(data(:,2))),3,'first');% + indices(k);
    
    zeroV = find(diff(sign(PPdispVol-mean(PPdispVol))),3,'first') + indices(k);
    
    volume = PPdispVol(zeroV(1):zeroV(3));
    pressure = data(zeroP(1):zeroP(3),2)*6894.76; % psi to Pa
    
    Xq = linspace(1,length(pressure),length(volume));
    Vq = interp1(1:length(pressure),pressure,Xq, 'linear','extrap');
    
    figure
    hold on
    titleText = sprintf("PV Diagram for \\Delta T = %s ^oC", fileValues(fileNum));
    title(titleText)
    plot(volume, Vq)
    xlabel("Volume [m^3]")
    ylabel("Pressure [Pa]")
    
    
    time = data(:,1);
    time = time(1:length(pressure));
    current = data(:,7);
    current = current(1:length(pressure));
%     current = mean(current);
    V = 48;
%     Qin(k) = V*current*time(end)

    Qin(k) = V*trapz(time,current);
    
    % Qin = 2.0635;
    
    minVol = find(volume == min(volume),1,'first');
    maxVol = find(volume == max(volume),1,'first');
    Wout(k) = trapz(volume(minVol:maxVol), Vq(minVol:maxVol));
    Win(k) = trapz(volume(1:minVol),Vq(1:minVol)) + trapz(volume(maxVol:end),Vq(maxVol:end));
    netWork(k) = trapz(volume, Vq);%Wout(k) - Win(k);
    
    efficiency(k) = (netWork(k) / Qin(k)) * 100;

    Qout(k) = Qin(k) - netWork(k);
end
