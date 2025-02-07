function [avgISP, stdISP] = analyzeISP(meanMWater, stdMWater, nSims)
% Function that finds average ISP for a bottle rocket
%   Inputs: None
%   Outputs: Average ISP, ISP standard deviation

    sampFreq = 1652; % Sampled at 1.652 kHz
    maxTime = 5; % Maximum allowable time for data start and stop calcs
    relError = 0.05; % Relative error between thrust averages for start and stop calcs
    g0 = 9.81; % Freefall acceleration on Earth
    mWater = meanMWater + stdMWater*randn(nSims,1); % Mass of water propellant (1000 g)

    avgISP = zeros(nSims, 1);
    stdISP = zeros(nSims, 1);

    %% Data Extraction
    numFiles = 16;

    for i = 1:numFiles
        files{i} = dir(['..\Static Test Stand Data\Fixed Mass\LA_Test_FixedMass_Trial',num2str(i)]);
        files{i}.path = ['..\Static Test Stand Data\Fixed Mass\LA_Test_FixedMass_Trial',num2str(i)];
    end

    %% Data Processing
    for i = 1:numFiles
       files{i}.data = 4.44822*load(files{i}.path); %Convert from lbf to N
       files{i}.thrust = files{i}.data(:,3);
       files{i}.peakThrust = max(files{i}.thrust);
       files{i}.time = (linspace(0,length(files{i}.thrust)/sampFreq,length(files{i}.thrust)))'; % Compute time based on sampling frequency

       % Calculate a 5-value moving average for detecting trend changes
       thrustAvg = (files{i}.thrust(1:end-4) + files{i}.thrust(2:end-3) + files{i}.thrust(3:end-2) + files{i}.thrust(4:end-1) + files{i}.thrust(5:end))/5;

       % Check for abrupt changes in the data and if time is less than 5 seconds
       conditionStart = ischange(thrustAvg) & files{i}.time(1:end-4) < maxTime;

       % Check to see if data changes by less than relError*100% and if time is less than 5 seconds
       conditionStop = thrustAvg(1:end-1)./thrustAvg(2:end) - 1 > relError & files{i}.time(1:end-5) < maxTime;

       start(i) = find(conditionStart, 1, 'first');
       stop(i) = find(conditionStop, 1, 'last');

       files{i}.thrustTime = files{i}.time(stop(i))- files{i}.time(start(i)); 
    end

    for k = 1:nSims
        %fprintf("Setting up simulation %d of %d... Standby...\n", k, nSims);
        
        for i = 1:numFiles
            files{i}.waterFlow = mWater(k)/files{i}.thrustTime; % Calculate mass flow rate of the water
            files{i}.impulse = trapz(files{i}.time(start(i):stop(i)), files{i}.thrust(start(i):stop(i)) - files{i}.waterFlow.*files{i}.time(start(i):stop(i)));
            files{i}.isp = files{i}.impulse./(mWater(k)*g0);

            isp(i,1) = files{i}.isp;
        end
        
        avgISP(k) = mean(isp);
        stdISP(k) = std(isp);
    end
end

