%% Homework 7 - Mach Number and Flow Regimes
clc

%Temperature vector
TvectorLow = -40:20:100; %Lower temperature vector (separated by 20 degC
TvectorHigh = 200:100:1300; %Higher temperature vector separated by 100 degC
Tvector = [TvectorLow  TvectorHigh]; %Concatenate the lower and higher temperature vectors into a single temperature vector

%Sound speed vector
b = 1.4; %Adiabatic coefficient
R = 287.058; %Specific gas constant
SoundTvector = Tvector + 273.16; %Convert temperature vector from degC to K
SoundSpeedvector = sqrt(b*R*SoundTvector); %Calculate sound speeds based on the temperature vector

%Random sound speed
rng(uint64(now*1000)); %Seed the random number generator
index = randi([1,20], 1, 1); %Generate a random array index between 1 and 20
SoundSpeed = SoundSpeedvector(index) %Select a random sound speed based on the random index

%Random vehicle speed
VehicleSpeed = 10 + (1000-10)*rand %Generate a random vehicle speed between 10 and 1000

%Calculate the mach number
MachNumber = VehicleSpeed / SoundSpeed 

%Array of flight regimes
Regimes = ["Incompressible", "Subsonic", "Transonic", "Sonic", "Supersonic", "Hypersonic"];

%Output logic
if(MachNumber < 0.3)
    fprintf("%s flight regime and MachNumber is %f \n", Regimes(1), MachNumber) %Report the regime as incompressible
elseif(MachNumber >= 0.3 && MachNumber < 0.8)
    fprintf("%s flight regime and MachNumber is %f \n", Regimes(2), MachNumber) %Report the regime as subsonic
elseif(MachNumber >= 0.8 && MachNumber < 1)
    fprintf("%s flight regime and MachNumber is %f \n", Regimes(3), MachNumber) %Report the regime as transonic
elseif(MachNumber == 1)
    fprintf("%s flight regime and MachNumber is %f \n", Regimes(4), MachNumber) %Report the regime as sonic
elseif(MachNumber > 1 && MachNumber < 5)
    fprintf("%s flight regime and MachNumber is %f \n", Regimes(5), MachNumber) %Report the regime as supersonic
else
    fprintf("%s flight regime and MachNumber is %f \n", Regimes(6), MachNumber) %Report the regime as hypersonic
end
