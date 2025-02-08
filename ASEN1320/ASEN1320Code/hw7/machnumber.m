% Create a vector of size 20 filled with temperature values as specified 
% Hard-wiring of temperature values is not allowed (no point)

Tvec1 = -40:20:80;
Tvec2 = 100:100:1300;
Tvector = [Tvec1,Tvec2]; %%%%% CHECK THIS (2pt)

% Compute sound speed
% See Lecture 24 for multication of vector by a scalar
R = 287.058; % The specific gas constant for dry air [J/kg/K]
b = 1.4;     % The adiabatic coefficient 
SoundSpeedvector = sqrt(b*R*(Tvector+ 273.16)); %%%%% CHECK THIS (2pt)

% Generate a "random" array index from 1 to 20
% See Lecture 24 for randi
rng(uint64(now*1000))  % Set the random number seed with the current time
index = randi([1,20]); 
SoundSpeed = SoundSpeedvector(index); %%%%% CHECK THIS (2pt)

% Generate a "random" speed of an air vehicle object from 10 m/s from 1000 m/s
% See Lecture 24 for rand
VehicleSpeed = 10 + (1000-10)*rand; %%%%% CHECK THIS (2pt)

% Compute Mach number
MachNumber = VehicleSpeed/SoundSpeed; %%%%% CHECK THIS (2pt)

% Set up a string array
% See https://www.mathworks.com/help/matlab/ref/string.html
Regimes = ["Incompressible","Subsonic","Transonic","Sonic","Supersonic","Hypersonic"]; %%%%% CHECK THIS (2pt)

% See Lecture 25 for if-statement
% Recitation Week9

if (MachNumber < 0.3) 
   fprintf('%s flight regime and MachNumber is %f\n', Regimes(1), MachNumber)  %%%%% CHECK THIS OUTPUT (4 pt) if possible
elseif (MachNumber >= 0.3 && MachNumber < 0.8) 
   fprintf('%s flight regime and MachNumber is %f\n', Regimes(2), MachNumber)  
elseif (MachNumber >= 0.8 && MachNumber < 1) 
   fprintf('%s flight regime and MachNumber is %f\n', Regimes(3), MachNumber)  
elseif (MachNumber == 1) 
   fprintf('%s flight regime and MachNumber is %f\n', Regimes(4), MachNumber) 
elseif (MachNumber > 1 && MachNumber < 5) 
   fprintf('%s flight regime and MachNumber is %f\n', Regimes(5), MachNumber) 
else
   fprintf('%s flight regime and MachNumber is %f\n', Regimes(6), MachNumber) 
end
