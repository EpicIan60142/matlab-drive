%Setup Constants
m  = 2 * 1.00784 * 1.660538921e-27; % H2 unit mass [kg]
kb = 1.38064852e-23;                % Boltzmann’s constant [J/K]
na = 6.02214076e+23;                %Avogadro's constant
c  = 1000.0;                        % Specific heat capacity [J/kg/K]
 
%User input
n    = input('Enter the number of particles: ');
vmin = input('Enter the minimum velocity in m/s: ');
vmax = input('Enter the maximum velocity in m/s: ');
Ts   = input('Enter the temperature of surrounding air in K: ');
 
rng('default'); 
velocity = rand(n,1)*(vmax-vmin) + vmin;
vsqsum   = sum(velocity.^2);
 
Tg = m * vsqsum / n / (3 * kb);       %Eq (1)
fprintf('The gas temperature in K: %d\n', round(Tg) )

Q = round(c * na * m * (Ts - Tg));    %Eq (2) 
fprintf('The transferred heat from the surrounding air to the gas in J: %d\n', Q);

IsHeatTransfer = Q > 0;
fprintf('Heat is transferred from the surrounding air to the gas: %d (1=true, 0=false)\n', ...
        IsHeatTransfer)
