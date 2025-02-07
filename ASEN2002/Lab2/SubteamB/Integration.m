function [Cl, Cd] = Integration(Cps, normalized_Locations, Final)
% THIS FUNCITON INTEGRATES THE COEFFICIENTS OF PRESSURE AND USES THAT TO 
% GET US  COEFFIENCETS OF AXIAL AND NORMAL STRESS WHICH ARE CONVERTED TO
% COEFFIENCTS OF LIFT AND DRAG
% inputs-- the Coefficients of pressure, the normalized locations of the
% ports relative to the chord line and the struct final which has the
% angles of attack

% find out delata x and delta y from our normalized locations matrix
delta_x = diff(normalized_Locations(1,:));
delta_y = diff(normalized_Locations(2,:));

% for loop looping through each different angle of attack and velocity
% matched with it
for i = 1:length(Cps)
    Cp_avg = zeros(1,length(Cps(1,:))-1); 
    % loops through each port 
    for j = 1:length(Cps(1,:))-1
        Cpis = Cps(i,:);
        Cp_avg(j) = (Cpis(j)+Cpis(j+1))/2;
    
    end
    
    % get coefficients of normal and axial at a point by taking the average of 2
    % coeeficents of pressure and multiplying by delta x or y -- EQN IN 7.4
    Cns = Cp_avg .* delta_x;
    Cas = Cp_avg .* delta_y;
    
    % find the total coefficient of normal and axial by summing them all
    Cn(i) = -1.0*sum(Cns);
    Ca(i) = sum(Cas);
    
    % use the coeffecients of normal and axial and convert them to the
    % coefficents of lift and drag using EQN in 7.4 with angle of attack
    Cl(i) = Cn(i)*cosd(Final.AngAttack(i)) - Ca(i)*sind(Final.AngAttack(i));
    Cd(i) = Ca(i)*cosd(Final.AngAttack(i)) + Cn(i)*sind(Final.AngAttack(i));

end


end
