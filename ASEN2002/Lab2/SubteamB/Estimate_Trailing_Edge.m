function [Cp_TE] = Estimate_Trailing_Edge(Cp_no_approx,Port_Locations)
% THIS FUNCTION ESTIMATES THE COEFFICIENT OF PRESSURE FOR THE TRAILING EDGE
% USING LINEAR EXTRAPOLATION OF BOTH THE TOP AND BOTTOM COEFFICIENTS OF
% PRESSURE 
% inputs -- the Cps for multiple port locations and those port locations

% find the chord length using the trailing edge distance
chord = max(Port_Locations(1,:));
% normalize each port location with respect to the chord line to non
% dimensionalize them -- makes it easier to calculate things
normalized_Locations = Port_Locations/chord;

% fit a line to the top coeffient of pressures -- linear extrapolation to get an
% approximate coeffiencet of pressure at the trailing edge
upper_approx = @(x) (((Cp_no_approx(9)-Cp_no_approx(8))/(normalized_Locations(1,9)-normalized_Locations(1,8)))*(x-(normalized_Locations(1,9)))) + Cp_no_approx(9);

% fit a line to the bottom coeffiecent of pressures -- linear extrapolation to get an
% approximate coefficent of pressure at the trailing edge
%lower_approx = @(x) (((Cp_no_approx(11)-Cp_no_approx(10)/(normalized_x(1,11)-normalized_x(1,10)))*(x-(normalized_x(1,11))) + Cp_no_approx(11)));
lower_approx = @(x) (((Cp_no_approx(10)-Cp_no_approx(11)/(normalized_Locations(1,10)-normalized_Locations(1,11)))*(x-(normalized_Locations(1,10))) + Cp_no_approx(10)));

upper = upper_approx(1);
lower = lower_approx(1);

% take the average of those two trailing edge coeeficent of pressure values
% to get the approximate coeffiecent of pressure at the trailing edge 
Cp_TE = (upper + lower)/2;
end
