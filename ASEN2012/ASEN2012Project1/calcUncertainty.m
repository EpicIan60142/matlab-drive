function sigma_Cs = calcUncertainty(Cc, Cs, Mc, sigma_Mc, Ms, sigma_Ms, T0, sigma_T0, T1, sigma_T1, T2, sigma_T2)
% Function that calculates the fnal specific heat uncertainty of the sample
% based on best estimates of other parameters.
%   Inputs: Calorimeter specific heat Cc, sample specific heat best estimate, 
%   calorimeter mass Mc, calorimeter mass uncertainty sigma_Mc, sample mass 
%   Ms, sample mass uncertainty sigma_Ms, initial calorimeter temperature T0,
%   initial calorimeter temperature uncertainty sigma_T0, initial sample 
%   temperature T1, initial sample temperature uncertainty sigma_T1, final 
%   calorimeter temperature T2, final calorimeter temperature uncertainty 
%   sigma_T2
%   Outputs: Sample specific heat uncertainty

%Calculate general rule uncertainties
sigma_Cs_Mc = (Cc*(T2-T0))/(Ms*(T1-T2));
sigma_Cs_T0 = -(Mc*Cc)/(Ms*(T1-T2));
sigma_Cs_Ms = -(Mc*Cc*(T2-T0))/((T1-T2)*(Ms)^2);
sigma_Cs_T2 = (Mc*Cc*(T1-T0))/(Ms*(T1-T2)^2);
sigma_Cs_T1 = -(Mc*Cc*(T2-T0))/(Ms*(T1-T2)^2);

%Calculate final uncertainty for the specific heat of sample B
sigma_Cs = Cs*sqrt((sigma_Cs_Mc*sigma_Mc)^2+(sigma_Cs_T0*sigma_T0)^2+(sigma_Cs_Ms*sigma_Ms)^2+(sigma_Cs_T2*sigma_T2)^2+(sigma_Cs_T1*sigma_T1)^2);

end
