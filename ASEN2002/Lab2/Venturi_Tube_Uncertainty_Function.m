function sigma_Venturi = Venturi_Tube_Uncertainty_Function(V_Venturi, Tatm, Patm, R, A1, A2, sigma_deltaP, sigma_Patm)
% This function is used to calculate the uncertainty associated with the
% Venturi_Tube velocity.
%   
% Uncertainty is calculated through the use of the general method along
% with partial derivatives of deltaP, Tatm, and Patm.

deltaP = 1/2*(Patm./(R*Tatm))*(1-(A2/A1)^2).*V_Venturi.^2; % Calculation of deltaP via equation found in appendix

syms a b c
    %a = deltaP
    %b = Tatm
    %c = Patm

EQN = sqrt((2.*a.*R.*b)./(c.*(1-(A2./A1).^2)));

sigma_Tatm = 0.25;

p_deltaP(a,b,c) = diff(EQN,a);
p_deltaP_val = p_deltaP(deltaP,Tatm,Patm);

p_Tatm(a,b,c) = diff(EQN,b); 
p_Tatm_val = p_Tatm(deltaP,Tatm,Patm);

p_Patm(a,b,c) = diff(EQN,c); 
p_Patm_val = p_Patm(deltaP,Tatm,Patm);

partials = [p_deltaP_val, p_Tatm_val, p_Patm_val];
sigmas = [sigma_deltaP, sigma_Tatm, sigma_Patm];

sigma = zeros(1,length(V_Venturi));
for i = 1:length(V_Venturi)
    temp = vpa(sqrt(sum((partials(i).*sigmas).^2)),5);
    sigma(i) = temp(1);
end
%sigma = vpa(sqrt(sum((partials.*sigmas).^2)),5) % Something funky is happening with the matrixes, when ran in program it spits out 3 sigmas, will see if I can fix this

sigma_Venturi = sigma/100;

end

