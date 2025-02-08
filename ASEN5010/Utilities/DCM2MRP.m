function sigma = DCM2MRP(mat, short)

% q = DCM2CRP(mat)
% 
% sig = q/(1+sqrt(1+q'*q));
% 
% sigMag = norm(sig)^2
% 
% if short
%     if sigMag <= 1
%         sigma = sig;
%     else
%         sigma = -sig/(sigMag^2);
%     end
% else
%     if sigMag <= 1
%         sigma = -sig/(sigMag^2);
%     else
%         sigma = sig;
%     end
% end

zeta = sqrt(trace(mat) + 1);

if zeta == 0
    EP = DCM2EP(mat);
    sig = [EP(2)/(1+EP(1)); EP(3)/(1+EP(1)); EP(4)/(1+EP(1))];
else
    sig = (1/(zeta*(zeta + 2)))*[mat(2,3) - mat(3,2); mat(3,1) - mat(1,3); mat(1,2) - mat(2,1)];
end

sigMag = norm(sig)^2;

if short
    if sigMag <= 1
        sigma = sig;
    else
        sigma = -sig/(sigMag^2);
    end
else
    if sigMag <= 1
        sigma = -sig/(sigMag^2);
    else
        sigma = sig;
    end
end

end