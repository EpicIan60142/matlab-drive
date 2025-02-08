function q = DCM2CRP(mat)

% det1 = (1 + mat(1,1)) / (mat(2,3) - mat(3,2));
% if det1 < 0
%     q1 = det1 + sqrt(det1^2 - 1);
% else
%     q1 = det1 - sqrt(det1^2 - 1);
% end
% 
% det2 = (1 + mat(2,2)) / (mat(3,1) - mat(1,3));
% if det2 < 0
%     q2 = det2 + sqrt(det2^2 - 1);
% else
%     q2 = det2 - sqrt(det2^2 - 1);
% end
% 
% det3 = (1 + mat(3,3)) / (mat(1,2) - mat(2,1));
% if det3 < 0
%     q3 = det3 + sqrt(det3^2 - 1);
% else
%     q3 = det3 - sqrt(det3^2 - 1);
% end
% 
% q = [q1; q2; q3];

zeta = sqrt(trace(mat) + 1);

q = 1/(zeta^2)*[mat(2,3) - mat(3,2); mat(3,1) - mat(1,3); mat(1,2) - mat(2,1)];


end