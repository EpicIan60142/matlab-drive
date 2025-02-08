function EP = DCM2EP(C)
    % Sheppard's method
    q0 = sqrt(0.25*(1 + trace(C)));
    q1 = sqrt(0.25*(1 - trace(C) + 2*C(1,1)));
    q2 = sqrt(0.25*(1 - trace(C) + 2*C(2,2)));
    q3 = sqrt(0.25*(1 - trace(C) + 2*C(3,3)));
    
    q = [q0, q1, q2, q3];
    
    [~, idx] = max(q);
    
    if idx == 1 % q0 was largest, can divide by it safely
        q1 = (C(2,3)-C(3,2))/(4*q0);
        q2 = (C(3,1)-C(1,3))/(4*q0);
        q3 = (C(1,2)-C(2,1))/(4*q0);
        % q0 = (C(2,3)-C(3,2))/(4*q1); % Check sign on q0
    elseif idx == 2 % q1 was largest, can divide by it safely
        q0 = (C(2,3)-C(3,2))/(4*q1);
        q2 = (C(1,2)+C(2,1))/(4*q1);
        q3 = (C(3,1)+C(1,3))/(4*q1);
        % q1 = (C(2,3)-C(3,2))/(4*q0); % Check sign on q1
    elseif idx == 3 % q2 was largest, can divide by it safely
        q0 = (C(3,1)-C(1,3))/(4*q2);
        q1 = (C(1,2)+C(2,1))/(4*q2);
        q3 = (C(2,3)+C(3,2))/(4*q2);
        % q2 = (C(3,1)-C(1,3))/(4*q0); % Check sign on q2
    else
        q0 = (C(1,2)-C(2,1))/(4*q3);
        q1 = (C(3,1)+C(1,3))/(4*q3);
        q2 = (C(2,3)+C(3,2))/(4*q3);
        % q3 = (C(1,2)-C(2,1))/(4*q0); % Check sign on q3
    end
    
    EP = [q0, q1, q2, q3]';

end
