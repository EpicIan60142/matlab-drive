function [c_l, c_dw] = DiamondAirfoil(M, alpha, epsilon1, epsilon2)
% DiamondAirfoil - Function that solves for the sectional coefficient of
% lift and coefficient of wave drag for an arbitrary diamond airfoil.
%
%   Inputs:
%       M: Freestream mach number
%       alpha: Angle of attack (deg)
%       epsilon1: Leading edge half-angle (deg)
%       epsilon2: Trailing edge half-angle (deg)
%
%   Outputs:
%       c_l: Sectional coefficient of lift
%       c_dw: Sectional coefficient of wave drag
%
%   Author: Ian Faber, 04/19/2023
%   Collaborators: Maggie Wussow

gamma = 1.4; % Specific heat ratio for air

bowShock = false; % Bow shock flag

M1 = M; % For notation's sake

% Starting static to total pressure ratio
[~, ~, p1p01, ~, ~] = flowisentropic(gamma, M1, 'mach');

l1c = 1/(cosd(epsilon1) + (sind(epsilon1)/tand(epsilon2))); % L1/c
l2c = 1/(cosd(epsilon2) + (sind(epsilon2)/tand(epsilon1))); % L2/c

% Determine case based on AoA and airfoil geometry
if alpha >= 0 % Positive or Zero AoA
    if epsilon1 < alpha
        state = 1;
    elseif epsilon1 == alpha
        state = 2;
    elseif epsilon1 > alpha
        state = 3;
    end
else % Negative AoA
    if abs(alpha) > epsilon1
        state = 4;
    elseif abs(alpha) == epsilon1
        state = 5;
    elseif abs(alpha) < epsilon1
        state = 6;
    end
end

for i = 1 % Trick for breaking out of case statement if a bow shock occurs
    switch state
        case 1 % Shock between 1 and 4, fans between 1 and 2, 2 and 3, 4 and 5
            % Shock between 1 and 4
                theta14 = alpha + epsilon1;
                beta14 = ObliqueShockBeta(M1, theta14, gamma, 'Weak');
    
                % If output is imaginary or negative, non-physical result
                if ~isreal(beta14) || beta14 < 0 
                    bowShock = true;
                    break;
                end

                [p4p1, p04p01, M4] = shock(beta14, theta14, M1, gamma);
    
            % Fan between 1 and 2
                theta12 = alpha - epsilon1;
                [p2p02, M2] = fan(theta12,M1, gamma);
    
            % Fan between 2 and 3
                theta23 = epsilon1 + epsilon2;
                [p3p03, ~] = fan(theta23, M2, gamma);
    
            % Fan between 4 and 5
                if M4 > 1  % Make sure M4 is supersonic, otherwise no fan
                    theta45 = epsilon1 + epsilon2;
                    [p5p05, ~] = fan(theta45, M4, gamma);
                end
    
            % Lift/Drag pressure ratios
                p2p1 = p2p02/p1p01; % isentropic, p02 = p01
                p3p1 = p3p03/p1p01; % isentropic, p03 = p01
                %p4p1 from shock
                if M4 > 1 % Fan calc was valid
                    p5p1 = p5p05*1*p04p01/p1p01; % p04 = p05, p5/p1 = p5/p05 * p05/p04 * p04/p01 * p01/p1
                else % No fan occurred
                    p5p1 = p4p1;
                end
    
        case 2 % Shock between 1 and 4, fans between 2(1) and 3, 4 and 5
            % Shock between 1 and 4
                 theta14 = 2*epsilon1;
                 beta14 = ObliqueShockBeta(M1, theta14, gamma, 'Weak');

                 % If output is imaginary or negative, non-physical result
                 if ~isreal(beta14) || beta14 < 0
                    bowShock = true;
                    break;
                end

                 [p4p1, p04p01, M4] = shock(beta14, theta14, M1, gamma);

            % Fan between 2 and 3
                theta23 = alpha + epsilon2;
                M2 = M1; % No shock or fan between 1 and 2
                [p3p03, ~] = fan(theta23, M2, gamma); 

            % Fan between 4 and 5
                if M4 > 1 % Make sure M4 is supersonic, otherwise no fan
                    theta45 = epsilon1 + epsilon2;
                    [p5p05, ~] = fan(theta45, M4, gamma);
                end

            % Lift/Drag pressure ratios
                p2p1 = 1; % No shock or fan between 1 and 2
                p3p1 = p3p03/p1p01;
                % p4p1 from shock
                if M4 > 1 % Fan calc was valid
                    p5p1 = p5p05*1*p04p01/p1p01;
                else % No fan occurred
                    p5p1 = p4p1;
                end
    
        case 3 % Shocks between 1 and 4, 1 and 2, fans between 2 and 3, 4 and 5
            % Shock between 1 and 4
                theta14 = alpha + epsilon1;
                beta14 = ObliqueShockBeta(M1, theta14, gamma, 'Weak');

                % If output is imaginary or negative, non-physical result
                if ~isreal(beta14) || beta14 < 0
                    bowShock = true;
                    break;
                end

                [p4p1, p04p01, M4] = shock(beta14, theta14, M1, gamma);

            % Shock between 1 and 2
                theta12 = epsilon1 - alpha;
                beta12 = ObliqueShockBeta(M1, theta12, gamma, 'Weak');

                % If output is imaginary or negative, non-physical result
                if ~isreal(beta12) || beta12 < 0
                    bowShock = true;
                    break;
                end

                [p2p1, p02p01, M2] = shock(beta12, theta12, M1, gamma);

            % Fan between 2 and 3
                if M2 > 1 % Make sure M2 is supersonic, otherwise no fan
                    theta23 = epsilon1 + epsilon2;
                    [p3p03, ~] = fan(theta23, M2, gamma);
                end

            % Fan between 4 and 5
                if M4 > 1 % Make sure M4 is supersonic, otherwise no fan
                    theta45 = epsilon1 + epsilon2;
                    [p5p05, ~] = fan(theta45, M4, gamma);
                end

            % Lift/Drag pressure ratios
                % p2p1 from shock
                if M2 > 1 % Fan calc was valid
                    p3p1 = p3p03*1*p02p01/p1p01;
                else % No fan occurred
                    p3p1 = p2p1;
                end

                % p4p1 from shock
                if M4 > 1 % Fan calc was valid
                    p5p1 = p5p05*1*p04p01/p1p01;
                else % No fan occurred
                    p5p1 = p4p1;
                end
    
        case 4 % Shock between 1 and 2, fans between 1 and 4, 2 and 3, 4 and 5
            % Shock between 1 and 2
                theta12 = abs(alpha) + epsilon1;
                beta12 = ObliqueShockBeta(M1, theta12, gamma, 'Weak');

                % If output is imaginary or negative, non-physical result
                if ~isreal(beta12) || beta12 < 0
                    bowShock = true;
                    break;
                end

                [p2p1, p02p01, M2] = shock(beta12, theta12, M1, gamma);

            % Fan between 1 and 4
                theta14 = abs(alpha) - epsilon1;
                [p4p04, M4] = fan(theta14, M1, gamma);

            % Fan between 2 and 3
                if M2 > 1 % Make sure M2 is supersonic, otherwise no fan
                    theta23 = epsilon1 + epsilon2;
                    [p3p03, ~] = fan(theta23, M2, gamma);
                end

            % Fan between 4 and 5
                theta45 = epsilon1 + epsilon2;
                [p5p05, ~] = fan(theta45, M4, gamma);

            % Lift/Drag pressure ratios
                % p2p1 from shock
                if M2 > 1 % Fan calc was valid
                    p3p1 = p3p03*1*p02p01/p1p01;
                else % No fan occurred
                    p3p1 = p2p1;
                end
                p4p1 = p4p04/p1p01;
                p5p1 = p5p05/p1p01;
    
        case 5 % Shock between 1 and 2, fans between 2 and 3, 4(1) and 5
            % Shock between 1 and 2
                theta12 = 2*epsilon1;
                beta12 = ObliqueShockBeta(M1, theta12, gamma, 'Weak');

                % If output is imaginary or negative, non-physical result
                if ~isreal(beta12) || beta12 < 0
                    bowShock = true;
                    break;
                end

                [p2p1, p02p01, M2] = shock(beta12, theta12, M1, gamma);
           
            % Fan between 2 and 3
                if M2 > 1 % Make sure M2 is supersonic, otherwise no fan
                    theta23 = epsilon1 + epsilon2;
                    [p3p03, ~] = fan(theta23, M2, gamma);
                end

            % Fan between 4 and 5
                theta45 = abs(alpha) + epsilon2;
                M4 = M1; % No shock or fan between 1 and 4
                [p5p05, ~] = fan(theta45, M4, gamma);

            % Lift/Drag pressure ratios
                % p2p1 from shock
                if M2 > 1 % Fan calc was valid
                    p3p1 = p3p03*1*p02p01/p1p01;
                else % No fan occurred
                    p3p1 = p2p1;
                end
                p4p1 = 1; % No shock or fan between 1 and 4
                p5p1 = p5p05/p1p01;
    
        case 6 % Shocks between 1 and 2, 1 and 4, fans between 2 and 3, 4 and 5
            % Shock between 1 and 2
                theta12 = abs(alpha) + epsilon1;
                beta12 = ObliqueShockBeta(M1, theta12, gamma, 'Weak');

                % If output is imaginary or negative, non-physical result
                if ~isreal(beta12) || beta12 < 0
                    bowShock = true;
                    break;
                end

                [p2p1, p02p01, M2] = shock(beta12, theta12, M1, gamma);

            % Shock between 1 and 4
                theta14 = epsilon1 - abs(alpha);
                beta14 = ObliqueShockBeta(M1, theta14, gamma, 'Weak');

                % If output is imaginary or negative, non-physical result
                if ~isreal(beta14) || beta14 < 0
                    bowShock = true;
                    break;
                end

                [p4p1, p04p01, M4] = shock(beta14, theta14, M1, gamma);

            % Fan between 2 and 3
                if M2 > 1 % Make sure M2 is supersonic, otherwise no fan
                    theta23 = epsilon1 + epsilon2;
                    [p3p03, ~] = fan(theta23, M2, gamma);
                end

            % Fan between 4 and 5
                if M4 > 1 % Make sure M4 is supersonic, otherwise no fan
                    theta45 = epsilon1 + epsilon2;
                    [p5p05, ~] = fan(theta45, M4, gamma);
                end

            % Lift/Drag pressure ratios
                % p2p1 from shock
                if M2 > 1 % Fan calc was valid
                    p3p1 = p3p03*1*p02p01/p1p01;
                else % No fan occurred
                    p3p1 = p2p1;
                end

                % p4p1 from shock
                if M4 > 1 % Fan calc was valid
                    p5p1 = p5p05*1*p04p01/p1p01;
                else % No fan occurred
                    p5p1 = p4p1;
                end
                
        otherwise
            fprintf("What???? \n") % Should never reach this :P
    end
    
    % Calculate coefficients depending on the sign of alpha
    if alpha >= 0
        c_l = (2/(gamma*M1^2))*(l1c*(p4p1*cosd(alpha+epsilon1) - p2p1*cosd(alpha-epsilon1)) + l2c*(p5p1*cosd(alpha-epsilon2) - p3p1*cosd(alpha+epsilon2)));
        c_dw = (2/(gamma*M1^2))*(l1c*(p4p1*sind(alpha+epsilon1) - p2p1*sind(alpha-epsilon1)) + l2c*(p5p1*sind(alpha-epsilon2) - p3p1*sind(alpha+epsilon2)));
    else
        alpha = abs(alpha);
        c_l = (2/(gamma*M1^2))*(l1c*(p4p1*cosd(alpha-epsilon1) - p2p1*cosd(alpha+epsilon1)) + l2c*(p5p1*cosd(alpha+epsilon2) - p3p1*cosd(alpha-epsilon2)));
        c_dw = (2/(gamma*M1^2))*(l1c*(p2p1*sind(alpha+epsilon1) - p4p1*sind(alpha-epsilon1)) + l2c*(p3p1*sind(alpha-epsilon2) - p5p1*sind(alpha+epsilon2)));
    end

end

% If result from oblique shock was non-physical, report it as a bow shock
if bowShock
    fprintf("Bow Shock occurred at M = %.0f and alpha = %.2f deg!\n", M, alpha)
    
    % Update outputs so function doesn't break
    c_l = nan;
    c_dw = nan;
end

end


