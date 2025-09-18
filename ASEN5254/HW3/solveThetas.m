function thetas = solveThetas(point, testC3)
% Function that runs through the inverse kinematics procedure to solve for
% thetas for the 3 link manipulator in HW 3
%   Inputs:
%       - point: Coordinates of the end effector in the inertial frame
%       - testC3: Vector of cos(theta3) to calculate thetas for
%   Outputs:
%       - thetas: Calculated thetas for the provided cos(theta3) vector
%
%   By: Ian Faber, 09/18/2025
%

% Extract x and y
x = point(1);
y = point(2);

% Function handles for fsolve
s2Func = @(c2) sqrt(1-c2^2);
c2Func = @(c2,c3,s3) (x^2 + y^2 - 144*c3 - 209 + 144*s2Func(c2)*s3)/(144*c3 + 128) - c2;

s1Func = @(c1) sqrt(1-c1^2);
c1Func = @(c1,s2,c2,s3,c3) (x - s1Func(c1)*(s2*(-9*c3 - 8) -9*c2*s3))/(c2*(9*c3 + 8) - 9*s2*s3 + 8) - c1;

% Loop through all cos(theta3)'s
thetas.theta1 = [];
thetas.theta2 = [];
thetas.theta3 = [];
for k = 1:length(testC3)
        % Choose cos(theta3)
    c3 = testC3(k);

        % Positive and negative sin(theta3)
    for kk = 1:2
        if kk == 1
            s3 = sqrt(1-(c3^2));
                % Rootsolve for cos(theta2)
            c2 = fsolve(@(c2)c2Func(c2, c3, s3), 1);

            if abs(c2) > 1
                continue;
            end

                % Positive and negative sin(theta2)
            for idx = 1:2
                if idx == 1
                    s2 = sqrt(1-(c2^2));
                        % Rootsolve for cos(theta1)
                    c1 = fsolve(@(c1)c1Func(c1, s2, c2, s3, c3), 0);

                    if abs(c1) > 1
                        continue;
                    end

                        % Positive and negative sin(theta1)
                    for idx2 = 1:2
                        if idx2 == 1
                            s1 = sqrt(1-(c1^2));
                        else
                            s1 = -sqrt(1-(c1^2));
                        end
                    end
                else
                    s2 = -sqrt(1-(c2^2));
                        % Rootsolve for cos(theta1)
                    c1 = fsolve(@(c1)c1Func(c1, s2, c2, s3, c3), 0);

                    if abs(c1) > 1
                        continue;
                    end

                        % Positive and negative sin(theta1)
                    for idx2 = 1:2
                        if idx2 == 1
                            s1 = sqrt(1-(c1^2));
                        else
                            s1 = -sqrt(1-(c1^2));
                        end
                    end
                end
            end
        else
            s3 = -sqrt(1-(c3^2));

                % Rootsolve for cos(theta2)
            c2 = fsolve(@(c2)c2Func(c2, c3, s3), 0);

            if abs(c2) > 1
                continue;
            end

                % Positive and negative sin(theta2)
            for idx = 1:2
                if idx == 1
                    s2 = sqrt(1-(c2^2));

                        % Rootsolve for cos(theta1)
                    if (c2*(9*c3 + 8) - 9*s2*s3 + 8) == 0
                        continue;
                    end

                    c1 = fsolve(@(c1)c1Func(c1, s2, c2, s3, c3), 0);

                    if abs(c1) > 1
                        continue;
                    end

                        % Positive and negative sin(theta1)
                    for idx2 = 1:2
                        if idx2 == 1
                            s1 = sqrt(1-(c1^2));
                        else
                            s1 = -sqrt(1-(c1^2));
                        end
                    end
                else
                    s2 = -sqrt(1-(c2^2));
                        % Rootsolve for cos(theta1)
                    if (c2*(9*c3 + 8) - 9*s2*s3 + 8) == 0
                        continue;
                    end

                    c1 = fsolve(@(c1)c1Func(c1, s2, c2, s3, c3), 0);
                    
                    if abs(c1) > 1
                        continue;
                    end

                        % Positive and negative sin(theta1)
                    for idx2 = 1:2
                        if idx2 == 1
                            s1 = sqrt(1-(c1^2));
                        else
                            s1 = -sqrt(1-(c1^2));
                        end
                    end
                end
            end
        end
    end

    try
            % Invalidate answer if imaginary angles
        if any(abs(imag([s1, c1, s2, c2, s3, c3])) > 0)
            continue;
        end
    
            % Calculate outputs
        theta1 = atan2(s1, c1);
        theta2 = atan2(s2, c2);
        theta3 = atan2(s3, c3);
    
            % Assign outputs
        thetas.theta1 = [thetas.theta1; theta1];
        thetas.theta2 = [thetas.theta2; theta2];
        thetas.theta3 = [thetas.theta3; theta3];
    catch
        continue;
    end

end


