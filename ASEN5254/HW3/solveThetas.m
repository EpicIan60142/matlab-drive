function thetas = solveThetas(point, testC1)
% Function that runs through the inverse kinematics procedure to solve for
% thetas for the 3 link manipulator in HW 3
%   Inputs:
%       - point: Coordinates of the end effector in the inertial frame
%       - testC1: Vector of cos(theta1) to calculate thetas for
%   Outputs:
%       - thetas: Calculated thetas for the provided cos(theta3) vector
%
%   By: Ian Faber, 09/18/2025
%

% Extract x and y
x = point(1);
y = point(2);

% Assign empty outputs
thetas.theta1 = [];
thetas.theta2 = [];
thetas.theta3 = [];

% Loop over cos(theta1) values
for k = 1:length(testC1)
        % Extract cos(theta1)
    c1 = testC1(k);

        % Calculate sin(theta1)
    for kk = 1:2
            % Positive or negative sin(theta1)
        if kk == 1
            s1 = sqrt(1-c1^2);
        else
            s1 = -sqrt(1-c1^2);
        end

            % Calculate xPrime and yPrime
        xPrime = x*c1 + y*s1 + 8*c1;
        yPrime = -x*s1 + y*c1 - 8*s1;

            % Calculate cos(theta3)
        c3 = (xPrime^2 + yPrime^2 - 145)/144;

            % Calculate sin(theta3)
        for idx = 1:2
                % Positive or negative sin(theta3)
            if idx == 1
                s3 = sqrt(1-c3^2);
            else
                s3 = -sqrt(1-c3^2);
            end

                % Calculate cos(theta2)
            c2 = (xPrime*(9*c3 + 8) + 9*yPrime*s3)/(xPrime^2 + yPrime^2);

                % Calculate sin(theta2)
            for idx2 = 1:2
                if idx2 == 1
                    s2 = sqrt(1-c2^2);
                else
                    s2 = -sqrt(1-c2^2);
                end

                    % Assign outputs
                theta1 = atan2(s1, c1);
                theta2 = atan2(s2, c2);
                theta3 = atan2(s3, c3);

                %     % Check result
                xCheck = cos(theta1)*(cos(theta2)*(9*cos(theta3) + 8) - 9*sin(theta2)*sin(theta3) + 8) ...
                       + sin(theta1)*(sin(theta2)*(-9*cos(theta3) - 8) - 9*cos(theta2)*sin(theta3));

                yCheck = cos(theta1)*(sin(theta2)*(9*cos(theta3) + 8) + 9*cos(theta2)*sin(theta3)) ...
                       + sin(theta1)*(cos(theta2)*(9*cos(theta3) + 8) - 9*sin(theta2)*sin(theta3) + 8);

                epsilon = 1e-10;
                if (xCheck <= point(1) + epsilon && xCheck >= point(1) - epsilon && yCheck <= point(2) + epsilon && yCheck >= point(2) - epsilon)
                    thetas.theta1 = [thetas.theta1; theta1];
                    thetas.theta2 = [thetas.theta2; theta2];
                    thetas.theta3 = [thetas.theta3; theta3];
                end
            end
        end
    end
end

end


