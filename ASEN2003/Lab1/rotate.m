function newCoords = rotate(oldCoords, yaw, pitch, roll)
    % Yaw rotates about z axis
    % Pitch rotates about y axis, 
    % Roll rotates about x axis
    % Angles taken in degrees
    
    pitch = pi*(pitch/180);
    yaw = pi*(yaw/180);
    roll = pi*(roll/180);

    rotateMatrix = [
                    cos(yaw)*cos(pitch), cos(yaw)*sin(pitch)*sin(roll) - sin(yaw)*cos(roll), cos(yaw)*sin(pitch)*cos(roll) + sin(yaw)*sin(roll);
                    sin(yaw)*cos(pitch), sin(yaw)*sin(pitch)*sin(roll) + cos(yaw)*cos(roll), sin(yaw)*sin(pitch)*cos(roll) - cos(yaw)*sin(roll);
                    -sin(pitch)        , cos(pitch)*sin(roll)                              , cos(pitch)*cos(roll)
                   ];
    
    newCoords = rotateMatrix*oldCoords;
end
