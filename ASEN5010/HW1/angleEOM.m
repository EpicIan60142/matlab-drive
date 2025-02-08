function [dAngle] = angleEOM(t, angle)
% Function that integrates 3-2-1 euler angles according to the euler angle
% kinematic differential equations IN DEGREES

yaw = angle(1);
theta = angle(2);
phi = angle(3);

w = [sin(0.1*t); 0.01; cos(0.1*t)]*20; % deg/s

dAngle = (1/cosd(theta))*[
                            0,           sind(phi),              cosd(phi);
                            0,           cosd(phi)*cosd(theta),  -sind(phi)*cosd(theta);
                            cosd(theta), sind(phi)*sind(theta),  cosd(phi)*sind(theta)
                         ]*w;






