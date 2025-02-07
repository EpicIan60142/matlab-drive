function [time, theta, w] = UWData(filename)
    data = readmatrix(filename);
    theta = data(:,2);
    good = theta>0.5 & theta<15;
    time = data(good,1);
    theta = data(good,2);
    w = data(good,3);
end