function [Ip, FB] = princInertia(Ib)

[FB, Ip] = eig(Ib);
FB = FB';

end