function dMRP = MRPEOM(t, MRP)

sig1 = MRP(1);
sig2 = MRP(2);
sig3 = MRP(3);

sig = [sig1; sig2; sig3];

sigNorm = norm(sig)^2;

% if norm(sig) >= 1
%     sig = -sig/sigNorm;
% end

sigNorm = norm(sig)^2;

w = deg2rad([sin(0.1*t); 0.01; cos(0.1*t)]*20); % rad/s

dMRP = 0.25*((1-sigNorm)*eye(3) + 2*tilde(sig) + 2*(sig*sig'))*w;

end