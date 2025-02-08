function sigOut = addMRP(sig1, sig2)

% sig1 = sig', sig2 = sig''

sig1Mag = norm(sig1);
sig2Mag = norm(sig2);

denom = 1 + sig1Mag^2*sig2Mag^2 - 2*dot(sig1, sig2);
num = (1-sig1Mag^2)*sig2 + (1-sig2Mag^2)*sig1 - 2*cross(sig2,sig1);

sigOut = num./denom;

end