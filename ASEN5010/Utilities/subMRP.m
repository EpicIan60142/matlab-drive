function sigOut = subMRP(sig, sig1)

% sig = sig, sig1 = sig'

sigMag = norm(sig);
sig1Mag = norm(sig1);

denom = 1 + sig1Mag^2*sigMag^2 + 2*dot(sig1, sig);
num = (1-sig1Mag^2)*sig - (1-sigMag^2)*sig1 + 2*cross(sig,sig1);

sigOut = num./denom;

% May need to invert final result, track frames

end