function BN = OLAE(Vb, Vn, w)
% Vb is a list of body vectors, Vn is a list of inertial vectors, and w is
% their weights

% Normalize unit vectors
for k = 1:size(Vb,2)
    Vb(:,k) = Vb(:,k)/norm(Vb(:,k));
end

for k = 1:size(Vn,2)
    Vn(:,k) = Vn(:,k)/norm(Vn(:,k));
end

% Calculate s and d vectors
for k = 1:size(Vb,2)
    sVec(:,k) = Vb(:,k) + Vn(:,k);
    dVec(:,k) = Vb(:,k) - Vn(:,k);
end

d = [];
for k = 1:size(dVec,2)
    d = [d; dVec(:,k)];
end

S = [];
for k = 1:size(sVec,2)
    S = [S; tilde(sVec(:,k))];
end

weights = [];
for k = 1:length(w)
    weights = [weights; w(k)*ones(3,1)];
end
W = diag(weights);

q = ((S'*W*S)^-1)*S'*W*d;

BN = CRP2DCM(q);

end