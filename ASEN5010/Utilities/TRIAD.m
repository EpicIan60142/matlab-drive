function BN = TRIAD(v1b, v2b, v1n, v2n)
% Assumes v1 is the more accurate sensor

t1b = v1b;
t2b = cross(v1b, v2b);
t2b = t2b/norm(t2b);
t3b = cross(t1b, t2b);

t1n = v1n;
t2n = cross(v1n, v2n);
t2n = t2n/norm(t2n);
t3n = cross(t1n, t2n);

BT = [t1b, t2b, t3b];
NT = [t1n, t2n, t3n];

BN = BT*NT';

end