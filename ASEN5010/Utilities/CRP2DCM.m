function mat = CRP2DCM(q)

mat = (1/(1 + q'*q)) * ( (1-q'*q)*eye(3) + 2*(q*q') - 2*tilde(q));

end