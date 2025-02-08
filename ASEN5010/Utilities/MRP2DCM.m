function mat = MRP2DCM(sigma)

mat = eye(3) + (1/(1+norm(sigma)^2)^2)*(8*tilde(sigma)*tilde(sigma) - 4*(1-norm(sigma)^2)*tilde(sigma));

end