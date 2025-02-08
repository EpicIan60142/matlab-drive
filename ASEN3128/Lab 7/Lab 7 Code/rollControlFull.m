function ail_perturb = rollControlFull(phi_c, phi, p, ka, kp)
    ail_perturb = ka*(kp*(phi_c - phi) - p);
end