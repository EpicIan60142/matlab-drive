function rud_perturb = YawDamperControl(r_c, r, kr)
    rud_perturb = (r_c - r)*kr;
end