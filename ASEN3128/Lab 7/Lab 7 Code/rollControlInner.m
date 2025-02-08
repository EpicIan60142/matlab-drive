function ail_perturb = rollControlInner(p_c, p, ka)
    ail_perturb = ka*(p_c - p);
end