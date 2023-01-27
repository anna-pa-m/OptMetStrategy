function [smin, smax, ss, ps, extra] = two_states_environment(s0, s1, rho, delta, nrs, sigma, lb1, ub1, lb2, ub2)
    smin = 3e-4;
    smax = 1.0;
    ss = logspace(log10(smin),log10(smax), nrs+1);
    pd1 = makedist('Normal', 'mu', s0, 'sigma', sigma);
    pd2 = makedist('Normal', 'mu', s1, 'sigma', sigma);
    ps = rho * pdf(truncate(pd1, lb1, ub1),ss) + (1-rho) * pdf(truncate(pd2, lb2, ub2), ss);
    z = 0.0;
    for i = 1:numel(ss)-1
        z = z + ps(i) * (ss(i+1)-ss(i));
    end
    ps = ps ./ z;
    extra = '_s0_' + string(s0) + '_s1_' + string(s1);

end
