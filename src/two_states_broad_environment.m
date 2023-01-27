function [smin, smax, ss, ps, extra] = two_states_broad_environment(s0, s1, rho, nrs, sigma, unip)
    smin = 3e-4;
    smax = 1.0;
    ss = logspace(log10(smin),log10(smax), nrs+1);
    ps = rho .* normpdf(ss, s0, sigma) + (1-rho) .* normpdf(ss, s1,sigma);
    ps = ps + unip .* unifpdf(ss);
    z = 0.0;
    for i = 1:numel(ss)-1
        z = z + ps(i) * (ss(i+1)-ss(i));
    end
    ps = ps ./ z;
    extra = '_s0_' + string(s0) + '_s1_' + string(s1);

end
