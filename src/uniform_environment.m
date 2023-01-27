function [smin, smax, ss, ps, extra] = uniform_environment(s0, nrs)
    
    smin = s0 - 5 / 100 * s0;
    smax = s0 + 5 / 100 * s0;
    ss = logspace(log10(smin),log10(smax), nrs+1);
    ps = unifpdf(ss, smin, smax);
    extra = '_s0_' + string(s0);
end