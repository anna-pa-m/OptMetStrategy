function [smin, smax, ss, ps, extra] = exp_environment(s0, delta, nrs)


    smin=s0.*power(10,-2.0.*delta);
    smax=s0.*power(10, 1.0.*delta);
    ss = logspace(log10(smin),log10(smax), nrs+1);
    for i = 1:numel(ss)
        ps(i) = exp(-ss(i)/s0);
    end
    norm = 0.;
    for i = 1:numel(ss)-1
        s(i) = 0.5*(ss(i+1) + ss(i));
        ds(i) = ss(i+1) - ss(i);
        % normalization
        norm = norm + ps(i)*ds(i);
    end
    ps = ps ./ norm;
    extra = '_s0_' + string(s0);

end