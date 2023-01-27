function [pq, peps, qlin, epslin] = get_q_eps_dist(a, ss, da, ds, q, eps, nqs, nes, ps, pags)

    disp('Computing p*(q) and p*(eps)');
    qmin = min(q);
    qmax = max(q);
    epsmin = min(eps);
    epsmax = max(eps);
    qlin = linspace(qmin,qmax,nqs);
    dq = (qmax - qmin)/numel(qlin);
    epslin = linspace(epsmin, epsmax, nes);
    deps = (epsmax - epsmin)/numel(epslin);

   %%% to be done
    for i=1:numel(ss)-1
        for j=1:numel(a)
            if isnan(pags(j,i))
                continue
            end
            pas(j,i) = pags(j,i) .* ps(i);
            prob_as(j,i) = pas(j,i) * da * ds(i);
            %%%% histogram q
            for l=1:numel(qlin)
                if (q(j) >= qlin(l) && q(j) < qlin(l) + dq)
                    Delta_q(j,i,l) = 1./dq;
                else
                    Delta_q(j,i,l) = 0.;
                end
            end
            for l=1:numel(epslin)
                if (eps(j) >= epslin(l) && eps(j) < epslin(l) + deps)
                    Delta_eps(j,i,l) = 1./deps;
                else
                    Delta_eps(j,i,l) = 0.;
                end
            end
        end
    end
    
    for l=1:numel(qlin)
        pq(l) = 0;
        for i=1:numel(ss)-1
            for j=1:numel(a)
                if isnan(pags(j,i))
                    continue
                end
                pq(l) = pq(l) + ds(i).*da.*pas(j,i) * Delta_q(j,i,l); 
            end
        end
    end
    for l=1:numel(epslin)
        peps(l) = 0;
        for i=1:numel(ss)-1
            for j=1:numel(a)
                if isnan(pags(j,i))
                    continue
                end
                peps(l) = peps(l) + ds(i).*da.*pas(j,i) * Delta_eps(j,i,l); 
            end
        end
    end
    


end