function [plam, Istar, fstar, aav] = get_mu_dist(a, ss, da, ds, f, pags, qa, lam, dlam, ps)

    disp('Computing I*, \mu* and histogram for p*(\mu)');
    Istar = 0;
    fstar = 0;
    for i=1:numel(ss)-1
        aav(i) = 0.;
        for j=1:numel(a)
            if isnan(pags(j,i))
                aav(i) = NaN;
                pas(j,i) = NaN;
                continue
            end
            aav(i) = aav(i) + a(j).*pags(j,i)*da;
            % I = \int ds da P(a,s) log [P(a|s)/P(a)]
            if ~isnan(log2(pags(j,i)./qa(j))) & ~isinf(log2(pags(j,i)./qa(j)))
                Istar = Istar+ds(i).*ps(i).*da.*pags(j,i).*log2(pags(j,i)./qa(j));
            end
            pas(j,i) = pags(j,i) .* ps(i);
            prob_as(j,i) = pas(j,i) * da * ds(i);
            % < \mu >_{p(mu, s)}
            fstar = fstar + ds(i).*da.*pas(j,i).*f(j,i);
            %%%%%%%%%%%%%%%%%%%%%%%% histogram \mu
            for l=1:numel(lam)
                if (f(j,i) >= lam(l) && f(j,i) < lam(l)+dlam)
                    Delta(j,i,l)=1./dlam;
                else
                    Delta(j,i,l)=0.;
                end
            end
        end
    end
    disp('I* ='); disp(Istar);
    disp('\mu* ='); disp(fstar);
    disp('Computing p*(\mu)');

    normplam=0.;
    for l=1:numel(lam)
        plam(l)=0;
        for i=1:numel(ss)-1
            for j=1:numel(a)
                if isnan(pags(j,i))
                    continue
                end
                plam(l)=plam(l)+ds(i).*da.*pas(j,i) * Delta(j,i,l); % P(s,\alpha)
            end
        end
        normplam=normplam+plam(l).*dlam;
    end

    avglam=0.;
    for l=1:numel(lam)
        avglam=avglam+lam(l).*plam(l).*dlam;
    end
    disp('avg \mu'); disp(avglam);


end