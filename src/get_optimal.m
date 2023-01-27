function [qa, pags, av_alpha] = get_optimal(a, ss, da, ds, ps, f, beta)
    
    
    disp('Beta '); disp(beta);
    for j=1:numel(a)
        qa(j) = 1.;
        for i = 1:numel(ss)-1
            pags(j,i) = 1.0;
        end
    end

    disp('Finding optimal channel, iter n.');
    th = 1e-2;
    err_pa = 1;
    err_qa = 1;
    k = 1;
    while err_pa > th & err_qa > th
        err_qa = 0;
        err_pa = 0;
        for i=1:numel(ss)-1
            z(i)=0;
            % N(s ; \beta)
            for j=1:numel(a)
                if isnan(f(j,i))
                    pags(j,i) = NaN;
                    continue
                end
                z(i)=z(i)+da.*qa(j).*exp(beta.*f(j,i));
            end
            %disp(z(i))
        end
        % compute P(\alpha, s) and update Q(\alpha)
        for j=1:numel(a)
            qanew(j)=0.;
            for i=1:numel(ss)-1
                if isnan(f(j,i))
                    continue
                end
                err_pa = max(err_pa, abs(qa(j).*exp(beta.*f(j,i))./z(i) - pags(j,i)));
                pags(j,i) = qa(j).*exp(beta.*f(j,i))./z(i);
                qanew(j) = qanew(j)+ds(i).*ps(i).*pags(j,i);
            end
            err_qa = max(err_qa, abs(qa(j) - qanew(j)));
            qa(j) = qanew(j);
        end
        k = k + 1;
        if(mod(k, 1000) == 0)
            fprintf("steps: %d err Q: %f err P: %f\n", k, err_qa, err_pa);
        end
    end

    %for s = 1:numel(ss)-1
    %    assert((sum(pags(:,s) .* (a(2) -a(1))) - 1) < 1e-3)
    %end

    av_alpha = 0;
    for i = 1:numel(a)
        av_alpha = av_alpha + da * qa(i) * a(i);
    end


end