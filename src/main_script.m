%%% 
%%% Optimal metabolic strategies for microbial growth    
%%% in stationary random environment
%%% Authors: Anna Paola Muntoni and Andrea De Martino
%%%


clear all
close all

%%%% Plot properties
set(0, 'DefaultAxesXgrid', 'off');
set(0, 'DefaultAxesFontName', 'Helvetica');
set(0, 'DefaultAxesFontSize', 20);
set(0, 'DefaultColorbarFontSize', 20);
set(0, 'DefaultColorbarFontName', 'Helvetica');
set(0, 'DefaultLineLineWidth', 3.0);
myfontsize = 30;
%mycolormap = "YlOrBr";
%mycolormap2 = "RdYlBu";
save_fig = 1;
%%%%


%%% number of values taken by the stress, the degree of freedom x and
%%% the growth rate \mu(x,s)
nrs = 500;
nra = 500;
nrl = 200;
nrq = 100;
nre = 100;

%%% 
phi = 0.48;
q0 = 0.5;
x0 = 0.05;
w = 0.169;
%%% maintenance
m = 0.0;
mc = phi/x0;


% nutrients requirement for r: respiration f: fermentation
qr = 1;
qf = 10;
% expenditure for r and f
xr = 1.;
xf = 0.1;
% sc is the constant appearing in hat(x) (s). 

sc = (xr-xf)./(qf-qr);
disp('s_c ='); disp(sc);

amin = 0.0001;
amax = 1 -amin;
a = linspace(amin,amax,nra);
da = (amax-amin)./numel(a);
nu = 3./2;
en = 1./(nu-1);

%%% define the distribution of the stress level
%%% possible values here uniform, bimodal, exponential and power law
%%% distribution
type_stress = 'exp';
switch type_stress
    case 'uniform'
        beta_all = logspace(-2, 2, 11);
        s0 = 0.6;
        [smin, smax, ss, ps, extra] = uniform_environment(s0, nrs);
    case 'exp'
        beta_all = logspace(-1, 2, 9);
        s0 = 0.5;
        delta = 1.0;
        [smin, smax, ss, ps, extra] = exp_environment(s0, delta, nrs);
    case 'power_law'
        beta_all = logspace(-2, 2, 9);
        s0 = 0.1;
        delta = 5;
        [smin, smax, ss, ps, extra] = power_law_environment(s0, delta, nrs);
    case 'two_states'
        beta_all = logspace(-1, 2, 11);
        s0 = 0.05;
        lb1 = 0.0;
        ub1 = 0.1;
        s1 = 0.8;
        lb2 = 0.7;
        ub2 = 0.9;
        sigma = 1e-2; 
        rho = 0.95;
        [smin, smax, ss, ps, extra] = two_states_environment(s0, s1, rho, delta, nrs, sigma, lb1, ub1, lb2, ub2);
    case 'two_states_broad'
        s0 = 0.05;
        s1 = 0.8;
        sigma = 1e-2;
        rho = 0.95;
        beta_all = logspace(-1, 2, 11);
        unip = 0.8;
        [smin, smax, ss, ps, extra] = two_states_broad_environment(s0, s1, rho, nrs, sigma, unip);
end

extra = extra + '_m_' + string(m);

%%% values for \beta
beta_all(end+1) = 0;
beta_all(end+1) = 300;
beta_all = sort(beta_all);
%%%

%%% colormaps
%cmap_lambda = flip(brewermap(numel(beta_all), "RdYlBu"));
%cmap_qa = flip(brewermap(numel(beta_all), "RdYlBu"));
cmap_lambda = turbo(length(beta_all));
cmap_qa = turbo(length(beta_all));
alphac = 0.7;
%%%
%%
%%% compute the map \mu(x,s)
for i = 1:numel(ss)-1
    s(i) = 0.5*(ss(i+1) + ss(i));
    % best value of x given s, hat(x)
    ahatth(i) = power(s(i),en)./(power(s(i),en)+power(sc,en));
    ds(i) = ss(i+1) - ss(i);
    for j=1:numel(a)
        a(j) = amin+(j-1)*da;
        q(j) = (qr+power((1-a(j)),nu)*(qf-qr)); % nutrient req for hat(x)
        eps(j) = (xf+power(a(j),nu)*(xr-xf)); % expenditure for hat(x)
        % this is \mu(x, s)
        f(j,i) = (phi-(s(i)*q0+x0)*m)/(w+s(i)*q(j)+eps(j));
        if f(j,i) <= 0
            f(j,i) = NaN;
        end
    end
end
%%
fmin = min(min(f));
fmax = max(max(f));
disp('\mu min ='); disp(fmin);
disp('\mu max ='); disp(fmax);
lam = linspace(fmin,fmax,nrl);
dlam = (fmax-fmin)./numel(lam);

ferm = zeros(numel(beta_all),1);
av_lambda = zeros(numel(beta_all),1);
for b = 1:numel(beta_all)
    beta = beta_all(b);
    %%% get the optimal distribution p*(x) and p*(x | s)
    [qa, pags, av_alpha] = get_optimal(a, ss, da, ds, ps, f, beta);
    %%% get the distribution for the growth rate, mean value for x, optimal
    %%% I* and \mu*
    [plam, Istar, fstar, aav] = get_mu_dist(a, ss, da, ds, f, pags, qa, lam, dlam, ps);
    %%% get the distribution for q and eps
    [pq, peps, qlin, epslin] = get_q_eps_dist(a, ss, da, ds, q, eps, nrq, nre, ps, pags);
    
    %%% save all results to plot them later
    I_all(b) = Istar;
    f_all(b) = fstar;
    av_lambda(b) = fstar;
    ferm(b) = 1.0 - av_alpha;
    fmax_all(b) = fmax;
    aav_all{b} = aav;
    ahatth_all{b} = ahatth;

    lam_b(b,:) = lam;
    plam_b(b,:) = plam;
    
    qlin_b(b,:) = qlin;
    pq_b(b,:) = pq;

    epslin_b(b,:) = epslin;
    peps_b(b,:) = peps;
    
    a_b(b,:) = a;
    qa_b(b,:) = qa;
    pags_all{b} = pags;

end


%% Plot 
figure('units','normalized','outerposition',[0 0 1 1])

for b = 1:length(beta_all)
    beta = beta_all(b);
    area(lam_b(b,:), plam_b(b,:), 'FaceColor', cmap_lambda(b,:) , ...
        'FaceAlpha', alphac, 'LineStyle', '-', 'EdgeColor', cmap_lambda(b,:));
    hold on
end
xlim([fmin, fmax])
xlabel("$\mu$",'fontsize',myfontsize, 'Interpreter','Latex');
ylabel("$p^{\star}(\mu | \beta)$",'fontsize',myfontsize, 'Interpreter', 'Latex');
hold on;
colormap(cmap_lambda);

h = colorbar('Ticks', linspace(0,1,numel(round(beta_all(1:2:end),2))), ...
    'TickLabels', round(beta_all(1:2:end),2), 'Location', 'Southoutside');
ylabel(h, '\beta', 'fontsize', myfontsize);

pbaspect([1.5 1 1])
if save_fig
    filename = type_stress + extra + "_p_mu_all.pdf";
    print(filename, '-dpdf', '-bestfit');
end
pause(1)
close()

figure('units','normalized','outerposition',[0 0 1 1])
for b = 1:length(beta_all)
    beta = beta_all(b);
    area(a_b(b,:), qa_b(b,:), 'FaceColor', cmap_qa(b,:) , ...
        'FaceAlpha', alphac, 'LineStyle', '-', 'EdgeColor', cmap_qa(b,:));
    hold on
end
colormap(cmap_qa);
h = colorbar('Ticks', linspace(0,1,numel(round(beta_all(1:2:end),2))), ...
    'TickLabels', round(beta_all(1:2:end),2), 'Location', 'Southoutside');
ylabel(h, '\beta', 'fontsize', myfontsize);
xlabel("$x$",'fontsize', myfontsize, 'Interpreter', 'Latex');
ylabel("$p^{\star}(x)$",'fontsize', myfontsize, 'Interpreter', 'Latex');
hold on;
pbaspect([1.5 1 1])
if save_fig
    filename = type_stress + extra + "_p_x_all.pdf";
    print(filename, '-dpdf', '-bestfit');
end
close all

%%
vmin = Inf;
vmax = -Inf;
vmin_log = Inf;
vmax_log = -Inf;

for b = 1:numel(beta_all)
    pags = pags_all{b};
    vmin = min(vmin, min(min(pags(pags > 0))));
    vmax = max(vmax, max(max(pags)));
end

for b = 1:numel(beta_all)
    for i = 1:size(pags,1)
        for j = 1:size(pags,2)
            if isnan(pags_all{b}(i,j))
                continue
            end
            pags_all{b}(i,j) = max(vmin, pags_all{b}(i,j));
        end
    end
end

for b = 1:numel(beta_all)
    pags = pags_all{b};
    vmin_log = min(vmin_log, min(min(log10(pags))));
    vmax_log = max(vmax_log, max(max(log10(pags))));
end

for b = 1:numel(beta_all)
    beta = beta_all(b);
    pags = pags_all{b};
    lam =  lam_b(b,:);
    plam = plam_b(b,:);
    qlin = qlin_b(b,:);
    pq = pq_b(b,:);
    epslin = epslin_b(b,:);
    peps = peps_b(b,:);
    a = a_b(b,:);
    qa = qa_b(b,:);
    aav = aav_all{b};
    figure()

    contourf(log10(ss(1:end-1)), a, log10(pags), ...
        'HandleVisibility','off', 'LineStyle','None', 'LineWidth',1.2);
    %colormap(brewermap([],mycolormap2))
    colormap(turbo)
    caxis([vmin_log vmax_log])
    hold on
    xlabel ("$\log_{10}(s)$",'fontsize', myfontsize, 'Interpreter','latex');
    ylabel ("$x$",'fontsize', myfontsize, 'Interpreter', 'Latex');
    h = colorbar;
    ylabel(h,"$\log p^{\star}(x|s)$",'fontsize',myfontsize, 'Interpreter', 'Latex');
    hold on
    plot(log10(ss(1:end-1)),ahatth, '-.', 'Color','g', 'linewidth', 5.0);
    hold on
    plot(log10(ss(1:end-1)),aav, '--', 'Color', 'm', 'linewidth', 5.0);
    legend( '$\hat{x}$', '$< x >_{p(x|s)}$', ...
        'Interpreter', 'Latex', 'NumColumns', 2, 'Location','southoutside');
    legend boxoff
    pbaspect([1.5 1 1])
    hold off
    if save_fig
        filename = type_stress + extra + ...
            "_p_x_given_s_beta_" + sprintf('%1.2e', beta) + ".pdf";
        print(filename, '-dpdf', '-bestfit')
        pause(1)
    end
    close

    figure()
    area(qlin, pq, 'FaceColor', cmap_lambda(b,:) , ...
        'FaceAlpha', alphac, 'LineStyle', '-', 'EdgeColor', cmap_lambda(b,:));
    hold on
    plot(qlin, pq, 'Color', cmap_lambda(b,:));
    xlabel("$q$",'fontsize',myfontsize, 'Interpreter','Latex');
    ylabel("$p^{\star}(q| \beta)$",'fontsize',myfontsize, 'Interpreter', 'Latex');
    %xlim([qmin, qmax])
    pbaspect([1.5 1 1])

    if save_fig
        filename = type_stress + extra + "_beta_" + ...
            sprintf('%1.2e', beta) + "_p_q.pdf";
        print(filename, '-dpdf', '-bestfit');
        pause(1)
    end
    close

    figure()
    area(epslin, peps, 'FaceColor', cmap_lambda(b,:) , ...
        'FaceAlpha', alphac, 'LineStyle', '-', 'EdgeColor', cmap_lambda(b,:));
    hold on
    plot(epslin, peps, 'Color', cmap_lambda(b,:));
    xlabel("$\epsilon$",'fontsize',myfontsize, 'Interpreter','Latex');
    ylabel("$p^{\star}( \epsilon | \beta)$",'fontsize',myfontsize, 'Interpreter', 'Latex');
    %xlim([epsmin, epsmax])
    pbaspect([1.5 1 1])

    if save_fig
        filename = type_stress + extra + "_beta_" + ...
            sprintf('%1.2e', beta) + "_p_eps.pdf";
        print(filename, '-dpdf', '-bestfit');
        pause(1)
    end
    close

    figure()

    area(lam, plam, 'FaceColor', cmap_lambda(b,:) , ...
        'FaceAlpha', alphac, 'LineStyle', '-', 'EdgeColor', cmap_lambda(b,:));
    hold on
    plot(lam, plam, 'Color', cmap_lambda(b,:));
    xlabel("$\mu$",'fontsize',myfontsize, 'Interpreter','Latex');
    ylabel("$p^{\star}(\mu | \beta)$",'fontsize',myfontsize, 'Interpreter', 'Latex');
    xlim([fmin, fmax])
    pbaspect([1.5 1 1])

    if save_fig
        filename = type_stress + extra + "_beta_" + ...
            sprintf('%1.2e', beta) + "_p_mu.pdf";
        print(filename, '-dpdf', '-bestfit');
        pause(1)
    end
    close

    figure()
    area(a, qa, 'FaceColor', cmap_qa(b,:) , ...
        'FaceAlpha', alphac, 'LineStyle', '-', 'EdgeColor', cmap_qa(b,:));
    hold on
    plot(a, qa, 'Color', cmap_qa(b,:));
    hold on
    xlabel("$x$",'fontsize',myfontsize, 'Interpreter','Latex');
    ylabel("$p^{\star}(x)$",'fontsize',myfontsize, 'Interpreter', 'Latex');
    hold on;
    pbaspect([1.5 1 1])

    if save_fig
        filename = type_stress + extra + "_beta_" + ...
            sprintf('%1.2e', beta) + "_p_x.pdf";
        print(filename, '-dpdf', '-bestfit');
        pause(1)
    end
    close

    
end

figure()
ss_plot = linspace(max(0,smin), smax, 10000);
area(ss, ps, 'FaceAlpha',0.5, 'LineStyle','None');
if strmatch(type_stress, 'power_law')
    semilogx(ss, ps, 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',8.0);
end
xlim([smin smax]);
ylim([0, max(ps) + 1])
xlim([smin smax]);
set(gca, 'FontSize', 35);
ax = gca;
xlabel('$s$','Interpreter', 'Latex', 'FontSize', 40);
ylabel('$p\left(s\right)$', 'Interpreter', 'Latex', 'FontSize', 40);
pbaspect([1.5 1 1])
filename = type_stress + extra + "_environment.pdf";
if save_fig
        print(filename, '-dpdf', '-bestfit');
        pause(1)
end
close

figure()
subplot(1,2,1)
contourf(log10(ss(1:end-1)),a,f,40, 'HandleVisibility','Off', 'LineStyle','None');
%colormap(flip(brewermap([], "Spectral")))
colormap(turbo)
hold on
hold on;
plot(log10(ss(1:end-1)),ahatth,"white", 'LineWidth',8.0); 
set(gca, 'FontSize', 35);

hold on
xlabel ("$\log_{10}\left(s \right) $", 'Interpreter', 'Latex', 'Fontsize', 40);
xlim([log10(smin) log10(smax)]);
ylabel ("$x$", 'Interpreter', 'Latex', 'FontSize', 40);
originalSize1 = get(gca, 'Position');

h = colorbar();
ylabel(h,"$\mu(x,s)$",'Fontsize', 40, 'Interpreter', 'Latex');
caxis([fmin fmax]);

set(gca,'Position', originalSize1)
pbaspect([1.5 1 1])
filename = type_stress + extra + "_mu_map.pdf";
if save_fig
    print(filename, '-dpdf', '-bestfit');
    pause(1)
end
close

