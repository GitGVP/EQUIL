%% Variational epsilon convergence against equil_ana
% Reproduces the layout and colour ordering of Figs. 3--5 of
% Van Parys et al. (2026) with the variational solver.

nsim = 20;
eps_max = 0.4;
values_eps = flip(logspace(log10(1e-3),log10(eps_max),nsim));

% These settings are intentionally stricter than the routine defaults.
% Angular moments are formed with mean(...) in the variational assembly.
[L,LX] = equilVariationalSol('beta',0.8,'s0',4, ...
    'm',8,'nq',8,'spline_p',4,'om_pts',128,'Ns',3, ...
    'NLtol',3e-13,'nk',45,'debug',0);
LX.eps_val = eps_max;
LX.Sbc = [-0.8;0.2;0];
LY = equilVariationalY(L,LX);
assert(LY.isconverged);

% equil_ana and its expensive NLO elongation solve are independent of
% epsilon here.  Evaluate both exactly once, outside the continuation.
Lana = L;
Lana.P.do_shift_NLO = true;
analytic = equil_ana(Lana,LX);
delta0 = analytic.delta_ana;
delta1 = interp1(analytic.r_fine,analytic.delta1_ana,L.r_q,'pchip');
S20 = analytic.S2_ana;
S21 = interp1(analytic.r_fine,analytic.S2_1_fine,L.r_q,'pchip');
t20 = analytic.t2_ana;

nr = L.Nq;
delta_num = zeros(nsim,nr);
delta1_num = zeros(nsim,nr);
S2_num = zeros(nsim,nr);
S2_1_num = zeros(nsim,nr);
t2_num = zeros(nsim,nr);
errors = zeros(nsim,5);
residuals = zeros(nsim,1);
% The article scripts use the Euclidean norm of the sampled radial values
% while labelling it L2.  Retain that convention here to reproduce their
% vertical scale; quadrature weighting changes the prefactor, not orders.
wnorm = @(v)norm(v(:));

for k = 1:nsim
    if k > 1
        LX.eps_val = values_eps(k);
        LX.x = LY.x;
        LX.local_B_guess = LY.local_B_quadrature;
        LY = equilVariationalY(L,LX);
    end
    assert(LY.isconverged, ...
        'Variational solve failed at epsilon=%g.',values_eps(k));

    delta_num(k,:) = LY.delta(2:end-1).';
    S2_num(k,:) = LY.S(2:end-1,1,1).';
    t2_num(k,:) = LY.t2(2:end-1).';
    delta1_num(k,:) = (delta_num(k,:)-delta0.')/values_eps(k);
    S2_1_num(k,:) = (S2_num(k,:)-S20.')/values_eps(k);
    errors(k,:) = [ ...
        wnorm(delta_num(k,:).'-delta0), ...
        wnorm(delta_num(k,:).'-delta0-values_eps(k)*delta1), ...
        wnorm(S2_num(k,:).'-S20), ...
        wnorm(S2_num(k,:).'-S20-values_eps(k)*S21), ...
        wnorm(t2_num(k,:).'-t20)];
    residuals(k) = norm(LY.residual);
    fprintf('epsilon %.3e: |R|=%.3e\n',values_eps(k),residuals(k));
end

% Paper fits use the asymptotic interval below epsilon_1=0.02.
fit_mask = values_eps <= 0.02;
orders = fit_orders(values_eps,errors,fit_mask);
fprintf(['Analytical-reference orders: Delta LO %.3f, Delta NLO %.3f, ', ...
    'S2 LO %.3f, S2 NLO %.3f, t2 LO %.3f\n'],orders);

% Independent numerical-reference check.  A quadratic fit through the four
% smallest-epsilon solutions estimates the numerical epsilon->0 LO and NLO
% profiles without using equil_ana.
nref = 4;
ref_rows = nsim-nref+1:nsim;
V = [ones(nref,1),values_eps(ref_rows).',values_eps(ref_rows).'.^2];
delta_coeff = V\delta_num(ref_rows,:);
S2_coeff = V\S2_num(ref_rows,:);
delta0_numref = delta_coeff(1,:).';
delta1_numref = delta_coeff(2,:).';
S20_numref = S2_coeff(1,:).';
S21_numref = S2_coeff(2,:).';
errors_numref = zeros(nsim,4);
for k = 1:nsim
    epsilon = values_eps(k);
    errors_numref(k,:) = [ ...
        wnorm(delta_num(k,:).'-delta0_numref), ...
        wnorm(delta_num(k,:).'-delta0_numref-epsilon*delta1_numref), ...
        wnorm(S2_num(k,:).'-S20_numref), ...
        wnorm(S2_num(k,:).'-S20_numref-epsilon*S21_numref)];
end
numref_fit = fit_mask;
numref_fit(ref_rows) = false;
orders_numref = fit_orders(values_eps,errors_numref,numref_fit);
fprintf(['Numerical-reference orders: Delta LO %.3f, Delta NLO %.3f, ', ...
    'S2 LO %.3f, S2 NLO %.3f\n'],orders_numref);
fprintf('Relative NLO-target difference: Delta %.3e, S2 %.3e\n', ...
    wnorm(delta1_numref-delta1)/max(wnorm(delta1),eps), ...
    wnorm(S21_numref-S21)/max(wnorm(S21),eps));

% Small epsilon is purple; large epsilon is red, as in the article.
log_fraction = (log10(values_eps)-log10(min(values_eps))) ...
    /(log10(max(values_eps))-log10(min(values_eps)));
colors = hsv2rgb([0.75*(1-log_fraction(:)), ...
                  0.8*ones(nsim,1),0.92*ones(nsim,1)]);

%% Fig. 3 counterpart: LO and LO+NLO errors on the same panels
figure;
tiledlayout(1,2,'TileSpacing','compact','Padding','none');
nexttile;
plot_convergence_panel(values_eps,errors(:,1),errors(:,2),colors, ...
    fit_mask,orders(1:2), ...
    '$\|\hat\Delta_{\rm num}-\hat\Delta_{\rm ana}\|_{L^2}$',[]);
nexttile;
plot_convergence_panel(values_eps,errors(:,3),errors(:,4),colors, ...
    fit_mask,orders(3:4), ...
    '$\|\hat S_{2,{\rm num}}-\hat S_{2,{\rm ana}}\|_{L^2}$',0.155313);

%% Fig. 4 counterpart: shift and its NLO profile
figure;
tiledlayout(1,2,'TileSpacing','compact','Padding','none');
nexttile;
plot_profile_panel_article(L.r_q,delta0,delta_num,colors,'$\hat\Delta_0$');
h = plot(L.r_q,delta0+eps_max*delta1,'--','LineWidth',2, ...
    'Color',[1,0.35,0.5]);
legend(h,{'$\hat\Delta_0+\epsilon\hat\Delta_1$ $(\epsilon=0.4)$'}, ...
    'Interpreter','latex','FontSize',12,'Box','off','Location','northwest');
nexttile;
plot_profile_panel_article(L.r_q,delta1,delta1_num,colors,'$\hat\Delta_1$');
hnum = plot(L.r_q,delta1_numref,':','LineWidth',2.2,'Color',[0.3,0.3,0.3]);
hana = plot(nan,nan,'k--','LineWidth',3.5);
legend([hana,hnum],{'analytical','numerical $\epsilon\to0$'}, ...
    'Interpreter','latex','FontSize',12,'Box','off','Location','southwest');

%% Fig. 5 counterpart: elongation and its NLO profile
figure;
tiledlayout(1,2,'TileSpacing','compact','Padding','none');
nexttile;
plot_profile_panel_article(L.r_q,S20,S2_num,colors,'$\hat S_2$');
h = plot(L.r_q,S20+eps_max*S21,'--','LineWidth',2, ...
    'Color',[1,0.35,0.5]);
legend(h,{'$\hat S_{2,0}+\epsilon\hat S_{2,1}$ $(\epsilon=0.4)$'}, ...
    'Interpreter','latex','FontSize',12,'Box','off','Location','southwest');
nexttile;
plot_profile_panel_article(L.r_q,S21,S2_1_num,colors,'$\hat S_{2,1}$');
hnum = plot(L.r_q,S21_numref,':','LineWidth',2.2,'Color',[0.3,0.3,0.3]);
hana = plot(nan,nan,'k--','LineWidth',3.5);
legend([hana,hnum],{'analytical','numerical $\epsilon\to0$'}, ...
    'Interpreter','latex','FontSize',12,'Box','off','Location','southwest');


function orders = fit_orders(epsilon,errors,mask)
    orders = zeros(1,size(errors,2));
    for j = 1:size(errors,2)
        p = polyfit(log(epsilon(mask)),log(errors(mask,j).'),1);
        orders(j) = p(1);
    end
end

function plot_convergence_panel(epsilon,error0,error1,colors,mask,orders,title_text,epsilon2)
    hold on; grid on
    for k = 1:numel(epsilon)
        loglog(epsilon(k),error0(k),'o','LineStyle','none', ...
            'MarkerSize',6,'MarkerFaceColor',colors(k,:), ...
            'MarkerEdgeColor','k');
        loglog(epsilon(k),error1(k),'square','LineStyle','none', ...
            'MarkerSize',6,'MarkerFaceColor',colors(k,:), ...
            'MarkerEdgeColor','k');
    end
    xf = logspace(log10(min(epsilon(mask))),log10(max(epsilon(mask))),100);
    p0 = polyfit(log(epsilon(mask)),log(error0(mask)),1);
    p1 = polyfit(log(epsilon(mask)),log(error1(mask)),1);
    loglog(xf,exp(polyval(p0,log(xf))),'k-','LineWidth',1.5);
    loglog(xf,exp(polyval(p1,log(xf))),'k-','LineWidth',1.5);
    xline(0.02,'k--');
    if ~isempty(epsilon2)
        xline(epsilon2,'r--');
    end
    h0 = loglog(nan,nan,'ok','MarkerFaceColor','k');
    h1 = loglog(nan,nan,'sk','MarkerFaceColor','k');
    legend([h0,h1], ...
        {sprintf('$O(\\epsilon^{%.2f})$',orders(1)), ...
         sprintf('$O(\\epsilon^{%.2f})$',orders(2))}, ...
        'Interpreter','latex','FontSize',12,'Box','off','Location','southeast');
    xlabel('$\epsilon$','Interpreter','latex','FontSize',14);
    title(title_text,'Interpreter','latex','FontSize',14);
    xlim([min(epsilon),max(epsilon)]);
    set(gca,'XScale','log','YScale','log','LooseInset',get(gca,'TightInset'));
end

function plot_profile_panel_article(r,analytical,numerical,colors,ylabel_text)
    hold on
    plot(r,analytical,'k--','LineWidth',3.5);
    for k = 1:size(numerical,1)
        plot(r,numerical(k,:),'Color',colors(k,:));
    end
    grid on
    xlabel('$\hat r$','Interpreter','latex','FontSize',12);
    ylabel(ylabel_text,'Interpreter','latex','FontSize',12);
    axis tight
    set(gca,'LooseInset',get(gca,'TightInset'));
end
