%% Section 4.2: Isotropic Plasmas
% This script can be run as is to reproduce Figs. 3, 4, 5 of the article
% Van Parys et al., 
% "Investigation of finite aspect ratio effects in axisymmetric 
% magnetic equilibria with toroidal rotation and pressure anisotropy"



base_args = {'beta',0.8, 's0', 4, 'debug', 4, ...
            'residuals_fun', @residuals_iso_static, 'equation_of_state',@isotropic, ...
             'jacobian_fun', @jacobian_iso_static, 'nk', 30, ...
              'm',15, 'om_pts', 200, 'Nb', 1, 'mach20',0,'Ns',3};
         
eps_max = 0.4;

[L, LX] = equilSol(base_args{:});
LX.eps_val =eps_max;

LX.Sbc(1) = -0.8;
LX.Sbc(2) = 0.2;

LY = equilY(L, LX);
LX.x = LY.x;L.P.do_shift_NLO = true;
LY = equilY(L, LX);
%% scan_eps_from equilibrium part
% Basically scan_eps_from_equilibrium written out, but for source control
L0 = L; LX0 = LX; LY0 = LY;


nsim = 20;

eps_max = LX0.eps_val;
values_eps = logspace(log10(0.001), log10(eps_max), nsim);  
L2error_all      = zeros(1, nsim);
L20error_all     = zeros(1, nsim);
L2t2error_all    = zeros(1, nsim);
L2S2error_all    = zeros(1, nsim);
L2S2_1_error_all = zeros(1, nsim);
H1_error_all = zeros(1, nsim);
fprintf('--- Single-run sweep, eps in [%g, %g], nsim = %d ---\n', values_eps(1), values_eps(end), nsim);

Nq = L0.Nq;
[delta_num, ...
    delta1_num, t2_num, S2_num, deltap_num, S2_1_num] = deal(zeros(nsim, Nq));

delta0_ana = LY0.delta_ana;
delta1_ana = interp1(LY0.r_fine, LY0.delta1_ana, L0.r_q, 'spline');
S2_0_ana = LY0.S2_ana;
S2_1_ana = interp1(LY0.r_fine, LY0.S2_1_fine, L0.r_q, 'spline');
t2_ana = LY0.t2_ana;
t2p_ana = LY0.t2p_ana;

t2_num(1,:)    = LY0.t2(2:end-1);
delta_num(1,:)  = LY0.delta(2:end-1);
delta1_num(1,:) = (LY0.delta(2:end-1) - delta0_ana) / LX0.eps_val;
S2_num(1,:)  = squeeze(LY0.S(2:end-1,1,1));
S2_1_num(1,:)   = (squeeze(LY0.S(2:end-1,1,1)) - S2_0_ana) / LX0.eps_val;

L2error_all(1)      = norm(delta1_ana - delta1_num(1,:).');
L20error_all(1)     = norm(delta0_ana - delta_num(1,:).');
L2t2error_all(1)    = norm(t2_ana - t2_num(1,:).');
L2S2error_all(1)    = norm(S2_0_ana - S2_num(1,:).');
L2S2_1_error_all(1) = norm(S2_1_ana - S2_1_num(1,:).');
H1_error_all(1) = norm(t2p_ana - t2_num(1,:).');
 
deltap_num(1,:) = LY0.deltap(2:end-1);

% set previous LY for warm-start in next iterations
LY_prev = LY0;
LX = LX0; L = L0; L.P.hot_restart = true; L.P.do_shift_NLO = false;                     
for ii = 1:nsim-1
    eps_val = values_eps(end-ii);    % like original script: decreasing eps
    LX.eps_val = eps_val;
    LX.x = LY_prev.x;                % warm start from previous result
    % run solver
    LY = equilY(L, LX);

    fprintf('Run %d/%d, eps = %.1e, |res| = %.4e, bp = %.1f, li = %.1f\n', ...
            ii+1, nsim, LX.eps_val, LY.res_norms(end), LY.bp, LY.li);

    LY_prev = LY;

    % store results at index ii+1
    delta_num(ii+1,:)  = LY.delta(2:end-1);
    delta1_num(ii+1,:) = (LY.delta(2:end-1) - delta0_ana) / LX.eps_val;
    t2_num(ii+1,:)    = LY.t2(2:end-1);
    S2_num(ii+1,:)  = squeeze(LY.S(2:end-1,1,1));
    deltap_num(ii+1,:) = LY.deltap(2:end-1);
    S2_1_num(ii+1,:) = (squeeze(LY.S(2:end-1,1,1))- S2_0_ana) / LX.eps_val;

    L2error_all(1,ii+1)      = norm(delta1_ana - delta1_num(ii+1,:).');
    L20error_all(1,ii+1)     = norm(delta0_ana - delta_num(ii+1,:).');
    L2t2error_all(1,ii+1)    = norm(t2_ana -  t2_num(ii+1,:).');
    L2S2error_all(1,ii+1)    = norm(S2_0_ana - S2_num(ii+1,:).');
    L2S2_1_error_all(1,ii+1) = norm(S2_1_ana - S2_1_num(ii+1,:).');
    H1_error_all(1,ii+1) = norm(t2p_ana - t2_num(ii+1,:).');
end
% reverse eps order to match original plotting if needed
values_eps = flip(values_eps);

% prepare colormap based on eps (ascending)
eps_sorted = sort(values_eps, 'ascend');


n_vals = numel(eps_sorted);
h = linspace(0,0.75,n_vals)';          % hue sweep
s = 0.8*ones(n_vals,1);                 % saturation = 1
v = 0.92*ones(n_vals,1);                 % value = 1

cmap = hsv2rgb([h s v]);
log_eps = log10(eps_sorted);
log_eps_norm = (log_eps - min(log_eps)) / (max(log_eps) - min(log_eps));
col_inds = round(1 + log_eps_norm * (size(cmap,1) - 1));

%% Convergence plot

figure;
t = tiledlayout(1,2,'TileSpacing','compact','Padding','none');
ax = nexttile(1);

h1 = plot_error_panel_2(L20error_all, values_eps, 0.02, cmap, col_inds, ...
    {'o'}, '$\hat\Delta_0$');
hold on;
h2 = plot_error_panel_2(L2error_all .* values_eps, values_eps, 0.02, cmap, col_inds, ...
    {'square'}, '$\hat\Delta_0+\epsilon\hat\Delta_1$');

legend([h1 h2], 'Interpreter','latex', ...
       'FontSize',14, 'Box','off');

title('$\|\hat\Delta_{\mathrm{num}}-\hat\Delta_{\mathrm{ana}}\|_{L^2}$', ...
      'Interpreter','latex','FontSize',14);
xlabel('$\epsilon$', 'Interpreter', 'latex','FontSize',14)

axis(ax,'tight');
set(ax,'LooseInset',get(ax,'TightInset'));
ylim([1e-4 5.48])
xlim([1e-3 max(values_eps)])
yticks([1e-4, 1e-3 1e-2 1e-1 1e0])
ax = nexttile(2);

h1 = plot_error_panel_2(L2S2error_all, values_eps, 0.02, cmap, col_inds, ...
    {'o'}, '$\hat S_{2,0}$');
hold on;
h2 = plot_error_panel_2(L2S2_1_error_all .* values_eps, values_eps, 0.02, cmap, col_inds, ...
    {'square'},  '$\hat S_{2,0}+\epsilon\hat S_{2,1}$');
hold on;
plot([0.155313 0.155313], [1e-4 1e2], '--', 'Color', 'r')
legend([h1 h2], 'Interpreter','latex', ...
       'FontSize',14, 'Box','off');


title('$\|\hat S_{2,\mathrm{num}}-\hat S_{2,\mathrm{ana}}\|_{L^2}$', ...
      'Interpreter','latex','FontSize',14);

xlabel('$\epsilon$', 'Interpreter', 'latex','FontSize',14)

axis(ax,'tight');
set(ax,'LooseInset',get(ax,'TightInset'));
ylim([6.8e-6 4.67])
xlim([1e-3 max(values_eps)])
yticks([1e-5 1e-4, 1e-3 1e-2 1e-1 1e0])


%% Profiles plot
figure;
t = tiledlayout(1,2,'TileSpacing','compact','Padding','none');
ax = nexttile(1);
plot_profile_panel(L0.r_q, delta0_ana, delta_num, cmap, col_inds, '$\Delta_0$ profiles', '$\hat\Delta_0$');
delete(ax.Title)
hold on;
h = plot(L0.r_q, delta0_ana + eps_max * delta1_ana...
    ,'--', 'LineWidth',2,'Color', [1 0.35 0.5]);
legend(h, {'$\hat\Delta_0 + \epsilon \hat\Delta_1$ $(\epsilon=0.4)$'}, 'Interpreter','latex','FontSize',12,'Box','off', 'Location', 'Northwest')
axis(ax,'tight'); set(ax,'LooseInset',get(ax,'TightInset'));
ax = nexttile(2);
plot_profile_panel(L0.r_q,delta1_ana, delta1_num, cmap, col_inds, '$\hat \Delta_1$ profiles', '$\hat \Delta_{1}$')
delete(ax.Title)
legend( {'$\textrm{analytical}$'}, 'Interpreter','latex','FontSize',12,'Box','off', 'Location', 'Southwest')
axis(ax,'tight'); set(ax,'LooseInset',get(ax,'TightInset'));

% elongation
figure;
t = tiledlayout(1,2,'TileSpacing','compact','Padding','none');
ax = nexttile(1);
plot_profile_panel(L0.r_q, S2_0_ana, S2_num, cmap, col_inds, '$\hat S_2$ profiles', '$\hat S_2$');
delete(ax.Title)
hold on;
h = plot(L0.r_q, S2_0_ana + eps_max * S2_1_ana, ...
              '--','LineWidth',2,'Color',[1 0.35 0.5]);
    legend(h, {'$\hat S_{2,0} + \epsilon \hat S_{2,1}$ $(\epsilon=0.4)$'}, 'Interpreter','latex','FontSize',12,'Box','off', 'Location', 'Southwest')
axis(ax,'tight'); set(ax,'LooseInset',get(ax,'TightInset'));
ax = nexttile(2);
plot_profile_panel(L0.r_q, S2_1_ana, S2_1_num, cmap, col_inds, '$\hat S_2$ profiles', '$\hat S_{2,1}$');
delete(ax.Title)
legend( {'$\textrm{analytical}$'}, 'Interpreter','latex','FontSize',12,'Box','off', 'Location', 'Southwest')
axis(ax,'tight'); set(ax,'LooseInset',get(ax,'TightInset'));