%% Test with outer loop
base_args = {'beta',0.8, 's0', 4, 'A0',1, 'debug', 4, ...
             'jacobian_fun', @jacobian_noEL_2, 'nk', 30, ...
              'm',15, 'om_pts', 200, 'Ns', 3};
         
% --- outer-loop configuration: ---
%outerloop_pname = 'om_pts';             % parameter name to pass to equilSol
%outerloop_pvals = linspace(200,2400,3); % numeric vector of values
%outerloop_legend_name = 'N_{\omega}';  % string used in legend (LaTeX)

% conclusion: no impact ever since I did trapz -> mean in the P residuals


%outerloop_pname = 'm';             % parameter name to pass to equilSol
%outerloop_pvals = linspace(15,31,3); % numeric vector of values
%outerloop_legend_name = 'm';  % string used in legend (LaTeX)

% an epsilon independent shift by very very little so these are all
% converged in m

%outerloop_pname = 'Ns';             % parameter name to pass to equilSol
%outerloop_pvals = [3]; % numeric vector of values
%outerloop_legend_name = 'N_S';  % string used in legend (LaTeX)

outerloop_pname = 'Nb';             % parameter name to pass to equilSol
outerloop_pvals = [1]; % numeric vector of values
outerloop_legend_name = 'N_B';  % string used in legend (LaTeX)

% no change at all, these are converged in N_S, actually in \Delta_1 we can
% see some change from N_s = 1 to N_S = 2 because now we capture the
% triangularity contribution, but then the squareness does not contribute.

n_om = numel(outerloop_pvals);

% eps grid
nsim = 20;
values_eps = logspace(log10(0.001),log10(0.4),nsim);

% storage: cell arrays per om_pts
delta1_num_all = cell(n_om,1);
delta1_ana_all = cell(n_om,1);
delta_num_all = cell(n_om,1);
delta0_ana_all = cell(n_om,1);
L2error_all = zeros(n_om, nsim);
L20error_all = zeros(n_om, nsim);
L2t2error_all = zeros(n_om, nsim);   % new
L2S2error_all = zeros(n_om, nsim);

% markers to vary with om_pts
marker_list = {'o','s','^','d','v','>','<' ,'p'}; % extend if needed

% loop over om_pts
for k = 1:n_om
    % parameters
    val = outerloop_pvals(k);
    args = [base_args, {outerloop_pname, val}];
    [L, LX] = equilSol(args{:});
    LX.eps_val = 0.4;

    LX.Sbc(1) = -0.8;
    LX.S2bc = -0.8;
    
    if L.P.Ns >1
        LX.S3bc = 0.2;
        LX.Sbc(2) = 0.2;
    end
    LY = equilY(L, LX);

    LX.x = LY.x;
    L.P.hot_restart = 'true';
    fprintf('--- Starting %s = %d (run %d of %d) ---\n', outerloop_pname, val, k, n_om);

    % initial equilibration using previous LX.x when available
    % first run uses LX.x from above

    % preallocate
    [delta0_ana, delta_num, delta1_ana,...
        delta1_num,t2_num,S2_num] = deal(zeros(nsim, L.Nq));


    % initial run (use last eps = 0.4 initial state already present)
    % run with the current LX.eps_val (already set to 0.4)
    LY = equilY(L, LX);
    LX.x = LY.x; % warm start

    % store first
    delta0_ana(1,:) = LY.delta_ana;
    delta_num(1,:) = LY.delta(2:end);
    delta1_ana(1,:) = interp1(LY.r_fine, LY.delta1_ana, L.r_q, 'spline');
    delta1_num(1,:) = (LY.delta(2:end) - LY.delta_ana) / LX.eps_val;
    L2error_all(k,1) = norm(delta1_ana(1,:) - delta1_num(1,:));
    L20error_all(k,1) = norm(delta0_ana(1,:) - delta_num(1,:));
    L2t2error_all(k,1) = norm(LY.t2_ana - LY.t2);
    L2S2error_all(k,1) = norm(LY.S2_ana - LY.S(2:end,:,1));
    
    t2_num(1,:) = LY.t2;
    S2_num(1,:) = squeeze(LY.S(2:end,:,1));

    % sweep eps decreasing (as in original script)
    for ii = 1:nsim-1
        LX.eps_val = values_eps(end-ii);
        LY = equilY(L, LX);
        fprintf('%s %d, Run %d, eps %.1e, |res| = %.4e, bp_LIU = %.1e, li_LIU = %.1e \n', ...
                  outerloop_pname, val, ii, LX.eps_val, LY.res_norms(end), LY.bp_liu, LY.li_liu);
        LX.x = LY.x;

        delta0_ana(ii+1,:) = LY.delta_ana;
        delta_num(ii+1,:) = LY.delta(2:end);
        delta1_ana(ii+1,:) = interp1(LY.r_fine, LY.delta1_ana, L.r_q, 'spline');
        delta1_num(ii+1,:) = (LY.delta(2:end) - LY.delta_ana) / LX.eps_val;
        L2error_all(k,ii+1) = norm(delta1_ana(ii+1,:) - delta1_num(ii+1,:));
        L20error_all(k,ii+1) = norm(delta0_ana(ii+1,:) - delta_num(ii+1,:));
        L2t2error_all(k,ii+1) = norm(LY.t2_ana - LY.t2);
        L2S2error_all(k,ii+1) = norm(LY.S2_ana - LY.S(2:end,:,1));
        
        % for profile plots
        t2_num(ii+1,:) = LY.t2;
        S2_num(ii+1,:) = squeeze(LY.S(2:end,:,1));
        
    end

    % save per-om results
    delta1_num_all{k} = delta1_num;
    delta1_ana_all{k} = delta1_ana;
    delta_num_all{k} = delta_num;
    delta0_ana_all{k} = delta0_ana;
end

% reverse eps order to match original plotting if needed
values_eps = flip(values_eps);

% prepare colormap based on eps (ascending)
eps_sorted = sort(values_eps, 'ascend');
n_vals = numel(eps_sorted);
cmap = jet(n_vals);
log_eps = log10(eps_sorted);
log_eps_norm = (log_eps - min(log_eps)) / (max(log_eps) - min(log_eps));
col_inds = round(1 + log_eps_norm * (size(cmap,1) - 1));

figure;
n_panels = 4;
errors = {L2error_all, L20error_all, L2t2error_all, L2S2error_all};
titles = {'$L^2_{\Delta_1}$ error', '$L^{2}_{\Delta_0}$ error', '$L^2_{t_2}$ error', '$L^2_{S_2}$ error'};
skip_counts = [8, 5, 5, 6];
for ip = 1:n_panels
    subplot(n_panels,1,ip);
    plot_error_panel(errors{ip}, values_eps, outerloop_pvals, cmap, col_inds, marker_list, titles{ip}, outerloop_legend_name, skip_counts(ip));
end

% rerun at high eps for bp, li kappa, delta
val = outerloop_pvals(k);
args = [base_args, {outerloop_pname, val}];
[L, LX] = equilSol(args{:});
LX.eps_val = 0.4;

LX.Sbc(1) = -0.8;
LX.S2bc = -0.8;

if L.P.Ns >1
    LX.S3bc = 0.2;
    LX.Sbc(2) = 0.2;
end
LY = equilY(L, LX);

suptitle_str = sprintf(['Parameters: $\\beta = %.1f$, $s_0 = %.1f$, $A_0 = %.1f$, ', ...
                        '$S_2^{BC} = %.1f$\n LIUQE parameters: ($\\epsilon = 0.4$): $\\beta_{p,LIU} = %.1f$, $l_{i,LIU} = %.1f$, $\\kappa = %.1f$, $\\delta = %.1f$'], ...
                        L.P.beta, L.P.s0, L.P.A0, LX.S2bc, LY.bp_liu, LY.li_liu, LY.kappa(end), LY.deltatrig(end));
sgtitle(suptitle_str, 'Interpreter', 'latex');



% since I reran, now I can compare with LO and NLO

% or better, all at the same time
figure;
subplot(4,1,1);
plot_profile_panel(L.r_q, LY.delta_ana, delta_num, cmap, col_inds, '$\Delta_0$ profiles', '$\hat\Delta_0$');
hold on; 
h = plot(L.r_q, LY.delta_ana + 0.4 * interp1(LY.r_fine, LY.delta1_ana, L.r_q, 'spline')...
    , 'LineWidth',2,'Color', [1 0.35 0.5]);
legend(h, {'$\Delta_0 + \epsilon \Delta_1$'}, 'Interpreter','latex','FontSize',12,'Box','off', 'Location', 'Northwest')
subplot(4,1,2);
plot_profile_panel(L.r_q, interp1(LY.r_fine, LY.delta1_ana, L.r_q, 'spline'), delta1_num, cmap, col_inds, '$\Delta_1$ profiles', '$\hat\Delta_1$')
subplot(4,1,3);
plot_profile_panel(L.r_q, LY.t2_ana, t2_num, cmap, col_inds, '$t_2$ profiles', '$t_2$')
subplot(4,1,4);
plot_profile_panel(L.r_q, LY.S2_ana, S2_num, cmap, col_inds, '$\hat S_2$ profiles', '$\hat S_2$')
hold on; 
h = plot(L.r_q, LY.S2_ana + 0.4 * interp1(LY.r_fine, LY.S2_1_fine, L.r_q, 'spline'), ...
          'LineWidth',2,'Color',[1 0.35 0.5]);
legend(h, {'$S_{2,0} + \epsilon S_{2,1}$'}, 'Interpreter','latex','FontSize',12,'Box','off')
sgtitle(suptitle_str, 'Interpreter', 'latex');
