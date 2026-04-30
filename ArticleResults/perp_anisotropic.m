%% Section 4.4: Strong perpendicular anisotropic discharges
% This script can be run as is to reproduce Figs. 8, 9, 10 of the article
% Van Parys et al., 
% "Investigation of finite aspect ratio effects in axisymmetric 
% magnetic equilibria with toroidal rotation and pressure anisotropy"


% the MEQ code is needed for these results
addpath('/home/vanparys/Documents/PhD/Codes/meq.main')


% 1) Get target (q,p)
[L, LX] = fbt('ana',1,0, 'iterq', 50, 'nr', 66, 'nz', 64, ...
              'ifield', true);
bp_target = 2;
LX.bp  = bp_target;
LX.bpD = bp_target;
LX.Wk  = bp_target * LX.Wk;
LY = fbtt(L, LX);
[Le, LXe] = equilSol('beta', 1.1, 's0', 5, 'debug', 4, ... 
                         'residuals_fun',@residuals_noRepl,...
                         'jacobian_fun', @jacobian_noRepl, 'Ns', 3,'Nb',10);
eps_val = 0.37225;
LXe.eps_val = eps_val;
LXe.Sbc(1) = 0; LXe.Sbc(2) = 0;
LYe = equilY(Le, LXe);

LXe_old = LXe;

niter = 10;
for iter = 1:niter
    q_interp = interp1(L.pQ.^2, 1./LY.iqQ, LYe.psiN(2:end-1), 'pchip');  % 150x1
    
    order_max = 10;
    orders = 2:2:order_max; 
    
    K = numel(orders);
    X = zeros(numel(Le.r_q), K);
    for k = 1:K
        X(:,k) = Le.r_q.^(orders(k));
    end
    y = q_interp - 1;                % 150x1
    
    coeffs = X \ y;                   % least-squares
    
    % make row form of coeffs for elementwise operations
    c_row = coeffs(:).';   % 1 x K
    ord_row = orders;      % 1 x K
    
    % q_fit and qp_fit preserve the shape of input rr (1xN or Nx1)
    % uses bsxfun for compatibility with older MATLAB versions
    q_fit = @(rr) reshape( 1 + ...
        sum( bsxfun(@times, bsxfun(@power, rr(:), ord_row), c_row ), 2 ), ...
        size(rr) );
    
    qp_fit = @(rr) reshape( ...
        sum( bsxfun(@times, bsxfun(@power, rr(:), ord_row - 1), ord_row .* c_row ), 2 ), ...
        size(rr) );
    
    q_fitted = q_fit(Le.r_q);             % 150x1, fitted values on Le.r_q
    
    
    
    beta_interp = 1/(eps_val^2  / ( mu0 / (LY.TQ(1)/LY.rA)^2)) * interp1(L.pQ.^2, LY.PQ, LYe.psiN(2:end-1), 'pchip');  % 150x1
    
    order_max = 6;
    orders = 0:2:order_max; 
    
    K = numel(orders);
    X = zeros(numel(Le.r_q), K);
    for k = 1:K
        X(:,k) = Le.r_q.^(orders(k));
    end
    y = beta_interp;                % 150x1
    
    coeffs = X \ y;                   % least-squares
    
    % make row form of coeffs for elementwise operations
    c_row = coeffs(:).';   % 1 x K
    ord_row = orders;      % 1 x K
    
    % q_fit and qp_fit preserve the shape of input rr (1xN or Nx1)
    % uses bsxfun for compatibility with older MATLAB versions
    beta_fit = @(rr) reshape(sum( bsxfun(@times, bsxfun(@power, rr(:), ord_row), c_row ), 2 ), ...
        size(rr) );
    
    betap_fit = @(rr) reshape(sum( bsxfun(@times, bsxfun(@power, rr(:), ord_row - 1), ord_row .* c_row ), 2 ), ...
        size(rr) );

    % blending parameters / schedule
    alpha_q = 1;                             % apply q fully from the start
    % ramp beta in over a few iterations: 0 on iter 1, then 1/3, 2/3, 1...
    alpha_beta = min(1, max(0, (iter-1)/3));

    if isfield(LXe_old,'qfun') && ~isempty(LXe_old.qfun)
        prev_qfun = LXe_old.qfun;
        prev_qpfun = LXe_old.qpfun;
    else
        % fallback: use current fitted values as "previous" (safe for first time)
        prev_qfun = q_fit;
        prev_qvec = q_fitted;
        prev_qpfun = qp_fit;
        prev_qpvec = qp_fit(Le.r_q);
    end

    % blended q functions (apply full q immediately — alpha_q = 1)
    LXe.qfun  = @(rr) (1-alpha_q)*prev_qfun(rr) + alpha_q*q_fit(rr);
    LXe.qpfun = @(rr) (1-alpha_q)*prev_qpfun(rr) + alpha_q*qp_fit(rr);


    % blended beta: do not overwrite on iter==1 (alpha_beta==0)
    if alpha_beta > 0
        % get previous beta vector (if it exists)
        if isfield(LXe_old,'kinetic_profiles') && isfield(LXe_old.kinetic_profiles,'beta') ...
                && ~isempty(LXe_old.kinetic_profiles.beta)
            % evaluate previous beta on Le.r_q
            prev_beta_vec = LXe_old.kinetic_profiles.beta(Le.r_q);
        else
            prev_beta_vec = zeros(size(Le.r_q)); % conservative default
        end

        new_beta_vec = beta_fit(Le.r_q);
        blended_beta_vec = (1-alpha_beta)*prev_beta_vec + alpha_beta*new_beta_vec;

        % store blended vector and function (use interpolation for function)
        LXe.kinetic_profiles.beta = @(rr) reshape(interp1(Le.r_q, blended_beta_vec, rr(:), 'pchip', 'extrap'), size(rr));
        % compute derivative betap (approx from betap_fit blended similarly)
        if exist('betap_fit','var')
            prev_betap_vec = zeros(size(Le.r_q));
            if isfield(LXe_old,'kinetic_profiles') && isfield(LXe_old.kinetic_profiles,'betap') ...
                    && ~isempty(LXe_old.kinetic_profiles.betap)
                prev_betap_vec = LXe_old.kinetic_profiles.betap(Le.r_q);
            end
            new_betap_vec = betap_fit(Le.r_q);
            blended_betap_vec = (1-alpha_beta)*prev_betap_vec + alpha_beta*new_betap_vec;
            LXe.kinetic_profiles.betap = @(rr) reshape(interp1(Le.r_q, blended_betap_vec, rr(:), 'pchip', 'extrap'), size(rr));
        end
    end

    % preserve other settings
    LXe.eps_val = eps_val;
    LXe.Sbc(1) = 0; LXe.S2bc = 0; LXe.Sbc(2) = 0; LXe.S3bc = 0;
    %LXe.Sbc(1) = eqExp.S(end,1) / (a_minor * eps_val); LXe.S2bc = eqExp.S(end,1) / (a_minor * eps_val); 
    %LXe.Sbc(2) = eqExp.S(end,2) / (a_minor * eps_val); LXe.S3bc = eqExp.S(end,2) / (a_minor * eps_val);
    LXe.x = LYe.x; Le.P.hot_restart = true;

    % Try running equilY. If it errors or returns non-finite values, reduce alpha_beta and retry.
    max_retries = 4;
    retry = 0;
    success = false;
    current_alpha_beta = alpha_beta;
    LXe_try = LXe;  % copy for attempts

    while ~success && retry <= max_retries
        try
            LYe_try = equilY(Le, LXe_try);
            % basic sanity-check: ensure finite fields (adjust as appropriate)
            if all(isfinite(LYe_try.x(:)))
                success = true;
                LYe = LYe_try;
                % update LXe_old for next iteration
                LXe_old = LXe_try;
            else
                error('nonfinite_result');
            end
        catch
            % reduce step and retry
            retry = retry + 1;
            current_alpha_beta = current_alpha_beta / 2;
            % re-blend beta with a smaller alpha (if beta was ever set)
            if current_alpha_beta > 0 && isfield(LXe_old,'kinetic_profiles')
                prev_beta_vec = LXe_old.kinetic_profiles.beta(Le.r_q);
                new_beta_vec = beta_fit(Le.r_q);
                blended_beta_vec = (1-current_alpha_beta)*prev_beta_vec + current_alpha_beta*new_beta_vec;
                LXe_try.kinetic_profiles.beta = @(rr) reshape(interp1(Le.r_q, blended_beta_vec, rr(:), 'pchip', 'extrap'), size(rr));
                % update betap similarly if exists
                if exist('betap_fit','var')
                    prev_betap_vec = zeros(size(Le.r_q));
                    if isfield(LXe_old,'kinetic_profiles') && isfield(LXe_old.kinetic_profiles,'betap')
                        prev_betap_vec = LXe_old.kinetic_profiles.betap(Le.r_q);
                    end
                    new_betap_vec = betap_fit(Le.r_q);
                    blended_betap_vec = (1-current_alpha_beta)*prev_betap_vec + current_alpha_beta*new_betap_vec;
                    LXe_try.kinetic_profiles.betap = @(rr) reshape(interp1(Le.r_q, blended_betap_vec, rr(:), 'pchip', 'extrap'), size(rr));
                end
            end
            % if retries exhausted, fall back to applying only q (no beta)
            if retry == max_retries
                LXe_try = LXe; % reset
                if isfield(LXe_try,'kinetic_profiles')
                    LXe_try.kinetic_profiles = rmfield(LXe_try.kinetic_profiles,'beta');
                    if isfield(LXe_try.kinetic_profiles,'betap'), LXe_try.kinetic_profiles = rmfield(LXe_try.kinetic_profiles,'betap'); end
                end
                try
                    LYe_try = equilY(Le, LXe_try);
                    if all(isfinite(LYe_try.x(:)))
                        success = true;
                        LYe = LYe_try;
                        LXe_old = LXe_try;
                    end
                catch
                    % give up on this iteration: keep previous LYe and proceed
                    warning('equilY_failed','equilY failed even after retries; keeping previous state.');
                    LYe = LYe; % unchanged
                    success = true;
                end
            end
        end
    end

    % proceed to next iteration
end

A0 = 3;
[L2, LX2] = equilSol('beta', 1.1, 's0', 5, 'debug', 4, ...
                         'residuals_fun',@residuals_noRepl,...
                         'jacobian_fun', @jacobian_noRepl, 'Ns', 3,'Nb',10,...
                         'equation_of_state',@biMaxwellian,'A0',A0,'Bc0',0.98,...
                         'hot_restart', true,'damping',0.7);
eps_val = 0.37225;
LX2.eps_val = eps_val;
LX2.Sbc(1) = 0;  LX2.Sbc(2) = 0; 
LX2.qfun = LXe.qfun;
LX2.qpfun = LXe.qpfun;
LX2.x = LYe.x;
LX2.kinetic_profiles.beta = @(r) LXe.kinetic_profiles.beta(r)/A0*1.5;
LX2.kinetic_profiles.betap = @(r) LXe.kinetic_profiles.betap(r)/A0*1.5;
LY2 = equilY(L2, LX2);

LX2.x = LY2.x;L2.P.do_shift_NLO = true;
LY2 = equilY(L2, LX2);
L2.P.damping =1;

%% scan_eps_from equilibrium part
% Basically scan_eps_from_equilibrium written out, but for source control
L0 = L2; LX0 = LX2; LY0 = LY2;


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

    fprintf('Run %d/%d, eps = %.1e, |res| = %.4e, bp_LIU = %.1f, li_LIU = %.1f\n', ...
            ii+1, nsim, LX.eps_val, LY.res_norms(end), LY.bp_liu, LY.li_liu);

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
ylim([1.9e-4 3.83])
xlim([1e-3 max(values_eps)])
yticks([1e-4, 1e-3 1e-2 1e-1 1e0])
ax = nexttile(2);

h1 = plot_error_panel_2(L2S2error_all, values_eps, 0.02, cmap, col_inds, ...
    {'o'}, '$\hat S_{2,0}$');
hold on;
h2 = plot_error_panel_2(L2S2_1_error_all .* values_eps, values_eps, 0.02, cmap, col_inds, ...
    {'square'},  '$\hat S_{2,0}+\epsilon\hat S_{2,1}$');
hold on;
plot([0.146189 0.146189], [1e-4 1e2], '--', 'Color', 'r')
legend([h1 h2], 'Interpreter','latex', ...
       'FontSize',14, 'Box','off');


title('$\|\hat S_{2,\mathrm{num}}-\hat S_{2,\mathrm{ana}}\|_{L^2}$', ...
      'Interpreter','latex','FontSize',14);

xlabel('$\epsilon$', 'Interpreter', 'latex','FontSize',14)

axis(ax,'tight');
set(ax,'LooseInset',get(ax,'TightInset'));
ylim([6.3e-7 2.44])
xlim([1e-3 max(values_eps)])
yticks([1e-6 1e-5 1e-4, 1e-3 1e-2 1e-1 1e0])


%% profiles
figure;
t = tiledlayout(1,2,'TileSpacing','compact','Padding','none');
ax = nexttile(1);
plot_profile_panel(L0.r_q, delta0_ana, delta_num, cmap, col_inds, '$\Delta_0$ profiles', '$\hat\Delta$');
delete(ax.Title)
legend( {'$\textrm{analytical}$'}, 'Interpreter','latex','FontSize',12,'Box','off', 'Location', 'Northwest')
axis(ax,'tight'); set(ax,'LooseInset',get(ax,'TightInset'));
ax = nexttile(2);
plot_profile_panel(L0.r_q,S2_1_ana, S2_1_num, cmap, col_inds, '$\hat \Delta_1$ profiles', '$\hat S_{2}/\epsilon$')
delete(ax.Title)
legend( {'$\textrm{analytical}$'}, 'Interpreter','latex','FontSize',12,'Box','off', 'Location', 'Southwest')
axis(ax,'tight'); set(ax,'LooseInset',get(ax,'TightInset'));