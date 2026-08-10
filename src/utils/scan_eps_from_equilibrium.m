function scan_eps_from_equilibrium(L0,LX0,LY0,varargin)
    if ~isfield(L0.P,'Nb')
        scan_eps_from_variational_equilibrium(L0,LX0,LY0,varargin{:});
        return
    end

    % defaults
    nsim_default = 10;
    skip_default = [0,0,0,0,0];

    % override in order, if provided
    switch numel(varargin)
        case 0
            nsim = nsim_default;
            skip_counts = skip_default;
        case 1
            nsim = varargin{1};
            skip_counts = skip_default;
        otherwise
            nsim = varargin{1};
            skip_counts = varargin{2};
    end

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
        delta1_num, t2_num, t2p_num, S2_num, deltap_num, S2_1_num] = deal(zeros(nsim, Nq));

    %[Rlfcs, Zlfcs] = deal(zeros(nsim, size(LY0.RR,2)));
    %Rlfcs(1,:) = LY0.RR(end,:); 
    %Zlfcs(1,:) = LY0.ZZ(end,:);
    delta0_ana = LY0.delta_ana;
    delta1_ana = interp1(LY0.r_fine, LY0.delta1_ana, L0.r_q, 'spline');
    S2_0_ana = LY0.S2_ana;
    S2_1_ana = interp1(LY0.r_fine, LY0.S2_1_fine, L0.r_q, 'spline');
    t2_ana = LY0.t2_ana;
    t2p_ana = LY0.t2p_ana;
    
    t2_num(1,:)    = LY0.t2(2:end-1);
    t2p_num(1,:)    = LY0.t2p(2:end-1);
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
    LX = LX0; L = L0; L.P.do_ana = false; L.P.do_shift_NLO = false;                     
    for ii = 1:nsim-1
        eps_val = values_eps(end-ii);
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
        t2p_num(ii+1,:)    = LY.t2p(2:end-1);
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
    cmap = jet(n_vals);
    log_eps = log10(eps_sorted);
    log_eps_norm = (log_eps - min(log_eps)) / (max(log_eps) - min(log_eps));
    col_inds = round(1 + log_eps_norm * (size(cmap,1) - 1));
    
    % plotting: adapt to single "outer" run
    figure;
    n_panels = 5;
    errors = {L2error_all, L20error_all, L2t2error_all, L2S2error_all, L2S2_1_error_all};
    titles = {'$L^2_{\Delta_1}$ error', '$L^{2}_{\Delta_0}$ error', '$L^2_{t_2}$ error', '$L^2_{S_2}$ error', '$L^2_{S_{2,1}}$ error'};
    outerloop_legend_name = ' ';
    for ip = 1:n_panels
        subplot(n_panels,1,ip);
        plot_error_panel(errors{ip}, values_eps, 0, cmap, col_inds, {'o'}, titles{ip}, outerloop_legend_name, skip_counts(ip));
    end
    
    figure;
    subplot(5,1,1);
    plot_profile_panel(L0.r_q, delta0_ana, delta_num, cmap, col_inds, '$\Delta_0$ profiles', '$\hat\Delta_0$');
    hold on; 
    h = plot(L0.r_q, delta0_ana + eps_max * delta1_ana...
        , 'LineWidth',2,'Color', [1 0.35 0.5]);
    legend(h, {'$\Delta_0 + \epsilon \Delta_1$'}, 'Interpreter','latex','FontSize',12,'Box','off', 'Location', 'Northwest')
    subplot(5,1,2);
    plot_profile_panel(L.r_q,delta1_ana, delta1_num, cmap, col_inds, '$\Delta_1$ profiles', '$\hat\Delta_1$')
    subplot(5,1,3);
    plot_profile_panel(L.r_q, t2_ana, t2_num, cmap, col_inds, '$t_2$ profiles', '$t_2$')
    subplot(5,1,4);
    plot_profile_panel(L.r_q, S2_0_ana, S2_num, cmap, col_inds, '$\hat S_2$ profiles', '$\hat S_2$')
    hold on; 
    h = plot(L0.r_q, S2_0_ana + eps_max * S2_1_ana, ...
              'LineWidth',2,'Color',[1 0.35 0.5]);
    legend(h, {'$S_{2,0} + \epsilon S_{2,1}$'}, 'Interpreter','latex','FontSize',12,'Box','off')
    subplot(5,1,5);
    plot_profile_panel(L0.r_q,S2_1_ana, S2_1_num, cmap, col_inds, '$\hat S_{2,1}$ profiles', '$\hat S_{2,1}$')

    % For the perp anis case
    
    
    figure(); 
    t = tiledlayout(2,2,'TileSpacing','compact','Padding','none'); % minimal spacing
    ax = nexttile(1);
    plot_error_panel(L20error_all, values_eps, 0, cmap, col_inds, {'o'}, '$L^{2}_{\Delta_0}$ error', outerloop_legend_name, 4);
    axis(ax,'tight');
    set(ax,'LooseInset',get(ax,'TightInset')); % removes extra padding
    ax = nexttile(3);
    plot_error_panel(L2error_all, values_eps, 0, cmap, col_inds, {'o'}, '$L^2_{\Delta_1}$ error', outerloop_legend_name, 4);
    axis(ax,'tight');
    set(ax,'LooseInset',get(ax,'TightInset')); % removes extra padding

    ax = nexttile(2);
    plot_profile_panel(L0.r_q, delta0_ana, delta_num, cmap, col_inds, '$\Delta_0$ profiles', '$\hat\Delta_0$');
    axis(ax,'tight'); set(ax,'LooseInset',get(ax,'TightInset'));
    
    ax = nexttile(4);
    plot_profile_panel(L0.r_q,S2_1_ana, S2_1_num, cmap, col_inds, '$\hat S_{2,1}$ profiles', '$\hat S_{2,1}$')
    axis(ax,'tight'); set(ax,'LooseInset',get(ax,'TightInset'));


    %% shift LO and NLO

    figure(); 
    t = tiledlayout(2,2,'TileSpacing','compact','Padding','none'); % minimal spacing
    ax = nexttile(1);
    plot_error_panel(L20error_all, values_eps, 0, cmap, col_inds, {'o'}, '$L^{2}_{\Delta_0}$ error', outerloop_legend_name, 4);
    axis(ax,'tight');
    set(ax,'LooseInset',get(ax,'TightInset')); % removes extra padding
    ax = nexttile(3);
    plot_error_panel(L2S2_1_error_all, values_eps, 0, cmap, col_inds, {'o'}, '$L^2_{S_{2,1}}$ error', outerloop_legend_name, 5);
    axis(ax,'tight');
    set(ax,'LooseInset',get(ax,'TightInset')); % removes extra padding

    ax = nexttile(2);
    plot_profile_panel(L0.r_q, delta0_ana, delta_num, cmap, col_inds, '$\Delta_0$ profiles', '$\hat\Delta_0$');
    axis(ax,'tight'); set(ax,'LooseInset',get(ax,'TightInset'));
    
    ax = nexttile(4);
    plot_profile_panel(L0.r_q,S2_1_ana, S2_1_num, cmap, col_inds, '$\hat S_{2,1}$ profiles', '$\hat S_{2,1}$')
    axis(ax,'tight'); set(ax,'LooseInset',get(ax,'TightInset'));



    % 

    % figure;
    % t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    % ax = nexttile(1);
    % plot_profile_panel(L.r_q, t2_ana, t2_num, cmap, col_inds, '$t_2$ profiles', '$t_2$')
    % axis(ax,'tight');
    % set(ax,'LooseInset',get(ax,'TightInset')); % removes extra padding
    % ax = nexttile(3);
    % plot_profile_panel(L.r_q, t2p_ana, t2p_num, cmap, col_inds, "$t_2'$ profiles", "$t_2'$")
    % axis(ax,'tight');
    % set(ax,'LooseInset',get(ax,'TightInset')); % removes extra padding
    % 
    % tmp = nexttile(2, [2 1]);   % create placeholder axes occupying right column
    % delete(tmp);                % remove it so we can place our own axes freely
    % %% Parameters for vertical stacking
    % tpos = t.Position;         % [x y width height] of tiledlayout
    % rx = tpos(1) + 0.5*tpos(3); % right column x
    % rw = 0.5*tpos(3);           % right column width
    % 
    % % fractions of right column height
    % conv_frac = 0.5;   % convergence plot height fraction
    % geom_frac = 0.35;  % geometries height fraction
    % cb_frac   = 0.08;  % colorbar height fraction
    % pad_v     = 0.02;  % vertical gaps between elements
    % 
    % % compute heights in normalized units
    % conv_h = conv_frac * tpos(4);
    % geom_h = geom_frac * tpos(4);
    % cb_h   = cb_frac   * tpos(4);
    % 
    % %% Position y coordinates
    % top_y = tpos(2) + tpos(4);         % top of right column
    % conv_y = top_y - conv_h - pad_v;   % convergence plot bottom
    % geom_y = conv_y - geom_h - pad_v;  % geometry plots bottom
    % cb_y   = tpos(2) + 0.01;           % colorbar bottom margin
    % 
    % %% Convergence axes
    % conv_x = rx + 0.05*rw;             % small horizontal padding
    % conv_w = rw - 0.04*rw;
    % 
    % ax_conv = axes('Units','normalized','Position',[conv_x, conv_y, conv_w, conv_h]);
    % hold(ax_conv,'on');
    % set(ax_conv, 'XScale', 'log', 'YScale', 'log');
    % grid(ax_conv,'on');
    % plot_error_panel(L2t2error_all, values_eps, 0, cmap, col_inds, {'o'}, '$L^{2}_{\Delta_0}$ error', outerloop_legend_name, 4);
    % plot_error_panel(H1error_all, values_eps, 0, cmap, col_inds, {'o'}, '$L^{2}_{\Delta_0}$ error', outerloop_legend_name, 4);
    % geom_x = conv_x;           % start left aligned
    % geom_w = conv_w / 10;      % width of each small axis
    % margin = 0.05;             % 5% margin to avoid clipping
    % for i = 1:10
    %     ax_geom = axes('Units','normalized','Position',[geom_x + (i-1)*geom_w, geom_y, geom_w, geom_h]);
    %     color_i = cmap(col_inds(11-i), :);        % small epsilon = red, large = blue
    %     plot(ax_geom, Rlcfs(1+2*i,:), Zlcfs(1+2*i,:), 'LineWidth',2, 'Color', color_i);
    % 
    %     axis(ax_geom,'equal'); 
    %     axis(ax_geom,'off');    
    %     xdata = Rlcfs(1+2*i,:);
    %     ydata = Zlcfs(1+2*i,:);
    %     xlim(ax_geom, [min(xdata)-margin*(max(xdata)-min(xdata)), max(xdata)+margin*(max(xdata)-min(xdata))]);
    %     ylim(ax_geom, [min(ydata)-margin*(max(ydata)-min(ydata)), max(ydata)+margin*(max(ydata)-min(ydata))]);
    % end
    % 
    % %% Horizontal colorbar below
    % cb = colorbar(ax_conv,'southoutside');  % attach to conv axes
    % cb.Orientation = 'horizontal';
    % cb.Units = 'normalized';
    % cb_w = conv_w;         % full width of right column
    % cb.Position = [conv_x, cb_y, cb_w, cb_h];
    % 
    % colormap(cmap_smooth);    % ensure smooth colormap is active
    % cb.Label.Interpreter = 'latex';
    % %cb.Label.String = '$\epsilon$';
    % cb.Ticks = linspace(0,1,5);
    % cb.TickLabels = arrayfun(@(x) sprintf('%.3f', 10^(min(log_eps)+(max(log_eps)-min(log_eps))*x)), cb.Ticks, 'UniformOutput', false);
 end



