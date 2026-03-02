%% Convergence study with Leading Order (LO) residuals for all profiles
% Two different ways of computing the L2 squared norm: diagonal matrix with
% gaussian weights or via matrix inversion of the projection from
% quadrature values to degrees of freedom space. (both give the same result)
% ------------------- 
% The various definitions of h (grid spacing) lead to the same results
% 

nsim = 13;
ms = [round(linspace(3,15,nsim))];
L2delta = zeros(nsim,1);
L2delta_q = zeros(nsim,1);

h_eff = zeros(nsim,1);
for ii=1:nsim
    
    [L, LX, LY] = equilSol('debug',4,'residuals_fun', @residuals_LO, ...
    'jacobian_fun',@jacobian_LO, 'Ns', 1, 'Nb', 10, 'm', ms(ii), 'nq', 12, 'spline_p',6);
    [~, ...
              ~, ~, ~,...
              ~, ~, ~, ~, ...
              ~, ~, ~, ~, ...
              ML2, ML2_q] = ... % for L2 error computation
              assemble_FE_matrices_bspline(L.r_nodes, L.P.Nb, L.P.Ns, L.P.nq, L.P.spline_p);
    err_t2 = LY.t2(2:end)-LY.t2_ana;
    err_delta = LY.delta(2:end)-LY.delta_ana;
    err_P = LY.P(2:end)-LY.P_ana;
    err_S2 = LY.S(2:end,1,1) - LY.S2_ana;
    L2delta(ii) = sqrt(err_delta.' * ML2 * err_delta + ...
        err_t2.' * ML2 * err_t2 + ...
        err_P.' * ML2 * err_P + ...
        err_S2.' * ML2 * err_S2);
    L2delta_q(ii) = sqrt(err_delta.' * ML2_q * err_delta + ...
        err_t2.' * ML2_q * err_t2 + ...
        err_P.' * ML2_q * err_P + ...
        err_S2.' * ML2_q * err_S2);
    %h_eff(ii) = 1/ numel(L.r_q);
    %h_eff(ii) = max(diff(L.r_q)); % use max difference, does not change anything from above
    R = reshape(L.r_q, L.P.nq, [])';           % m x nq
    elem_len = max(R,[],2) - min(R,[],2);  % length of each element
    h_geom = exp(mean(log(elem_len))); 
    h_eff(ii) = h_geom;
end


figure;
loglog(h_eff, L2delta, 'o', 'LineWidth',3);
grid on; hold on;
loglog(h_eff, L2delta_q, 'x', 'LineWidth',3);
p = polyfit(log(h_eff(:)), log(L2delta(:)), 1);
slope = p(1);
h_fit = linspace(min(h_eff), max(h_eff), 200);
L2_fit = exp(polyval(p, log(h_fit)));
loglog(h_fit, L2_fit, '-', 'LineWidth',3);
legend({'$M^t A^{-1} M$','$W$', sprintf('fit: slope = %.3g', slope)}, 'Location', 'SouthEast',...
    'Box','off', 'Interpreter','latex', 'FontSize',14);
title('$L_2$ error (all profiles)', 'Interpreter','latex', 'FontSize',14)
xlabel('$h$', 'Interpreter','latex', 'FontSize',14)


%% Do t2 only
% This is done with spline order = 3, the problematic residuals (delta, P)
% are removed (but in fact this does not affect the t2 convergence, i.e. this 
% % can be done using residuals_LO and jacobian_LO without problem.)
% ----------
% NOTE: in the choice of ms we skip even values, these lead to a fit with
% values that oscillate between overestimated error and underestimated (without 
% changing the obtained slope).


nsim = 12;
ms = [round(linspace(9,31,nsim))];
L2_t2_only = zeros(nsim,1);
H1_t2_only = zeros(nsim,1);
h_eff = zeros(nsim,1);
for ii=1:nsim
    
    [L, LX, LY] = equilSol('debug',4,'residuals_fun', @residuals_LO, ...
    'jacobian_fun',@jacobian_LO, 'Ns', 1, 'Nb', 10, 'm', ms(ii), 'nq', 6, 'spline_p',3);
    [~, ...
              ~, ~, ~,...
              ~, ~, ~, ~, ...
              ~, ~, ~, ~, ...
              ML2, ML2_q] = ... % for L2 error computation
              assemble_FE_matrices_bspline(L.r_nodes, L.P.Nb, L.P.Ns, L.P.nq, L.P.spline_p);
    err_t2 = LY.t2(2:end)-LY.t2_ana;
    err_t2p = LY.t2p-LY.t2p_ana;
    L2_t2_only(ii) = sqrt(err_t2.' * ML2 * err_t2);
    H1_t2_only(ii) = sqrt(err_t2p.' * ML2_q * err_t2p);    
    %h_eff(ii) = 1/ numel(L.r_q);
    %h_eff(ii) = max(diff(L.r_q)); % use max difference, does not change anything from above
    R = reshape(L.r_q, L.P.nq, [])';           % m x nq
    elem_len = max(R,[],2) - min(R,[],2);  % length of each element
    h_geom = exp(mean(log(elem_len))); 
    h_eff(ii) = h_geom;
end

figure;
loglog(h_eff, L2_t2_only, 'o', 'LineWidth',3);
grid on; hold on;
p = polyfit(log(h_eff(:)), log(L2_t2_only(:)), 1);
slope = p(1);
h_fit = linspace(min(h_eff), max(h_eff), 200);
L2_fit = exp(polyval(p, log(h_fit)));
loglog(h_fit, L2_fit, '-', 'LineWidth',3);
legend({'$W$', sprintf('fit: slope = %.3g', slope)}, 'Location', 'SouthEast',...
    'Box','off', 'Interpreter','latex', 'FontSize',14);
title('$L_2$ error ($t_2$)', 'Interpreter','latex', 'FontSize',14)
xlabel('$h$', 'Interpreter','latex', 'FontSize',14)


figure;
tiledlayout(1,2,"TileSpacing","compact","Padding","compact")
nexttile;

p = polyfit(log(h_eff(:)), log(L2_t2_only(:)), 1);
slope = p(1);
h_fit = linspace(min(h_eff), max(h_eff), 200);
L2_fit = exp(polyval(p, log(h_fit)));
loglog(h_fit, L2_fit, '-', 'LineWidth',3);
grid on; hold on;
loglog(h_eff, L2_t2_only, 'o', 'LineWidth',3,'k');
legend({sprintf('fit: slope = %.3g', slope)}, ...
    'Location','SouthEast','Box','off','Interpreter','latex','FontSize',14);
title('$\|t_2 - t_{2,h}\|_{L^2}$', 'Interpreter','latex','FontSize',14)
xlabel('$h$', 'Interpreter','latex','FontSize',14)
nexttile;
p = polyfit(log(h_eff(:)), log(H1_t2_only(:)), 1);
slope = p(1);
h_fit = linspace(min(h_eff), max(h_eff), 200);
H1_fit = exp(polyval(p, log(h_fit)));
loglog(h_fit, H1_fit, '-', 'LineWidth',3);
grid on; hold on;
loglog(h_eff, H1_t2_only, 'o', 'LineWidth',3,'k');
legend({sprintf('fit: slope = %.3g', slope)}, ...
    'Location','SouthEast','Box','off','Interpreter','latex','FontSize',14);
title('$|t_2 - t_{2,h}|_{H^1}$', 'Interpreter','latex','FontSize',14)
xlabel('$h$', 'Interpreter','latex','FontSize',14)

%% Plots of the profiles 
% with tanh forcing, problems at r=0 are much much worse, this is
% illustrated in the following example and plots of profiles and first
% derivatives

[L, LX] = equilSol('debug',4,'residuals_fun', @residuals_LO, ...
    'jacobian_fun',@jacobian_LO, 'Ns', 1, 'Nb', 10);%, 'nq', 15, 'spline_p',15, 'm', 10);
r_ped = 0.6;
beta_ped = 0.1;
%betappedfun = @(r) beta_ped * (1./ cosh(8*pi*(r+r_ped)).^2-1 ./cosh(8*pi*(r-r_ped)).^2);
betappedfun = @(r) - 4 * pi * beta_ped ./ (1+ cosh(8*pi*(r-r_ped)));
LX.kinetic_profiles.beta = @(r) zeros(size(r)); LX.kinetic_profiles.betap = betappedfun; 
LY = equilY(L, LX);

figure;
tiledlayout(4,1,"TileSpacing",'compact','Padding','compact')
nexttile; hold on;
plot(L.r_q, LY.t2_ana, 'LineWidth', 1.5);
plot(LY.r_plt, LY.t2, 'o');
xlabel('$\hat r$', 'FontSize',14, 'Interpreter','latex');
ylabel('$t_2$', 'FontSize',14, 'Interpreter','latex');
grid on;
nexttile; hold on;
plot(L.r_q, LY.delta_ana, 'LineWidth', 1.5);
plot(LY.r_plt, LY.delta, 'o');
xlabel('$\hat r$', 'FontSize',14, 'Interpreter','latex');
ylabel('$\hat \Delta$', 'FontSize',14, 'Interpreter','latex');
grid on;
nexttile; hold on;
plot(L.r_q, LY.P_ana, 'LineWidth', 1.5);
plot(LY.r_plt, LY.P, 'o');
xlabel('$\hat r$', 'FontSize',14, 'Interpreter','latex');
ylabel('$\hat P$', 'FontSize',14, 'Interpreter','latex');
grid on;
nexttile; hold on;
plot(L.r_q, LY.S2_ana, 'LineWidth', 1.5);
plot(LY.r_plt, LY.S(:,1,1), 'o');
xlabel('$\hat r$', 'FontSize',14, 'Interpreter','latex');
ylabel('$\hat S_2$', 'FontSize',14, 'Interpreter','latex');
grid on;


%
figure;
tiledlayout(4,1,"TileSpacing",'compact','Padding','compact')
nexttile; hold on;
plot(L.r_q, LY.t2p_ana, 'LineWidth', 1.5);
plot(L.r_q, LY.t2p, 'o');
xlabel('$\hat r$', 'FontSize',14, 'Interpreter','latex');
ylabel("$t_2'$", 'FontSize',14, 'Interpreter','latex');
grid on;
nexttile; hold on;
plot(L.r_q, LY.deltap_ana, 'LineWidth', 1.5);
plot(LY.r_plt, LY.deltap, 'o');
xlabel('$\hat r$', 'FontSize',14, 'Interpreter','latex');
ylabel("$\hat \Delta$'", 'FontSize',14, 'Interpreter','latex');
grid on;
nexttile; hold on;
plot(L.r_q, LY.Pp_ana, 'LineWidth', 1.5);
plot(LY.r_plt, LY.Pp, 'o');
xlabel('$\hat r$', 'FontSize',14, 'Interpreter','latex');
ylabel("$\hat P'$", 'FontSize',14, 'Interpreter','latex');
grid on;
nexttile; hold on;
plot(L.r_q, interp1(LY.r_fine,LY.S2p_fine, L.r_q), 'LineWidth', 1.5);
plot(LY.r_plt, LY.Sp(:,1,1), 'o');
xlabel('$\hat r$', 'FontSize',14, 'Interpreter','latex');
ylabel("$\hat S_2'$", 'FontSize',14, 'Interpreter','latex');
grid on;



% p-refinement
nsim = 12;

spline_ps = 3:5;
nP = numel(spline_ps);

% m-grid allowed to depend on p
ms = repmat(round(linspace(9,31,nsim)).',1,nP);
%ms(:,2) = ms(:,2) - 2;
%ms(:,3) = ms(:,3) - 4; 
%ms(:,4) = ms(:,4) - 6;

ms(:,2) = linspace(7,18,12);
ms(:,3) = linspace(7,18,12);
%ms(:,4) = linspace(5,16,12);

% number of last points to discard for each p
% (can be tuned later; same length as spline_ps)
n_skip_tail = [0 0 0];   % example

L2_t2_only = zeros(nsim,nP);
h_eff      = zeros(nsim,nP);
slopes     = zeros(nP,1);

figure; hold on; grid on;
set(gca,'XScale','log','YScale','log');

fit_handles = gobjects(nP,1);

for ip = 1:nP
    
    spline_p = spline_ps(ip);
    nq = 2*spline_p+1;
    
    for ii = 1:nsim
        
        [L, LX, LY] = equilSol( ...
            'debug',4, ...
            'residuals_fun', @residuals_LO, ...
            'jacobian_fun',  @jacobian_LO, ...
            'Ns',1,'Nb',10, ...
            'm', ms(ii,ip), ...
            'nq', nq, ...
            'spline_p', spline_p);
        
        [~,~,~,~,~,~,~,~,~,~,~,~,ML2,~] = ...
            assemble_FE_matrices_bspline( ...
                L.r_nodes, L.P.Nb, L.P.Ns, L.P.nq, L.P.spline_p);
        
        err_t2 = LY.t2(2:end) - LY.t2_ana;
        L2_t2_only(ii,ip) = sqrt(err_t2.' * ML2 * err_t2);
        
        % effective mesh size (geometric mean element length)
        R = reshape(L.r_q, L.P.nq, [])';     % m x nq
        elem_len = max(R,[],2) - min(R,[],2);
        h_eff(ii,ip) = exp(mean(log(elem_len)));
    end
    
    % indices to keep
    n_keep = nsim - n_skip_tail(ip);
    idx = 1:n_keep;
    
    % fit using only retained points
    pfit = polyfit(log(h_eff(idx,ip)), log(L2_t2_only(idx,ip)), 1);
    slopes(ip) = pfit(1);
    
    h_fit  = linspace(min(h_eff(idx,ip)), max(h_eff(idx,ip)), 200);
    L2_fit = exp(polyval(pfit, log(h_fit)));
    
    % plot retained data only (no legend entry)
    h_data = loglog(h_eff(idx,ip), L2_t2_only(idx,ip), 'o', 'LineWidth',2);
    h_data.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    % plot fit
    fit_handles(ip) = loglog(h_fit, L2_fit, '-', ...
        'LineWidth',2, 'Color', h_data.Color);
end

% legend: fits only
leg = cell(nP,1);
for ip = 1:nP
    leg{ip} = sprintf('$p=%d$, slope = %.3g', spline_ps(ip), slopes(ip));
end

legend(fit_handles, leg, ...
    'Location','SouthEast', ...
    'Interpreter','latex','FontSize',13,'Box','off');

title('$L_2$ error ($t_2$)', 'Interpreter','latex','FontSize',14);
xlabel('$h$', 'Interpreter','latex','FontSize',14);