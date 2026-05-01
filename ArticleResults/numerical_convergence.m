%% Section 2.3: Discretization & Numerical Implementation
% This script can be run as is to reproduce Fig. 1 of the article
% Van Parys et al., 
% "Investigation of finite aspect ratio effects in axisymmetric 
% magnetic equilibria with toroidal rotation and pressure anisotropy"

nsim = 12;
ms = round(linspace(9,31,nsim));
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
              assemble_FE_matrices_bspline_neumann(L.r_nodes, L.P.Nb, L.P.Ns, L.P.nq, L.P.spline_p);
    err_t2 = LY.t2(2:end-1)-LY.t2_ana;
    err_t2p = LY.t2p(2:end-1)-LY.t2p_ana;
    L2_t2_only(ii) = sqrt(err_t2.' * ML2 * err_t2);
    H1_t2_only(ii) = sqrt(err_t2p.' * ML2_q * err_t2p);    
    R = reshape(L.r_q, L.P.nq, [])';           % m x nq
    elem_len = max(R,[],2) - min(R,[],2);  % length of each element
    h_geom = exp(mean(log(elem_len))); 
    h_eff(ii) = h_geom;
end

figure;
tiledlayout(1,2,"TileSpacing","compact","Padding","compact")
nexttile;

p = polyfit(log(h_eff(:)), log(L2_t2_only(:)), 1);
slope = p(1);
h_fit = linspace(min(h_eff), max(h_eff), 200);
L2_fit = exp(polyval(p, log(h_fit)));
loglog(h_fit, L2_fit, '-', 'LineWidth',2, 'Color',[0.96, 0.30, 0.75]);
grid on; hold on;
loglog(h_eff, L2_t2_only, 'ok', 'LineWidth',2);
legend({sprintf('fit: slope = %.3g', slope)}, ...
    'Location','SouthEast','Box','off','Interpreter','latex','FontSize',14);
title('$\|t_2 - t_{2,h}\|_{L^2}$', 'Interpreter','latex','FontSize',14)
xlabel('$h$', 'Interpreter','latex','FontSize',14)
nexttile;
p = polyfit(log(h_eff(:)), log(H1_t2_only(:)), 1);
slope = p(1);
H1_fit = exp(polyval(p, log(h_fit)));
loglog(h_fit, H1_fit, '-', 'LineWidth',2, 'Color', [0.60, 0.20, 0.80]);
grid on; hold on;
loglog(h_eff, H1_t2_only, 'ok', 'LineWidth',2);
legend({sprintf('fit: slope = %.3g', slope)}, ...
    'Location','SouthEast','Box','off','Interpreter','latex','FontSize',14);
title('$|t_2 - t_{2,h}|_{H^1}$', 'Interpreter','latex','FontSize',14)
xlabel('$h$', 'Interpreter','latex','FontSize',14)