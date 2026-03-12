%% Case with rotation, illustrating non-default IO, aspect ratio scans


% Run with equation of state = isotropic_rotating
[L, LX, LY] = equilSol('equation_of_state', @isotropic_rotating, ...
    'mach20',1, 'debug',4);

% default: constant rotation mach number, quadratic pressure profile

% Plot the shift between pressure contours and flux surfaces
figure;axis equal;hold on;
colormap gray
contourf(LY.RR, LY.ZZ, ...
        LY.betapar, linspace(0,max(LY.betapar(:)),11), 'LineColor', 'k');
for i=floor(linspace(1,numel(LY.r_plt),7)); plot(LY.RR(i,:), LY.ZZ(i,:),'--','Color', [1,1,1,0.7],'LineWidth',1.5); end
xlabel('$R/R_0$', 'Interpreter','latex', 'FontSize',14)
ylabel('$Z/R_0$', 'Interpreter','latex', 'FontSize',14)
title('$\epsilon^{-2}\mu_0 p/B_0^2$', 'Interpreter','latex', 'FontSize',14)
colorbar();

% Scan in aspect ratio and compare to asymptotic expansion

% First rerun just to compute the Next to Leading Order elongation
LX.x = LY.x;L.P.do_shift_NLO = true;
LY = equilY(L, LX);
scan_eps_from_equilibrium(L,LX,LY,20)
% 3 plots: 
% 1. convergence rates for shift, diamagnetism and elongation
% 2. Profiles (colors for each simulation), with expansion prediction (black)
% and pink is LO + max_eps * NLO prediction
% 3. Only shift and and NLO elongation.