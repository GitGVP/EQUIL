%% Simple base case, to get familiar with I/O and plotting

% Runs with all default inputs besides debug = 4 which means that 
% Newton iterations are printed.

[L, LX, LY] = equilSol('debug',4);

figure;axis equal;hold on;
contour(LY.RR, LY.ZZ, ...
        [0;LY.psiN].*ones(size(LY.RR)), linspace(0,1,11), ...
        'LineWidth',2, 'ShowText',1);
for i=floor(linspace(1,numel(LY.omega_plt),11)); plot(LY.RR(:,i), LY.ZZ(:,i),'--','Color', [0,0,0,0.5],'LineWidth',0.5); end
xlabel('$R/R_0$', 'Interpreter','latex', 'FontSize',14)
ylabel('$Z/R_0$', 'Interpreter','latex', 'FontSize',14)
title('$\psi_N$', 'Interpreter','latex', 'FontSize',14)