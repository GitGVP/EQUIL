%% Variational counterpart of examples/rotating_case.m

[L,LX] = equilVariationalSol('equation_of_state',@isotropic_rotating, ...
    'mach20',1,'beta',0.3,'debug',4);
LX.eps_val = 0.3;
LY = equilVariationalY(L,LX);
assert(LY.isconverged);

figure; axis equal; hold on
colormap gray
contourf(LY.RR,LY.ZZ,LY.betapar, ...
         linspace(0,max(LY.betapar(:)),11),'LineColor','k');
for i = floor(linspace(1,numel(LY.r_plt),7))
    plot(LY.RR(i,:),LY.ZZ(i,:),'--','Color',[1,1,1,0.7], ...
         'LineWidth',1.5);
end
xlabel('$R/R_0$','Interpreter','latex','FontSize',14)
ylabel('$Z/R_0$','Interpreter','latex','FontSize',14)
title('$\epsilon^{-2}\mu_0p/B_0^2$','Interpreter','latex','FontSize',14)
colorbar

% The existing helper dispatches to the variational epsilon scan and uses
% the converged state only as a warm start for each new physical epsilon.
scan_eps_from_equilibrium(L,LX,LY,7)
