%% Variational counterpart of examples/simple_case.m

beta_final = 30;
eps_val = 0.4;
Ns = 10;
Sbc = zeros(Ns,1);

%% Variational solve: fixed-epsilon beta continuation is internal
[L,LX] = equilVariationalSol('beta',beta_final,'debug',4,'Ns',Ns);
LX.eps_val = eps_val;
LX.Sbc = Sbc;
LY = equilVariationalY(L,LX);

assert(LY.isconverged);
assert(LY.local_B_residual < 1e-10);
fprintf(['Variational: beta=%g, |R|=%.3e, ', ...
    'max|epsilon Delta''|=%.4f\n'],beta_final,norm(LY.residual), ...
    max(abs(eps_val*LY.deltap)));
%% Flux surfaces and constant-omega curves, including an axis zoom
figure;hold on;
omega_indices = unique(floor(linspace(1,numel(LY.omega_plt),11)));
contour(LY.RR,LY.ZZ,LY.psiN.*ones(size(LY.RR)), ...
    linspace(0,1,11),'LineWidth',1.5);
for k = omega_indices
    plot(LY.RR(:,k),LY.ZZ(:,k),'--', ...
        'Color',[0,0,0,0.5],'LineWidth',0.5);
end
xlabel('$R/R_0$','Interpreter','latex','FontSize',12);
ylabel('$Z/R_0$','Interpreter','latex','FontSize',12);
title('$\psi_N$','Interpreter','latex','FontSize',12);
axis equal; grid on; box on;


figure;
plot(LY.r_plt, LY.S(:,:,1),'.')
grid on; box on;