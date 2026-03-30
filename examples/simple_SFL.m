%% Simple circular case comparing theta_SFL grid vs. omega grid for (R,Z)
[L, LX] = equilSol('debug',4, 'Nb', 1, 'do_SFL', true,'q0',0.8);
LX.Sbc = zeros(L.P.Ns, 1);
LX.eps_val=0.1;
LY = equilY(L, LX);

% Y0 

figure;hold on;
plot(LY.r_plt, real(LY.Y0), '.')
plot(LY.r_plt, LX.eps_val * LY.Y0_LO, '.')
grid on;
xlabel('$\hat r$', 'Interpreter', 'latex', 'Fontsize',14)
legend({'$Y_0$', 'LO'}, 'Interpreter', 'latex', 'Fontsize',14, 'Box','off')

figure;
tiledlayout(1,2,"TileSpacing","tight","Padding","tight")
nexttile;hold on;
for i = 1:10:numel(LY.r_plt)
    plot(LY.RR_sfl(i,:), LY.ZZ_sfl(i,:), '--k');
end
for i=floor(linspace(1,numel(LY.omega_plt),21))
    plot(LY.RR_sfl(:,i), LY.ZZ_sfl(:,i),'--','Color', [0,0,0,1],'LineWidth',1); 
end
axis equal;
title('$\hat r =$ const, $\theta_{\mathrm{SFL}} =$ const', 'Interpreter', 'latex', 'Fontsize',14)
xlabel('$R/R_0$', 'Interpreter', 'latex', 'Fontsize',14)
ylabel('$Z/R_0$', 'Interpreter', 'latex', 'Fontsize',14)
nexttile;hold on;
for i = 1:10:numel(LY.r_plt)
    plot(LY.RR(i,:), LY.ZZ(i,:), '--k');
end
for i=floor(linspace(1,numel(LY.omega_plt),21))
    plot(LY.RR(:,i), LY.ZZ(:,i),'--','Color', [0,0,0,1],'LineWidth',1); 
end
title('$\hat r =$ const, $\omega =$ const', 'Interpreter', 'latex', 'Fontsize',14)
xlabel('$R/R_0$', 'Interpreter', 'latex', 'Fontsize',14)
axis equal;


figure;
tiledlayout(1,3,"TileSpacing","tight","Padding","tight")
nexttile;
contourf(LY.RR, LY.ZZ, LY.N_SFL)
axis equal;
xlabel('$R/R_0$', 'Interpreter', 'latex', 'Fontsize',14)
ylabel('$Z/R_0$', 'Interpreter', 'latex', 'Fontsize',14)
colorbar;
title('$N=g_{\vartheta\vartheta}/J_\vartheta$', 'Interpreter', 'latex', 'Fontsize',18)
nexttile;
contourf(LY.RR, LY.ZZ, LY.M_SFL)
axis equal;
xlabel('$R/R_0$', 'Interpreter', 'latex', 'Fontsize',14)
ylabel('$Z/R_0$', 'Interpreter', 'latex', 'Fontsize',14)
colorbar;
title('$M=g_{r\vartheta}/J_\vartheta$', 'Interpreter', 'latex', 'Fontsize',18)
nexttile;
contourf(LY.RR, LY.ZZ, LY.L_SFL.*LY.r_plt)
axis equal;
xlabel('$R/R_0$', 'Interpreter', 'latex', 'Fontsize',14)
ylabel('$Z/R_0$', 'Interpreter', 'latex', 'Fontsize',14)
colorbar;
title('$\hat rL=\hat rg_{rr}/J_\vartheta$', 'Interpreter', 'latex', 'Fontsize',18)

figure;
contourf(LY.RR, LY.ZZ, LY.N_ana)
axis equal;
xlabel('$R/R_0$', 'Interpreter', 'latex', 'Fontsize',14)
ylabel('$Z/R_0$', 'Interpreter', 'latex', 'Fontsize',14)
colorbar;
title('$N=g_{\vartheta\vartheta}/J_\vartheta$ (analytical)', 'Interpreter', 'latex', 'Fontsize',18)


figure;
contourf(LY.RR, LY.ZZ, LY.N_ana2)
axis equal;
xlabel('$R/R_0$', 'Interpreter', 'latex', 'Fontsize',14)
ylabel('$Z/R_0$', 'Interpreter', 'latex', 'Fontsize',14)
colorbar;
title('second $N=g_{\vartheta\vartheta}/J_\vartheta$ (analytical)', 'Interpreter', 'latex', 'Fontsize',18)


figure;hold on;
plot(LY.r_plt, real(LY.Nm1), '.')
plot(LY.r_plt, real(LY.Mm1), '.')
plot(LY.r_plt, imag(LY.Nm1), '-')
plot(LY.r_plt, imag(LY.Mm1), '-')
grid on;
xlabel('$\hat r$', 'Interpreter', 'latex', 'Fontsize',14)
legend({'$\Re{(N_{-1})}$','$\Re{(M_{-1})}$','$\Im{(N_{-1})}$','$\Im{(M_{-1})}$'}, 'Interpreter', 'latex', 'Fontsize',14, 'Box','off')



% Mm1

figure;hold on;
plot(LY.r_plt, imag(LY.Mm1), '.')
plot(LY.r_plt, LX.eps_val/2 * (LY.deltapp + LY.deltap ./ LY.r_plt + 1) .* LY.r_plt, '.')
grid on;
xlabel('$\hat r$', 'Interpreter', 'latex', 'Fontsize',14)
legend({'$\Im{(M_{-1})}$', 'LO'}, 'Interpreter', 'latex', 'Fontsize',14, 'Box','off')


figure;hold on;
plot(LY.r_plt, real(LY.Nm1), '.')
plot(LY.r_plt, LX.eps_val^2 * LY.r_plt .* LY.deltap, '.')
grid on;
xlabel('$\hat r$', 'Interpreter', 'latex', 'Fontsize',14)
legend({'$\Re{(N_{-1})}$', 'LO'}, 'Interpreter', 'latex', 'Fontsize',14, 'Box','off')




%% Solve for the upper-sideband

r=LY.r_fine;
[~, I] = min(sqrt((LX.qfun(r)-1).^2));
[~, J] = min(sqrt((LX.qfun(r)-2).^2));
[X2_i, X2_ip] = solve_Xi_eq(r(1:I), 0, 1, LX.qfun, LX.qpfun);
[X2_e, X2_ep] = solve_Xi_eq(r(I:J), 1, 0, LX.qfun, LX.qpfun); 
% TODO, solve this with Neumann outer BC

V2_fun = @(rr) interp1(LY.r_plt, LY.V2, rr, 'spline');

[X2_zeta, X2_zetap] = solve_Xi_eq(r(1:I), 0, 0, LX.qfun, LX.qpfun, V2_fun);
A = (LY.dxidrjump + X2_zetap(end)) / (X2_ep(1) - X2_ip(end));
xi_in = X2_zeta + A * X2_i;
xi_out = A * X2_e;
figure; hold on;
plot(r(1:I), xi_in, '.')
plot(r(I:J), xi_out, '.')

% check that we do have that jump
xi2ip = gradient(xi_in,r(1:I));
figure; hold on;
plot(r(1:I),  xi2ip-xi2ip(end))
plot(r(I:J-100), gradient(xi_out(1:end-100),r(I:J-100))-xi2ip(end))
plot(r, ones(size(r)) * LY.dxidrjump)


figure; hold on;
plot(r(1:I), X2_i, '.')
plot(r(I:J), X2_e, '.')
plot(r(1:I), X2_zeta, '.')