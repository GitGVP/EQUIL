%% Simple circular case comparing theta_SFL grid vs. omega grid for (R,Z)
[L, LX] = equilSol('debug',4, 'Nb', 1, 'do_SFL', true);
LX.Sbc = zeros(L.P.Ns, 1);
LY = equilY(L, LX);
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
