%% Optional up-down asymmetric shaping
% V(:,:,1) and V(:,:,2) are Fitzpatrick V_2 and V_3.  Setting Nh=0
% removes these profiles and recovers the symmetric variational system.

epsilon = 0.35;
Ns = 3;
Nh = 2;

[L,LX] = equilVariationalSol('beta',0,'Ns',Ns,'Nh',Nh, ...
    'm',4,'nq',10,'om_pts',96, ...
    'NLtol',1e-10,'debug',4, 'nk', 60);
LX.eps_val = epsilon;
LX.Sbc = [-0.4;0.1;0];
LX.Vbc = [0;0.6];
LY = equilVariationalY(L,LX);

% target_Vbc = [0.3;0];
% for fraction = [0.1,0.25,0.5,0.75,1]
%     LX.Vbc = fraction*target_Vbc;
%     LX.x = LY.x;
%     LX.local_B_guess = LY.local_B_quadrature;
%     LY = equilVariationalY(L,LX);
% end

%assert(LY.isconverged);
%assert(max(abs(squeeze(LY.V(end,1,:))-LX.Vbc)) < 1e-12);
up_down_asymmetry = max(hypot( ...
    LY.RR-fliplr(LY.RR),LY.ZZ+fliplr(LY.ZZ)),[],'all');
fprintf(['Asymmetric case: |R|=%.3e, min J/(epsilon^2 r)=%.4f, ', ...
    'up-down mismatch=%.3e\n'],norm(LY.residual), ...
    LY.min_J_over_eps2r,up_down_asymmetry);

figure;
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on; axis equal; grid on
contour(LY.RR,LY.ZZ,LY.psiN.*ones(size(LY.RR)), ...
    linspace(0,1,11),'LineWidth',1.5);
omega_indices = unique(floor(linspace(1,numel(LY.omega_plt),13)));
for k = omega_indices
    plot(LY.RR(:,k),LY.ZZ(:,k),'--','Color',[0,0,0,0.35], ...
        'LineWidth',0.6);
end
xlabel('$R/R_0$','Interpreter','latex');
ylabel('$Z/R_0$','Interpreter','latex');
title('up-down asymmetric flux surfaces');

nexttile; hold on; grid on
for is = 1:Ns
    plot(LY.r_plt,squeeze(LY.S(:,:,is)),'LineWidth',1.7, ...
        'DisplayName',sprintf('$S_%d$',is+1));
end
for iv = 1:Nh
    plot(LY.r_plt,squeeze(LY.V(:,:,iv)),'--','LineWidth',1.7, ...
        'DisplayName',sprintf('$V_%d$',iv+1));
end
xlabel('$\hat r$','Interpreter','latex');
ylabel('normalized shaping profile','Interpreter','latex');
legend('Interpreter','latex','Location','best','Box','off');
