%% TODO: Results section: parallel anisotropic rotating
% At the very least show the plasma diamagnetism, compared to the LO
% Maybe also compare NLO elongation

% So the idea is to show 3 cases: pressure-free,isotropic rotating,
% parallel-rotating, we'll show :
% 1) Plasma diamagnetism
% 2) shift of pressure contours ?
% 3) NLO elongation
% 4) delta ?



%% Do the 3 runs
% 1: parallel rotating 2: isotropic rotating 3: pressureless 4: parallel
% no rotation

eps_val = 0.35;
s0 = 2.5;
beta = 0.6;
mach20 = 0.7;

[L, LX] = equilSol('beta',beta, 's0', s0, 'debug', 4,...
    'residuals_fun',@residuals_vec,'jacobian_fun', @jacobian_vec, ...
    'equation_of_state', @parallel_rotating,'mach20',mach20,'Nb',10, 'Ns',2);
LX.eps_val = eps_val;
LY = equilY(L, LX);

[L2, LX2] = equilSol('beta',beta, 's0', s0, 'debug', 4,...
    'residuals_fun',@residuals_vec,'jacobian_fun', @jacobian_vec, ...
    'equation_of_state', @isotropic_rotating,'Ns', 2,'mach20',mach20,'Nb',1);
LX2.eps_val = eps_val;
LY2 = equilY(L2, LX2);

[L3, LX3] = equilSol('beta',0, 's0', s0, 'debug', 4,...
    'residuals_fun',@residuals_rotation,'jacobian_fun', @jacobian_rotation, ...
    'equation_of_state', @isotropic,'Ns', 2,'Nb',1);
LX3.eps_val = eps_val;
LY3 = equilY(L3, LX3);

[L4, LX4] = equilSol('beta',beta, 's0', s0, 'debug', 4,...
    'residuals_fun',@residuals_rotation,'jacobian_fun', @jacobian_rotation, ...
    'equation_of_state', @parallel_rotating,'mach20',0,'Nb',10, 'Ns',2);
LX4.eps_val = eps_val;
LY4 = equilY(L4, LX4);

figure;
subplot(3,1,1);
hold on;
plot(L3.r_q, LY3.gavg-1, '-','LineWidth',2) %no pressure first
plot(L.r_q, LY.gavg-1, '.') % anis rotating
plot(L2.r_q, LY2.gavg-1, '.') %isotropic rotating
plot(L4.r_q, LY4.gavg-1, '.')
ylabel('$(\langle RB_\phi\rangle-R_0B_0)/R_0B_0$','Interpreter','latex','FontSize',14);
legend({ '$p=0$','$p_\perp=0, \mathcal M^2=0.7',...
    'isotropic $\mathcal M^2=0.7$','$p_\perp=0, \mathcal M^2=0'},'box','off', 'Interpreter', 'latex',...
    'Fontsize',14, 'location','southwest')
grid on;
subplot(3,1,2);
hold on;
plot(LY3.r_plt, LY3.delta, '-','LineWidth',2)
plot(LY.r_plt, LY.delta, '.')
plot(LY2.r_plt, LY2.delta, '.')

plot(LY4.r_plt, LY4.delta, '.')
ylabel('$\hat\Delta$','Interpreter','latex','FontSize',14);
grid on;
subplot(3,1,3);
hold on;
plot(L3.r_q, squeeze(LY3.S(2:end,:,1))-LY3.S2_ana, '-','LineWidth',2)
plot(L.r_q, squeeze(LY.S(2:end,:,1))-LY.S2_ana, '.')
plot(L2.r_q, squeeze(LY2.S(2:end,:,1))-LY2.S2_ana, '.')

plot(L4.r_q, squeeze(LY4.S(2:end,:,1))-LY4.S2_ana, '.')
ylabel('$\hat S_2$','Interpreter','latex','FontSize',14);
grid on;


figure; axis equal;hold on;
contourf(LY4.RR, LY4.ZZ, LY4.betapar, 'LineColor', 'none')
for i=1:20:numel(LY4.r_plt)
    plot(LY4.RR(i,:), LY4.ZZ(i,:), 'w')
end

tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

%% Left column (three stacked plots) -- capture handles for legend
ax1 = nexttile(1);
hold(ax1,'on');
h1 = plot(ax1, L3.r_q, LY3.gavg-1, '-','LineWidth',2); % no pressure first
h2 = plot(ax1, L.r_q,  LY.gavg-1,  '.');
h3 = plot(ax1, L2.r_q, LY2.gavg-1, '.');
h4 = plot(ax1, L4.r_q, LY4.gavg-1, '.');
ylabel(ax1,'$(\langle RB_\phi\rangle-R_0B_0)/R_0B_0$','Interpreter','latex','FontSize',14);
grid(ax1,'on');

ax2 = nexttile(4);
hold(ax2,'on');
plot(ax2, LY3.r_plt, LY3.delta, '-','LineWidth',2)
plot(ax2, LY.r_plt,  LY.delta,  '.')
plot(ax2, LY2.r_plt, LY2.delta, '.')
plot(ax2, LY4.r_plt, LY4.delta, '.')
ylabel(ax2,'$\hat\Delta$','Interpreter','latex','FontSize',14);
grid(ax2,'on');

ax3 = nexttile(7);
hold(ax3,'on');
plot(ax3, L3.r_q, squeeze(LY3.S(2:end,:,1)) - LY3.S2_ana, '-','LineWidth',2)
plot(ax3, L.r_q,  squeeze(LY.S(2:end,:,1))   - LY.S2_ana, '.')
plot(ax3, L2.r_q, squeeze(LY2.S(2:end,:,1))  - LY2.S2_ana, '.')
plot(ax3, L4.r_q, squeeze(LY4.S(2:end,:,1))  - LY4.S2_ana, '.')
ylabel(ax3,'$\hat S_2$','Interpreter','latex','FontSize',14);
grid(ax3,'on');

%% Middle column tall contour (span rows 1:2, col 2) using LY4
ax4 = nexttile(2,[2 1]);  % shorter tall plot (2 rows high)
axis(ax4,'equal'); hold(ax4,'on');
contourf(ax4, LY4.RR, LY4.ZZ, (LY4.betapar+LY4.betaperp)/2,20, 'LineColor','none')
for i = 1:20:size(LY4.RR,1)
    plot(ax4, LY4.RR(i,:), LY4.ZZ(i,:), 'w')
end
title(ax4,'$(\beta_\parallel + \beta_\perp)/2$         $(\mathcal M^2=0)$','Interpreter','latex','FontSize',14);
colorbar(ax4,'eastoutside');
xlabel(ax4,'$R/R_0$','Interpreter','latex','FontSize',14);

%% Right column tall contour (span rows 1:2, col 3) using LY2
ax5 = nexttile(3,[2 1]);  % shorter tall plot (2 rows high)
axis(ax5,'equal'); hold(ax5,'on');
contourf(ax5, LY2.RR, LY2.ZZ, (LY2.betapar+LY2.betaperp)/2,20, 'LineColor','none')
for i = 1:20:size(LY2.RR,1)
    plot(ax5, LY2.RR(i,:), LY2.ZZ(i,:), 'w')
end
title(ax5,'$(\beta_\parallel + \beta_\perp)/2$         $(\mathcal M^2=0.7)$','Interpreter','latex','FontSize',14);
colorbar(ax5,'eastoutside');




% 6x3 layout. tall contours span rows 1..5. legend occupies row 6 cols 2..3.
tiledlayout(6,3,'TileSpacing','compact','Padding','compact');

%% Left column: three stacked small plots, each 2 rows high
ax1 = nexttile(1,[2 1]); % rows 1-2, col 1
hold(ax1,'on');
h1 = plot(ax1, L3.r_q, LY3.gavg-1, '-','LineWidth',3); % no pressure first
h2 = plot(ax1, L.r_q,  LY.gavg-1, '--','LineWidth',3);
h3 = plot(ax1, L2.r_q, LY2.gavg-1, '-.','LineWidth',3);
h4 = plot(ax1, L4.r_q, LY4.gavg-1, '--','LineWidth',3);
ylabel(ax1,'$(\langle RB_\phi\rangle-R_0B_0)/R_0B_0$','Interpreter','latex','FontSize',14);
grid(ax1,'on');

ax2 = nexttile(7,[2 1]); % rows 3-4, col 1
hold(ax2,'on');
plot(ax2, LY3.r_plt, LY3.delta, '-','LineWidth',3)
plot(ax2, LY.r_plt,  LY.delta, '--','LineWidth',3)
plot(ax2, LY2.r_plt, LY2.delta, '-.','LineWidth',3)
plot(ax2, LY4.r_plt, LY4.delta, '--','LineWidth',3)
ylabel(ax2,'$\hat\Delta$','Interpreter','latex','FontSize',14);
grid(ax2,'on');

ax3 = nexttile(13,[2 1]); % rows 5-6, col 1
hold(ax3,'on');
plot(ax3, L3.r_q, squeeze(LY3.S(2:end,:,1)) - LY3.S2_ana, '-','LineWidth',3)
plot(ax3, L.r_q,  squeeze(LY.S(2:end,:,1))   - LY.S2_ana, '--','LineWidth',3)
plot(ax3, L2.r_q, squeeze(LY2.S(2:end,:,1))  - LY2.S2_ana, '-.','LineWidth',3)
plot(ax3, L4.r_q, squeeze(LY4.S(2:end,:,1))  - LY4.S2_ana, '--','LineWidth',3)
ylabel(ax3,'$\hat S_2-\hat S_2^{LO}$','Interpreter','latex','FontSize',14);
xlabel(ax3,'$\hat r$','Interpreter','latex','FontSize',14);
grid(ax3,'on');

%% Middle column tall contour (rows 1..5, col 2) using LY4
ax4 = nexttile(2,[5 1]); % rows 1-5, col 2
axis(ax4,'equal'); hold(ax4,'on');
contourf(ax4, LY4.RR, LY4.ZZ, (LY4.betapar+LY4.betaperp)/2,20,'LineColor','none');
for i = 1:20:size(LY4.RR,1)
    plot(ax4, LY4.RR(i,:), LY4.ZZ(i,:), 'w')
end
title(ax4,'$(\beta_\parallel + \beta_\perp)/2$   $(\mathcal M^2=0)$','Interpreter','latex','FontSize',14);
colorbar(ax4,'eastoutside');
xlabel(ax4,'$R/R_0$','Interpreter','latex','FontSize',14);
ylabel(ax4,'$Z/R_0$','Interpreter','latex','FontSize',14);
%% Right column tall contour (rows 1..5, col 3) using LY2
ax5 = nexttile(3,[5 1]); % rows 1-5, col 3
axis(ax5,'equal'); hold(ax5,'on');
contourf(ax5, LY2.RR, LY2.ZZ, (LY2.betapar+LY2.betaperp)/2,20,'LineColor','none');
for i = 1:20:size(LY2.RR,1)
    plot(ax5, LY2.RR(i,:), LY2.ZZ(i,:), 'w')
end
title(ax5,'$(\beta_\parallel + \beta_\perp)/2$   $(\mathcal M^2=0.7)$','Interpreter','latex','FontSize',14);
colorbar(ax5,'eastoutside');
xlabel(ax5,'$R/R_0$','Interpreter','latex','FontSize',14);
%% Bottom row: legend tile centered and horizontal (row 6 cols 2-3)
axLeg = nexttile(17,[1 2]); % row6 col2-3
axis(axLeg,'off');

labels = {'$p=0$', '$p_\perp=0, \mathcal M^2=0.7$', 'isotropic $\mathcal M^2=0.7$', '$p_\perp=0, \mathcal M^2=0$'};
leg = legend(axLeg,[h1 h2 h3 h4], labels, ...
    'Orientation','horizontal','Interpreter','latex','Box','off','FontSize',12, 'NumColumns',2);

% center the legend inside the tile by adjusting positions
drawnow;                   % ensure positions are set
axPos  = axLeg.Position;   % [x y w h] in normalized figure units
legPos = leg.Position;     % [x y w h]
newX = axPos(1) + 0.5*(axPos(3) - legPos(3));
newY = axPos(2) + 0.5*(axPos(4) - legPos(4));
leg.Position = [newX newY legPos(3) legPos(4)];

% final style tweak
set([ax1 ax2 ax3 ax4 ax5],'FontSize',12);
