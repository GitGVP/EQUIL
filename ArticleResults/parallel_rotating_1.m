%% Section 4.3: Strong parallel anisotropy with near sonic flows
% This script can be run as is to reproduce Figs. 6, 7 of the article
% Van Parys et al., 
% "Investigation of finite aspect ratio effects in axisymmetric 
% magnetic equilibria with toroidal rotation and pressure anisotropy"

% Parameters (shaping, aspect ratio, shear, beta, rotation)
s2bc = -0.8;
s3bc = -0.3;
eps_val = 0.34;
s0 = 4.5;
beta = 0.65; 
mach20 = 1; 

% 1: parallel rotating 2: isotropic rotating 3: pressureless 4: parallel
% no rotation
[L, LX] = equilSol('beta',beta, 's0', s0, 'debug', 4,...
    'residuals_fun',@residuals_noRepl,'jacobian_fun', @jacobian_noRepl, ...
    'equation_of_state', @parallel_rotating,'mach20',mach20,'Nb',10, 'Ns',4);
LX.eps_val = eps_val;
LX.Sbc(1) = s2bc;
LX.Sbc(2) = s3bc;
LY = equilY(L, LX);

[L2, LX2] = equilSol('beta',beta, 's0', s0, 'debug', 4,...
    'residuals_fun',@residuals_noRepl,'jacobian_fun', @jacobian_noRepl, ...
    'equation_of_state', @isotropic_rotating,'Ns', 4,'mach20',mach20,'Nb',1);
LX2.eps_val = eps_val;

LX2.Sbc(1) = s2bc;
LX2.Sbc(2) = s3bc;
LY2 = equilY(L2, LX2);

[L3, LX3] = equilSol('beta',0, 's0', s0, 'debug', 4,...
    'residuals_fun',@residuals_noRepl,'jacobian_fun', @jacobian_noRepl, ...
    'equation_of_state', @isotropic,'Ns', 4,'Nb',10);
LX3.eps_val = eps_val;
LX3.Sbc(1) = s2bc;
LX3.Sbc(2) = s3bc;
LY3 = equilY(L3, LX3);

[L4, LX4] = equilSol('beta',beta, 's0', s0, 'debug', 4,...
    'residuals_fun',@residuals_noRepl,'jacobian_fun', @jacobian_noRepl, ...
    'equation_of_state', @parallel_rotating,'mach20',0,'Nb',10, 'Ns',4);
LX4.eps_val = eps_val;
LX4.Sbc(1) = s2bc;
LX4.Sbc(2) = s3bc;
LX4.x = LY3.x;
LY4 = equilY(L4, LX4);


figure;
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(1); 
hold(ax1,'on');
h1 = plot(ax1, L3.r_q, LY3.gavg-1, '-','LineWidth',3); % no pressure first
h2 = plot(ax1, L.r_q,  LY.gavg-1, '--','LineWidth',3);
h3 = plot(ax1, L2.r_q, LY2.gavg-1, '-.','LineWidth',3);
h4 = plot(ax1, L4.r_q, LY4.gavg-1, '--','LineWidth',3);
ylabel(ax1,'$(\langle RB_\phi\rangle-R_0B_0)/R_0B_0$','Interpreter','latex','FontSize',14);
xlabel(ax1,'$\hat r$','Interpreter','latex','FontSize',14);
grid(ax1,'on');
ax2 = nexttile(2);
hold(ax2,'on');
plot(ax2, LY3.r_plt, LY3.delta, '-','LineWidth',3)
plot(ax2, LY.r_plt,  LY.delta, '--','LineWidth',3)
plot(ax2, LY2.r_plt, LY2.delta, '-.','LineWidth',3)
plot(ax2, LY4.r_plt, LY4.delta, '--','LineWidth',3)
ylabel(ax2,'$\hat\Delta$','Interpreter','latex','FontSize',14);
xlabel(ax2,'$\hat r$','Interpreter','latex','FontSize',14);
grid(ax2,'on');
ax3 = nexttile(3); 
hold(ax3,'on');
plot(ax3, L3.r_q, (squeeze(LY3.S(2:end-1,:,1)) - LY3.S2_ana)/eps_val, '-','LineWidth',3)
plot(ax3, L.r_q,  (squeeze(LY.S(2:end-1,:,1))   - LY.S2_ana)/eps_val, '--','LineWidth',3)
plot(ax3, L2.r_q, (squeeze(LY2.S(2:end-1,:,1))  - LY2.S2_ana)/eps_val, '-.','LineWidth',3)
plot(ax3, L4.r_q, (squeeze(LY4.S(2:end-1,:,1))  - LY4.S2_ana)/eps_val, '--','LineWidth',3)
ylabel(ax3,'$(\hat S_2-\hat S_2^{LO})/\epsilon$','Interpreter','latex','FontSize',14);
xlabel(ax3,'$\hat r$','Interpreter','latex','FontSize',14);
grid(ax3,'on');
sgtitle('$\qquad$','Interpreter','latex')
drawnow;

tlPos = tl.Position;  
titleHeight = 0.06;
axLegPos = [tlPos(1), tlPos(2) + tlPos(4) - titleHeight, tlPos(3), titleHeight];
axLeg = axes('Position', axLegPos, 'Visible', 'off', 'HitTest', 'off');
lg = legend(axLeg, [h1 h2 h3 h4], ...
    {'$p=0\qquad$', '$p_\perp=0,\ \mathcal M^2=1\qquad$', 'isotropic $\mathcal M^2=1\qquad$', '$p_\perp=0,\ \mathcal M^2=0\qquad$'}, ...
    'Interpreter','latex','Box','off','FontSize',12);
lg.Orientation = 'horizontal';
lg.NumColumns = 4;
lg.AutoUpdate = 'off';    % prevents legend from changing if you modify plots later
drawnow;                  % ensure lg.Position is up-to-date
legPos = lg.Position;     % relative to figure
axPos  = axLeg.Position;  % relative to figure
newX = axPos(1) + (axPos(3) - legPos(3))/2;
newY = axPos(2) + (axPos(4) - legPos(4))/2;
lg.Position = [newX, newY, legPos(3), legPos(4)];


figure;
tiledlayout(1,3,'TileSpacing','compact','Padding','none');
load TCV_limiter.mat
rl = [Grl;Grl(1)]; zl = [Gzl;Gzl(1)];
R0 = 0.92;
colormap parula;
ax4 = nexttile(1); 
axis(ax4,'equal'); hold(ax4,'on');
plot(ax4,rl,zl,'color',0.3*[1 1 1],'linewidth',3);
contourf(ax4, R0*LY4.RR, R0*LY4.ZZ, LY4.betapar,7,'LineColor','none');
for i=floor(linspace(1,numel(LY4.r_plt),7)); plot(ax4,R0*LY4.RR(i,:), R0*LY4.ZZ(i,:),'--','Color', [0,0,0,0.7],'LineWidth',0.5); end
plot(ax4,R0*LY4.RR(end,:), R0*LY4.ZZ(end,:), 'k', 'LineWidth',4);
title(ax4,'$\beta_\parallel$   $(\beta_\perp=0,\mathcal M^2=0)$','Interpreter','latex','FontSize',14);
colorbar(ax4,'eastoutside');
xlabel(ax4,'$R$ [m]','Interpreter','latex','FontSize',14);
ylabel(ax4,'$Z$ [m]','Interpreter','latex','FontSize',14);

R0 = 0.95;
ax1 = nexttile(2); 
axis(ax1,'equal'); hold(ax1,'on');
plot(ax1,rl,zl,'color',0.3*[1 1 1],'linewidth',3);
contourf(ax1, R0*LY.RR, R0*LY.ZZ, LY.betapar,7,'LineColor','none');
for i=floor(linspace(1,numel(LY.r_plt),7)); plot(ax1,R0*LY.RR(i,:), R0*LY.ZZ(i,:),'--','Color', [0,0,0,0.7],'LineWidth',0.5); end
plot(ax1,R0*LY.RR(end,:), R0*LY.ZZ(end,:), 'k', 'LineWidth',4);
title(ax1,'$\beta_\parallel$   $(\beta_\perp=0,\mathcal M^2=1)$','Interpreter','latex','FontSize',14);
colorbar(ax1,'eastoutside');
xlabel(ax1,'$R$ [m]','Interpreter','latex','FontSize',14);
ylabel(ax1,'$Z$ [m]','Interpreter','latex','FontSize',14);

ax5 = nexttile(3); 
R0 = 0.95;
axis(ax5,'equal'); hold(ax5,'on');
plot(ax5,rl,zl,'color',0.3*[1 1 1],'linewidth',3);
contourf(ax5, R0*LY2.RR, R0*LY2.ZZ, (LY2.betapar+LY2.betaperp)/2,7,'LineColor','none');
for i=floor(linspace(1,numel(LY2.r_plt),7)); plot(ax5,R0*LY2.RR(i,:), R0*LY2.ZZ(i,:),'--','Color', [0,0,0,0.7],'LineWidth',0.5); end
plot(ax5,R0*LY2.RR(end,:), R0*LY2.ZZ(end,:), 'k', 'LineWidth',4);
title(ax5,'$\beta$   $(\mathcal M^2=1)$','Interpreter','latex','FontSize',14);
colorbar(ax5,'eastoutside');
xlabel(ax5,'$R$ [m]','Interpreter','latex','FontSize',14);
ylabel(ax5,'$Z$ [m]','Interpreter','latex','FontSize',14);
clim(ax4,[0 0.6]);
clim(ax1,[0 0.6]);
clim(ax5,[0 0.6]);