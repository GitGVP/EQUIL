%% Check whether anisotropic case gives the same as before now.
[L, LX] = equilSol('beta',0.6,'s0', 4, 'debug', 4,...
    'residuals_fun',@residuals_noRepl,'jacobian_fun', @jacobian_noRepl, ...
    'equation_of_state', @isotropic,'Ns', 3,'Nb',10);

LX.eps_val = 0.3;
LX.Sbc(1) = -0.8;
LX.S2bc = -0.8;
LX.S3bc = 0.2;
LX.Sbc(2) = 0.2;
LY = equilY(L, LX);

[L, LX] = equilSol('beta',0.3, 'A0',2.5,'s0', 4, 'debug', 4,...
    'residuals_fun',@residuals_noRepl,'jacobian_fun', @jacobian_noRepl, ...
    'equation_of_state', @biMaxwellian,'Ns', 3,'Nb',10, 'Bc0', 1.2,...
    'damping',1);
LX.x = LY.x; L.P.hot_restart = true;LX.eps_val = 0.3;
LY = equilY(L, LX);

tiledlayout(1,3,'TileSpacing','compact','Padding','compact')

% Compute common color limits for the last two plots
clim_common = [min([LY.betapar(:); LY.betaperp(:)]), max([LY.betapar(:); LY.betaperp(:)])];

% --- Plot 1 ---
nexttile; hold on;
contourf(LY.RR, LY.ZZ, LY.BB, 'LineColor', 'none');

for i=1:40:numel(LY.r_plt)
    plot(LY.RR(i,:), LY.ZZ(i,:), 'w')
end
contour(LY.RR, LY.ZZ, LY.BB, [L.P.Bc0 L.P.Bc0], 'r');
cb1 = colorbar;
axis equal;
xlabel('$R$','Interpreter','latex','FontSize',14);
ylabel('$Z$','Interpreter','latex','FontSize',14);
title('$B/B_0$','Interpreter','latex','FontSize',14); 

% --- Plot 2 ---
nexttile; hold on;
contourf(LY.RR, LY.ZZ, LY.betapar.*ones(size(LY.ZZ)), 'LineColor', 'none');
for i=1:40:numel(LY.r_plt)
    plot(LY.RR(i,:), LY.ZZ(i,:), 'w')
end
cb2 = colorbar;
caxis(clim_common); % apply common limits
axis equal;
xlabel('$R$','Interpreter','latex','FontSize',14);
ylabel('$Z$','Interpreter','latex','FontSize',14);
title('$\beta_\parallel$','Interpreter','latex','FontSize',14);

% --- Plot 3 ---
nexttile; hold on;
contourf(LY.RR, LY.ZZ, LY.betaperp.*ones(size(LY.ZZ)), 'LineColor', 'none');
for i=1:40:numel(LY.r_plt)
    plot(LY.RR(i,:), LY.ZZ(i,:), 'w')
end
cb3 = colorbar;
caxis(clim_common); % apply common limits
axis equal;
xlabel('$R$','Interpreter','latex','FontSize',14);
ylabel('$Z$','Interpreter','latex','FontSize',14);
title('$\beta_\perp$','Interpreter','latex','FontSize',14); 


%% Recheck that rotation indeed leads to 0 \Delta' on axis
% problem came from the pressure profile, no idea why and may need to
% investigate this more



figure;plot(L.r_q,LY.deltap,'.')


eps_val = 0.2;
s0 = 2;

[L, LX] = equilSol('beta',0.1, 's0', s0, 'debug', 4,...
    'residuals_fun',@residuals_vec,'jacobian_fun', @jacobian_vec, ...
    'equation_of_state', @isotropic_rotating,'mach20',0,'Nb',1, 'Ns',2);
LX.eps_val = eps_val;


%LX.Sbc(1) = 0;
%LX.S2bc = 0;
%LX.Sbc(2) = 0;
%LX.S3bc = 0;
LY = equilY(L, LX);

[L2, LX2] = equilSol('beta',0.1, 's0', s0, 'debug', 4,...
    'residuals_fun',@residuals_rotation,'jacobian_fun', @jacobian_rotation, ...
    'equation_of_state', @isotropic_rotating,'mach20',0,'Nb',1, 'Ns',2);
LX2.eps_val = eps_val;


%LX.Sbc(1) = 0;
%LX.S2bc = 0;
%LX.Sbc(2) = 0;
%LX.S3bc = 0;
LY2 = equilY(L2, LX2);


figure;
subplot(4,1,1);
hold on;
plot(L.r_q, LY.gavg-LY2.gavg, '-','LineWidth',2)
ylabel('$(\langle RB_\phi\rangle-R_0B_0)/R_0B_0$','Interpreter','latex','FontSize',14);
grid on;
subplot(4,1,2);
hold on;
plot(LY.r_plt, LY.delta-LY2.delta, '.')
grid on;
subplot(4,1,3);
hold on;
plot(L.r_q, LY.deltap-LY2.deltap, '.')
grid on;
subplot(4,1,4);
hold on;
plot(LY.r_plt, squeeze(LY.S(:,:,1))-squeeze(LY2.S(:,:,1)), '.')
grid on;


%% Actual rotation test

eps_val = 0.3;
s0 = 2;

[L, LX] = equilSol('beta',0.2, 's0', s0, 'debug', 4,...
    'residuals_fun',@residuals_vec,'jacobian_fun', @jacobian_vec, ...
    'equation_of_state', @parallel_rotating,'mach20',1,'Nb',1, 'Ns',2);
LX.eps_val = eps_val;


%LX.Sbc(1) = 0;
%LX.S2bc = 0;
%LX.Sbc(2) = 0;
%LX.S3bc = 0;
LY = equilY(L, LX);

figure;
subplot(4,1,1);
hold on;
plot(L.r_q, LY.gavg-1, '-','LineWidth',2)
ylabel('$(\langle RB_\phi\rangle-R_0B_0)/R_0B_0$','Interpreter','latex','FontSize',14);
grid on;
subplot(4,1,2);
hold on;
plot(LY.r_fine, LY.delta_fine, '-','LineWidth',2)
plot(LY.r_plt, LY.delta, '.')
grid on;
subplot(4,1,3);
hold on;
plot(LY.r_fine, LY.deltap_fine, '-','LineWidth',2)
plot(L.r_q, LY.deltap, '.')
grid on;
subplot(4,1,4);
hold on;
plot(LY.r_fine, LY.S2_fine, '-','LineWidth',2)
plot(LY.r_plt, squeeze(LY.S(:,:,1)), '.')
grid on;

figure; axis equal;hold on;
contourf(LY.RR, LY.ZZ, LY.betapar)
%for i=1:10:numel(L2.omega)
%    plot(LY2.RR(:,i), LY2.ZZ(:,i), 'w')
%end
for i=1:20:numel(LY.r_plt)
    plot(LY.RR(i,:), LY.ZZ(i,:), 'w')
end

%% Legacy

% Actually testing all the git stuff too
% Testing the rotation


% 
% %% TODO
% % Understand what is going on with rotation (unbalanced r pressure like in anisotropy)
% figure;hold on;
% plot(L2.r_q, LY2.t2, '.')
% plot(LY2.r_fine, LY2.t2_fine)
% 
% 
% figure; 
% subplot(3,1,1);
% contourf(LY.RR, LY.ZZ, LY.BB)
% colorbar;
% axis equal;
% subplot(3,1,2);
% contourf(LY2.RR, LY2.ZZ, LY2.BB)
% colorbar;
% axis equal;
% subplot(3,1,3);
% contourf(LY3.RR, LY3.ZZ, LY3.BB)
% colorbar;
% axis equal;
% 
% figure;hold on;
% plot(LY.r_plt, squeeze(LY.Bs(:,:,1)))
% plot(LY2.r_plt, squeeze(LY2.Bs(:,:,1)))
% plot(LY3.r_plt, squeeze(LY3.Bs(:,:,1)))
% legend({'isotropic', '$p_\perp=0$', '$p=0$'},'box','off', 'Interpreter', 'latex',...
%     'Fontsize',14, 'location','northwest')
% grid on;

% 

