%%
% Rewritten and reordered version
eps1 = 0.1;     % eps for first (left) panel
eps2 = 0.3;    % eps for second (right) panel
s0 = 2;
beta = 0.3;
mach20 = 0.5;

s2bc = -0.5;
s3bc = 0.1;

% 1) pressureless
% 2) isotropic rotating
% 3) parallel rotating
% 4) parallel, no rotation

% 1) pressureless
[L1, LX1] = equilSol('beta',0, 's0', s0, 'debug', 4, ...
    'residuals_fun',@residuals_hyb,'jacobian_fun', @jacobian_hyb, ...
    'equation_of_state', @isotropic,'Ns', 4,'Nb',10, 'nk',30);
LX1.eps_val = eps1;
LX1.Sbc(1) = s2bc;
LX1.S2bc = s2bc;
LX1.Sbc(2) = s3bc;
LX1.S3bc = s3bc;
LY1 = equilY(L1, LX1);

% Also compute with eps2 for right panel
LX1b = LX1; LX1b.eps_val = eps2;
LY1b = equilY(L1, LX1b);

% 2) isotropic rotating (original L2)
[L2, LX2] = equilSol('beta',beta, 's0', s0, 'debug', 4, ...
    'residuals_fun',@residuals_hyb,'jacobian_fun', @jacobian_hyb, ...
    'equation_of_state', @isotropic_rotating,'Ns', 3,'mach20',mach20,'Nb',10, 'nk',30);
LX2.eps_val = eps1;
LX2.Sbc(1) = s2bc;
LX2.S2bc = s2bc;
LX2.Sbc(2) = s3bc;
LX2.S3bc = s3bc;
LY2 = equilY(L2, LX2);

LX2b = LX2; LX2b.eps_val = eps2;%LX2b.x = LY2.x; L2.P.hot_restart = true;
LY2b = equilY(L2, LX2b);

% 3) parallel rotating 
[L3, LX3] = equilSol('beta',beta, 's0', s0, 'debug', 4, ...
    'residuals_fun',@residuals_hyb,'jacobian_fun', @jacobian_hyb, ...
    'equation_of_state', @parallel_rotating,'mach20',mach20,'Nb',10, 'Ns',3, 'nk',30);
LX3.eps_val = eps1;
LX3.Sbc(1) = s2bc;
LX3.S2bc = s2bc;
LX3.Sbc(2) = s3bc;
LX3.S3bc = s3bc;
LY3 = equilY(L3, LX3);

LX3b = LX3; LX3b.eps_val = eps2; %LX3b.x = LY3.x; L3.P.hot_restart = true;
LY3b = equilY(L3, LX3b);

% 4) parallel, no rotation 
[L4, LX4] = equilSol('beta',beta, 's0', s0, 'debug', 4, ...
    'residuals_fun',@residuals_hyb,'jacobian_fun', @jacobian_hyb, ...
    'equation_of_state', @parallel_rotating,'mach20',0,'Nb',10, 'Ns',3, 'nk',30);
LX4.eps_val = eps1;
LX4.Sbc(1) = s2bc;
LX4.S2bc = s2bc;
LX4.Sbc(2) = s3bc;
LX4.S3bc = s3bc;
LY4 = equilY(L4, LX4);

LX4b = LX4; LX4b.eps_val = eps2; %LX4b.x = LY4.x; L4.P.hot_restart = true;
LY4b = equilY(L4, LX4b);

% --- Plotting: two panels (1x2)
figure;
t = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% Left panel: original style (with fine lines), eps = eps1
ax1 = nexttile;
hold(ax1,'on');
p1 = plot(L1.r_q, (squeeze(LY1.S(2:end,1,1)) - LY1.S2_ana)/eps1, '.');
plot(LY1.r_fine, LY1.S2_1_fine, '-', 'Color', p1.Color);

p2 = plot(L2.r_q, (squeeze(LY2.S(2:end,1,1)) - LY2.S2_ana)/eps1, '.');
plot(LY2.r_fine, LY2.S2_1_fine, '-', 'Color', p2.Color);

p3 = plot(L3.r_q, (squeeze(LY3.S(2:end,1,1)) - LY3.S2_ana)/eps1, '.');
plot(LY3.r_fine, LY3.S2_1_fine, '-', 'Color', p3.Color);

p4 = plot(L4.r_q, (squeeze(LY4.S(2:end,1,1)) - LY4.S2_ana)/eps1, '.');
plot(LY4.r_fine, LY4.S2_1_fine, '-', 'Color', p4.Color);

ylabel('$(\hat S_2-\hat S_2^{LO})/\epsilon$','Interpreter','latex','FontSize',14);
xlabel('$\hat r$','Interpreter','latex','FontSize',14);
grid on;


title(ax1, sprintf('$\\epsilon = %.2f$', eps1), 'Interpreter','latex');

% Right panel: same quantities but eps = eps2, no fine-line plots
ax2 = nexttile;
hold(ax2,'on');
q1 = plot(L1.r_q, (squeeze(LY1b.S(2:end,1,1)) - LY1b.S2_ana)/eps2, '.');
q2 = plot(L2.r_q, (squeeze(LY2b.S(2:end,1,1)) - LY2b.S2_ana)/eps2, '.');
q3 = plot(L3.r_q, (squeeze(LY3b.S(2:end,1,1)) - LY3b.S2_ana)/eps2, '.');
q4 = plot(L4.r_q, (squeeze(LY4b.S(2:end,1,1)) - LY4b.S2_ana)/eps2, '.');
plot(LY1.r_fine, LY1b.S2_1_fine, '-', 'Color', p1.Color);
plot(LY2.r_fine, LY2b.S2_1_fine, '-', 'Color', p2.Color);
plot(LY3.r_fine, LY3b.S2_1_fine, '-', 'Color', p3.Color);
plot(LY4.r_fine, LY4b.S2_1_fine, '-', 'Color', p4.Color);
% use same marker colors as left panel for consistency
set(q1, 'Color', p1.Color);
set(q2, 'Color', p2.Color);
set(q3, 'Color', p3.Color);
set(q4, 'Color', p4.Color);

xlabel(ax2,'$\hat r$','Interpreter','latex','FontSize',14);
title(ax2, sprintf('$\\epsilon = %.2f$', eps2), 'Interpreter','latex');
grid(ax2,'on');
leg = legend([q1 q2 q3 q4], ...
    {'$p=0$', 'isotropic $\mathcal M^2=1$', '$p_\perp=0, \mathcal M^2=1$', '$p_\perp=0, \mathcal M^2=0$'}, ...
    'Interpreter','latex','Box','off','FontSize',12, 'Location', 'SouthEast', 'NumColumns',2);
% Optionally align y-limits between panels for direct comparison
yl = ylim(ax1);
ylim(ax2, yl);

% tidy up
hold(ax1,'off');
hold(ax2,'off');




% figure;hold on;
% for i = 1:20:size(LY2.RR,1)
%     plot(LY1b.RR(i,:), LY1b.ZZ(i,:), 'k')
% end
% axis equal