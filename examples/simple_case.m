%% Simple base case, to get familiar with I/O and plotting

% Runs with all default inputs besides debug = 4 which means that 
% Newton iterations are printed.

[L, LX, LY] = equilSol('debug',4, 'Nb', 1);

figure;axis equal;hold on;
contour(LY.RR, LY.ZZ, ...
        [0;LY.psiN].*ones(size(LY.RR)), linspace(0,1,11), ...
        'LineWidth',2, 'ShowText',1);
for i=floor(linspace(1,numel(LY.omega_plt),11)); plot(LY.RR(:,i), LY.ZZ(:,i),'--','Color', [0,0,0,0.5],'LineWidth',0.5); end
xlabel('$R/R_0$', 'Interpreter','latex', 'FontSize',14)
ylabel('$Z/R_0$', 'Interpreter','latex', 'FontSize',14)
title('$\psi_N$', 'Interpreter','latex', 'FontSize',14)

maxBs = squeeze(max(abs(LY.Bs)));
figure;hold on;
%plot(LY.r_plt, LY.Bs(:,1,1)/maxBs(1) ./ LY.r_plt.^2)
for i = 2:5
    plot(LY.r_plt, LY.Bs(:,1,i)/maxBs(i) ./ LY.r_plt.^(i-1))
end



[L, LX, LY] = equilSol('debug',4,'residuals_fun', @residuals_dummy, ...
    'jacobian_fun',@jacobian_dummy);


[L, LX] = equilSol('debug', 4);
LX.Sbc(1) = -0.25; LX.S2bc = -0.4;
LY = equilY(L, LX);


L.P.hot_restart = true; LX.x = LY.x; L.P.do_shift_NLO = true;
LY = equilY(L, LX);

S2_ana = LY.S2_ana; + LX.eps_val * interp1(LY.r_fine, LY.S2_1_fine, L.r_q, 'spline');
delta_ana = LY.delta_ana + LX.eps_val * interp1(LY.r_fine, LY.delta1_ana, L.r_q, 'spline');
RR = 1 - delta_ana.*LX.eps_val.^2 + ...
    LX.eps_val.^3.*LY.P(2:end).*cos(LY.omega_plt) + ...
    LX.eps_val.*L.r_q.*cos(LY.omega_plt) + ...
    LX.eps_val.^2.* cos(LY.omega_plt) .* S2_ana + ...
    LX.eps_val.^2.* cos(2*LY.omega_plt) .* LY.S3_ana;
ZZ = LX.eps_val.^3.*LY.P(2:end).*sin(LY.omega_plt) + ...
    LX.eps_val.*L.r_q.*sin(LY.omega_plt) - ...
    LX.eps_val.^2.* sin(LY.omega_plt) .* S2_ana - ...
    LX.eps_val.^2.* sin(2*LY.omega_plt) .* LY.S3_ana;

figure;axis equal; hold on;
for i=floor(linspace(1,numel(LY.r_plt),7)); plot(LY.RR(i,1:floor(end/2)), LY.ZZ(i,1:floor(end/2)),'Color', [0 0 0.5],'LineWidth',1); end
for i=floor(linspace(1,numel(L.r_q),7)); plot(RR(i,floor(end/2):end), ZZ(i,floor(end/2):end),'--','Color', 'k','LineWidth',1); end

%for i=floor(linspace(1,numel(LY.r_plt),7)); plot(LY.RR(i,:), LY.ZZ(i,:),'Color', [0 0 0.5],'LineWidth',1); end
%for i=floor(linspace(1,numel(L.r_q),7)); plot(RR(i,:), ZZ(i,:),'--','Color', 'k','LineWidth',1); end
xlabel('$R/R_0$', 'Interpreter', 'latex', 'Fontsize',12)
ylabel('$Z/R_0$', 'Interpreter', 'latex', 'Fontsize',12) 