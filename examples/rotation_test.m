% Actually testing all the git stuff too
% Testing the rotation

eps_val = 0.35;
s0 = 5;

[L, LX] = equilSol('beta',0.5, 's0', s0, 'debug', 4,...
    'residuals_fun',@residuals_rotation,'jacobian_fun', @jacobian_rotation, ...
    'equation_of_state', @isotropic,'Ns', 2,'mach20',0,'Nb',10);
LX.eps_val = eps_val;


LX.Sbc(1) = 0;
LX.S2bc = 0;
LX.Sbc(2) = 0;
LX.S3bc = 0;
LY = equilY(L, LX);


[L2, LX2] = equilSol('beta',0.5, 's0', s0, 'debug', 4,...
    'residuals_fun',@residuals_rotation,'jacobian_fun', @jacobian_rotation, ...
    'equation_of_state', @parallel_rotating,'Ns', 2,'mach20',0,'Nb',10);
LX2.eps_val = eps_val;


LX2.Sbc(1) = 0;
LX2.S2bc = 0;
LX2.Sbc(2) = 0;
LX2.S3bc = 0;
LY2 = equilY(L2, LX2);

[L3, LX3] = equilSol('beta',0, 's0', s0, 'debug', 4,...
    'residuals_fun',@residuals_rotation,'jacobian_fun', @jacobian_rotation, ...
    'equation_of_state', @isotropic,'Ns', 2,'mach20',0,'Nb',10);
LX3.eps_val = eps_val;


LX3.Sbc(1) = 0;
LX3.S2bc = 0;
LX3.Sbc(2) = 0;
LX3.S3bc = 0;
LY3 = equilY(L3, LX3);


figure;
subplot(2,1,1);
hold on;
plot(L.r_q, LY.t2, '-','LineWidth',2)
plot(L2.r_q, LY2.t2, '.')
plot(L3.r_q, LY3.t2, '.')
legend({'isotropic', '$p_\perp=0$', '$p=0$'},'box','off', 'Interpreter', 'latex',...
    'Fontsize',14, 'location','northwest')
grid on;
subplot(2,1,2);
hold on;
plot(LY.r_plt, LY.delta, '-','LineWidth',2)
plot(LY2.r_plt, LY2.delta, '.')
plot(LY3.r_plt, LY3.delta, '.')
grid on;

figure; axis equal;hold on;
contourf(LY2.RR, LY2.ZZ, LY2.betapar)
%for i=1:10:numel(L2.omega)
%    plot(LY2.RR(:,i), LY2.ZZ(:,i), 'w')
%end
for i=1:20:numel(LY2.r_plt)
    plot(LY2.RR(i,:), LY2.ZZ(i,:), 'w')
end


%% TODO
% Understand what is going on with rotation (unbalanced r pressure like in anisotropy)
% finish to understand the t2 LO and impact on diamagnetism
% update and recheck the scans in epsilon

%% Important: t2 is really not like the pressureless case in the parallel case! It is B_\phi
%% that should be like this! I should verify it.

% perfect agreement at low eps as expected for all 3
figure;hold on;
plot(L2.r_q, LY2.t2, '.')
plot(LY2.r_fine, LY2.t2_fine)


figure; 
subplot(3,1,1);
contourf(LY.RR, LY.ZZ, LY.BB)
colorbar;
axis equal;
subplot(3,1,2);
contourf(LY2.RR, LY2.ZZ, LY2.BB)
colorbar;
axis equal;
subplot(3,1,3);
contourf(LY3.RR, LY3.ZZ, LY3.BB)
colorbar;
axis equal;

figure;hold on;
plot(LY.r_plt, squeeze(LY.Bs(:,:,1)))
plot(LY2.r_plt, squeeze(LY2.Bs(:,:,1)))
plot(LY3.r_plt, squeeze(LY3.Bs(:,:,1)))
legend({'isotropic', '$p_\perp=0$', '$p=0$'},'box','off', 'Interpreter', 'latex',...
    'Fontsize',14, 'location','northwest')
grid on;


[L, LX] = equilSol('beta',0.25, 'A0',2,'s0', s0, 'debug', 4,...
    'residuals_fun',@residuals_rotation,'jacobian_fun', @jacobian_rotation, ...
    'equation_of_state', @biMaxwellian,'Ns', 2,'mach20',0,'Nb',10);
LX.eps_val = eps_val;


LX.Sbc(1) = 0;
LX.S2bc = 0;
LX.Sbc(2) = 0;
LX.S3bc = 0;
LY = equilY(L, LX);

