%% Variational and standard counterparts of examples/simple_case.m

beta_final = 20;
epsilon = 0.4;
Ns = 3;
Sbc = zeros(Ns,1);

%% Variational solve: fixed-epsilon beta continuation is internal
[L,LX] = equilVariationalSol('beta',beta_final,'debug',4,'Ns',Ns);
LX.eps_val = epsilon;
LX.Sbc = Sbc;
LY = equilVariationalY(L,LX);

assert(LY.isconverged);
assert(LY.local_B_residual < 1e-10);
fprintf(['Variational: beta=%g, |R|=%.3e, ', ...
    'max|epsilon Delta''|=%.4f\n'],beta_final,norm(LY.residual), ...
    max(abs(epsilon*LY.deltap)));

%% Standard strong-form solve: manual beta continuation at fixed epsilon
beta_stages = unique([0,0.1,1,3,10,beta_final],'stable');
LYold = [];
for stage = 1:numel(beta_stages)
    [Lold,LXold] = equilSol('beta',beta_stages(stage),'debug',0, ...
        'Ns',Ns,'Nb',1,'m',8,'nq',6,'spline_p',4,'om_pts',128, ...
        'residuals_fun',@residuals_iso_static, ...
        'jacobian_fun',@jacobian_iso_static, ...
        'equation_of_state',@isotropic,'do_ana',false,'NLtol',1e-11);
    LXold.eps_val = epsilon;
    LXold.Sbc = Sbc;
    if ~isempty(LYold)
        LXold.x = LYold.x;
    end
    LYold = equilY(Lold,LXold);
    assert(LYold.isconverged, ...
        'Standard solve failed at beta=%g.',beta_stages(stage));
end
fprintf(['Standard:    beta=%g, |R|=%.3e, ', ...
    'max|epsilon Delta''|=%.4f\n'],beta_final,LYold.res_norms(end), ...
    max(abs(epsilon*LYold.deltap)));

fields = {'t2','delta','P'};
labels = {'t2';'Delta';'P'};
absolute_error = zeros(5,1);
relative_error = zeros(5,1);
for k = 1:3
    variational = LY.(fields{k});
    standard = interp1(LYold.r_plt,LYold.(fields{k}),LY.r_plt,'pchip');
    difference = variational-standard;
    absolute_error(k) = max(abs(difference));
    relative_error(k) = norm(difference)/max(norm(standard),eps);
end
for k = 1:2
    variational = LY.S(:,:,k);
    standard = interp1(LYold.r_plt,LYold.S(:,:,k),LY.r_plt,'pchip');
    difference = variational-standard;
    absolute_error(3+k) = max(abs(difference));
    relative_error(3+k) = norm(difference)/max(norm(standard),eps);
end
comparison = table(labels,absolute_error(1:3),relative_error(1:3), ...
    'VariableNames',{'profile','max_abs_error','relative_L2_error'});
comparison = [comparison; table({'S2';'S3'},absolute_error(4:5), ...
    relative_error(4:5),'VariableNames',comparison.Properties.VariableNames)];
disp(comparison);
assert(all(relative_error(1:4) < [0.01;0.01;0.02;0.03]));
assert(absolute_error(5) < 0.015);
assert(abs(max(abs(epsilon*LY.deltap)) ...
    -max(abs(epsilon*LYold.deltap))) < 0.01);
axis_mask = LYold.r_plt > 0 & LYold.r_plt <= 0.05;
standard_axis_slope = LYold.r_plt(axis_mask) ...
    \squeeze(LYold.S(axis_mask,1,1));
fprintf(['S2''(0) from the constrained variational space: %.6f; ', ...
    'legacy S2/r fit on r<=0.05: %.6f\n'], ...
    LY.Sp(1,1,1),standard_axis_slope);


figure;
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
profile_data = {LY.delta,LYold.delta, ...
    LY.P,LYold.P, ...
    LY.S(:,:,1),LYold.S(:,:,1), ...
    LY.S(:,:,2),LYold.S(:,:,2)};
profile_labels = {'$\hat\Delta$','$\hat P$','$\hat S_2$','$\hat S_3$'};
for k = 1:4
    nexttile; hold on; grid on; box on
    plot(LY.r_plt,profile_data{2*k-1},'LineWidth',2, ...
        'DisplayName','variational');
    plot(LYold.r_plt,profile_data{2*k},'--','LineWidth',2, ...
        'DisplayName','standard');
    xlabel('$\hat r$','Interpreter','latex','FontSize',14);
    ylabel(profile_labels{k},'Interpreter','latex','FontSize',14);
    if k == 1
        legend('Location','best','Box','off','Interpreter','latex');
    end
end

%% Flux surfaces and constant-omega curves
figure; hold on; axis equal; grid on
YY = {LY,LYold};

for i = 1:2
    Y = YY{i};
    mask = Y.ZZ.*(-1)^(i+1) >= 0;
    R = Y.RR; Z = Y.ZZ; psi = Y.psiN.*ones(size(R));
    R(~mask) = NaN; Z(~mask) = NaN; psi(~mask) = NaN;

    contour(R,Z,psi,linspace(0,1,11), ...
        'LineWidth',1.5,'LineStyle',repmat('-',1,i));

    for k = unique(floor(linspace(1,numel(Y.omega_plt),11)))
        plot(R(:,k),Z(:,k),'--','Color',[0 0 0 .5],'LineWidth',.5)
    end
end

xlabel('$R/R_0$','Interpreter','latex','FontSize',12);
ylabel('$Z/R_0$','Interpreter','latex','FontSize',12);
title('Variational ($Z>0$) / standard strong form ($Z<0$)', ...
    'Interpreter','latex');


figure; hold on; axis equal; grid on
c1 = [0.12 0.20 0.38];   % midnight navy
c2 = [0.83 0.58 0.12];   % antique gold
plot(NaN,NaN,'-','Color',c1,'LineWidth',1.5);
plot(NaN,NaN,'--','Color',c2,'LineWidth',1.5);
contour(LY.RR,LY.ZZ,LY.psiN.*ones(size(LY.RR)),linspace(0,1,11), ...
    'Color',c1,'LineWidth',1.5);

contour(LYold.RR,LYold.ZZ,LYold.psiN.*ones(size(LYold.RR)),linspace(0,1,11), ...
    'Color',c2,'LineWidth',1.5,'LineStyle','--');

xlabel('$R/R_0$','Interpreter','latex','FontSize',12);
ylabel('$Z/R_0$','Interpreter','latex','FontSize',12);
legend('variational','strong form','Interpreter','latex', 'Box','off');
