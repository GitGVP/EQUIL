%% Here I want to check that we get similar results with anamak.
% try to match the pressure and TT' profiles


% Can compare shaping, q profile, pressure profiles etc.

figure;hold on;
plot(Le.r_q, LXe.q_vec, '.')
plot(L.pQ.^2, 1./LY.iqQ, 'x')
%plot(Le.r_q, 1 + 4.2 ./ (1+ (0.8 ./ Le.r_q).^5 ))
%meq_q = interp1(L.pQ(:).^2,  1./LY.iqQ, Le.r_q);
%plot(Le.r_q, meq_q, '--')
% L.pQ.^2 is LYe.PsiN

figure;hold on;
plot(LYe.psiN, LYe.kappa(2:end), '.')
plot(L.pQ.^2, LY.kappa, 'x')

figure;hold on;
plot(LYe.psiN, LYe.deltatrig(2:end), '.')
plot(L.pQ.^2, LY.delta, 'x')


% Match exactly the q profiles ?



%% fit the q profile

[L, LX, LY] = fgs('ana',1,0, 'iterq',50);
[Le, LXe] = equilSol('beta', 0, 's0',0.01, 'debug', 4);
LXe.eps_val = 0.1;


LXe.Sbc(1) = 0;
LXe.S2bc = -0;
LXe.Sbc(2) = 0;
LXe.S3bc = 0;

meq_qfun = @(r) ppval(spline(L.pQ(:).^2,  [0 1./LY.iqQ.' 0]),r);
meq_qpfun = @(r) ppval(fnder(spline(L.pQ(:).^2,  [0 1./LY.iqQ.' 0])),r);

%meq_qfun = @(r) ppval(makima(L.pQ(:).^2,1./LY.iqQ),r);
%meq_qpfun = @(r) ppval(fnder(makima(L.pQ(:).^2, 1./LY.iqQ)),r);


LXe.qfun = meq_qfun;
LXe.qpfun = meq_qpfun;
LXe.q_vec = meq_qfun(Le.r_q);
LXe.qp_vec = meq_qpfun(Le.r_q);

LYe = equilY(Le, LXe);


%%
addpath('/home/vanparys/Documents/PhD/Codes/meq')

% need to run FBT for different betas because the coil currents need to
% adapt
% , 'fbtagcon', {'Ip', 'bp', 'li'}

[L, LX] = fbt('ana',1,0, 'iterq',50, 'nr', 66, 'nz',64, 'ifield',true);
%LX = fbtxana(0,L);
bp_target = LYe.bp_liu;
%li_target = 1.1141;
%LX.li = li_target;
LX.bp = bp_target;
LX.bpD = bp_target;
%LX = fbtx(L, LX);
LX.Wk = bp_target * LX.Wk;
LY = fbtt(L, LX);

%s0 = 0.8, beta = 0.8 works, let's try less
% s0 = 1 beta = 0.6

% s0 = 4, beta = 1/16 works to eps = 0.37225 
% but actually beta is too small
[Le, LXe] = equilSol('beta',3, 's0', 0, 'debug', 4);
LXe.eps_val = 0.1;


LXe.Sbc(1) = 0;
LXe.S2bc = -0;
LXe.Sbc(2) = 0;
LXe.S3bc = 0;

LYe = equilY(Le, LXe);

Le.P.hot_restart = 'true';
Le.P.beta = 3.5;
LXe = equilX(Le);
LXe.x = LYe.x;
LYe = equilY(Le, LXe);


LXe.eps_val = 0.15;
LYe = equilY(Le, LXe);



for current_eps = [0.2 0.225 0.25 0.275 0.3] %0.325, 0.335 0.345 0.3515 0.365 0.37 0.37225]
    LXe.x = LYe.x;
    LXe.eps_val = current_eps;
    LYe = equilY(Le, LXe);

    fprintf('epsilon %.4e, li %.4e, bp %.4e\n', current_eps, LYe.li_liu, LYe.bp_liu)
end

% check if the eps is good
fprintf('vertical span of MEQ run (%.4e, %.4e)\nvertical span of equil run (%.4e, %.4e) \n',...
    min(LY.zq(:)), max(LY.zq(:)), min(LY.rA * LYe.ZZ(:)), max(LY.rA * LYe.ZZ(:)))

fprintf('horizontal span of MEQ run (%.4e, %.4e)\nhorizontal span of equil run (%.4e, %.4e) \n',...
    min(LY.rq(:)), max(LY.rq(:)), min(LY.rA * LYe.RR(:)), max(LY.rA * LYe.RR(:)))

figure; axis equal; hold on
FNfull = nan(size(L.rrx));                     
mask = L.zzx <= LY.zA;
FNfull(mask) = (LY.Fx(mask)-LY.FA)/(LY.FB-LY.FA);
contour(L.rrx, L.zzx, FNfull, linspace(0,1,11), 'LineWidth',2,'ShowText','on');
contour(LY.rA*LYe.RR(:,1:end/2),LY.zA +LY.rA*LYe.ZZ(:,1:end/2),...
    [0;LYe.psiN].*ones(size(LYe.RR(:,1:end/2))), linspace(0,1,11), 'LineWidth',2,'ShowText','on')
colorbar

figure;
plot(L.rx, (LY.Fx(17,:)-LY.FA)/(LY.FB-LY.FA), '.')
plot(L.rx, LY.Fx(17,:), '.')


figure;hold on;
plot(LYe.psiN, current_eps^2 / ( 4 * pi * 1e-7 ) * LXe.beta_vec * (LY.rBt/LY.rA).^2)
plot(L.pQ.^2, LY.PQ)


[Ltest, LXtest] = equilSol('beta', 3, 's0',0, 'debug', 4);
LXtest.eps_val = 1e-1;
LYtest = equilY(Ltest, LXtest);
fprintf('epsilon %.4e, li %.4e, bp %.4e\n', 1e-3, LYtest.li_liu, LYtest.bp_liu)


%% Check with shaping

addpath('/home/vanparys/Documents/PhD/Codes/meq')


[L, LX, LY] = fbt('ana',6,0);
% need to run FBT for different betas because the coil currents need to
% adapt
% , 'fbtagcon', {'Ip', 'bp', 'li'}

[L, LX] = fbt('ana',6,0, 'iterq',50, 'nr', 66, 'nz',64, 'ifield',true);
%LX = fbtxana(0,L);
bp_target = LYe.bp_liu;
%li_target = 1.1141;
%LX.li = li_target;
LX.bp = bp_target;
LX.bpD = bp_target;
LX.Ip = 2e5;
LX.IpD =  2e5;
%LX = fbtx(L, LX);
LX.Wk = bp_target * LX.Wk;
LY = fbtt(L, LX);


% check with shaping
% need to prepare initial guess
[Le, LXe] = equilSol('beta',1.5, 's0', 4, 'debug', 4,...
    'jacobian_fun', @jacobian_noEL_2, 'Ns', 10);
LXe.eps_val = 0.4;


LXe.Sbc(1) = -0.18;
LXe.S2bc = -0.18;
LXe.Sbc(2) = -0.23;
LXe.S3bc = -0.23;
LXe.Sbc(3) = -0.14;
LYe = equilY(Le, LXe);



% check solution against LO
figure;hold on;
plot(LYe.r_plt,LYe.delta, '.')
plot(Le.r_q, LYe.delta_ana)
plot(Le.r_q, LYe.delta_ana + LXe.eps_val * interp1(LYe.r_fine, LYe.delta1_ana, Le.r_q, 'spline'))

figure;hold on;
plot(LYe.r_plt, squeeze(LYe.S(:,:,1)))
plot(LYe.r_fine, LYe.S2_fine)
% could try to compare with NLO computed numerically.

figure;hold on;
plot(LYe.r_plt, squeeze(LYe.S(:,:,2)))
plot(LYe.r_fine, LYe.S3_fine)

figure;hold on;
for i=4:10
    plot(LYe.r_plt, squeeze(LYe.S(:,:,i)))
end


%% For the presentation
% Goal: have three plots of the flux surfaces, with various shifts.
% probably easier to make it on shot 1.

addpath('/home/vanparys/Documents/PhD/Codes/meq')
[L, LX] = fbt('ana',1,0, 'iterq',50, 'nr', 66, 'nz',64, 'ifield',true);
bp_target = 0;
LX.bp = bp_target;
LX.bpD = bp_target;
LX.Wk = bp_target * LX.Wk;
LY = fbtt(L, LX);

[Le, LXe] = equilSol('beta',0, 's0', 5, 'debug', 4,...
    'jacobian_fun', @jacobian_noEL_2, 'Ns', 2);
LXe.eps_val = 0.37225;


LXe.Sbc(1) = 0;
LXe.S2bc = 0;
LXe.Sbc(2) = 0;
LXe.S3bc = 0;
LYe = equilY(Le, LXe);


figure; axis equal; hold on
FNfull = nan(size(L.rrx));                     
mask = L.zzx <= LY.zA;
FNfull(mask) = (LY.Fx(mask)-LY.FA)/(LY.FB-LY.FA);
contour(L.rrx, L.zzx, FNfull, linspace(0,1,11), 'LineWidth',2,'ShowText','on');
contour(LY.rA*LYe.RR(:,1:end/2),LY.zA +LY.rA*LYe.ZZ(:,1:end/2),...
    [0;LYe.psiN].*ones(size(LYe.RR(:,1:end/2))), linspace(0,1,11), 'LineWidth',2,'ShowText','on')
colorbar

[Le2, LXe2] = equilSol('beta',0.5, 's0', 5, 'debug', 4,...
    'jacobian_fun', @jacobian_noEL_2, 'Ns', 2);
LXe2.eps_val = 0.37225;


LXe2.Sbc(1) = 0;
LXe2.S2bc = 0;
LXe2.Sbc(2) = 0;
LXe2.S3bc = 0;
LYe2 = equilY(Le2, LXe2);

[L2, LX2] = fbt('ana',1,0, 'iterq',50, 'nr', 66, 'nz',64, 'ifield',true);
bp_target = LYe2.bp_liu;
LX2.bp = bp_target;
LX2.bpD = bp_target;
LX2.Wk = bp_target * LX2.Wk;
LY2 = fbtt(L2, LX2);


figure; axis equal; hold on
FNfull2 = nan(size(L2.rrx));                     
mask = L2.zzx <= LY2.zA;
FNfull2(mask) = (LY2.Fx(mask)-LY2.FA)/(LY2.FB-LY2.FA);
contour(L2.rrx, L2.zzx, FNfull2, linspace(0,1,11), 'LineWidth',2,'ShowText','on');
contour(LY2.rA*LYe2.RR(:,1:end/2),LY2.zA +LY2.rA*LYe2.ZZ(:,1:end/2),...
    [0;LYe2.psiN].*ones(size(LYe2.RR(:,1:end/2))), linspace(0,1,11), 'LineWidth',2,'ShowText','on')
colorbar


% Compacted: loop over betas and plot n_beta-by-1 subplots
addpath('/home/vanparys/Documents/PhD/Codes/meq')

% user controls
betas = [0, 0.25, 0.5, 1];        % change/expand this vector to ramp beta
n_beta = numel(betas);
nx = struct('iterq',50,'nr',66,'nz',64,'ifield',true);  % fbt args

fig = figure('Name','beta ramp', 'Units','normalized', 'Position',[0.1 0.1 0.5 0.8]);
tiledlayout(fig, 1, n_beta, 'Padding', 'none', 'TileSpacing', 'none');

for i = 1:n_beta
    beta_val = betas(i);

    % 1) compute equilibrium solution for requested beta
    [Le, LXe] = equilSol('beta', beta_val, 's0', 5, 'debug', 4, ...
                         'jacobian_fun', @jacobian_noEL_2, 'Ns', 2);
    LXe.eps_val = 0.37225;
    LXe.Sbc(1) = 0; LXe.S2bc = 0; LXe.Sbc(2) = 0; LXe.S3bc = 0;
    LYe = equilY(Le, LXe);

    % 2) run forward solver and set target bp from equilibrium (or 0)
    [L, LX] = fbt('ana',1,0, 'iterq', nx.iterq, 'nr', nx.nr, 'nz', nx.nz, ...
                  'ifield', nx.ifield);
    if isfield(LYe,'bp_liu') && ~isempty(LYe.bp_liu)
        bp_target = LYe.bp_liu;
    else
        bp_target = 0;
    end
    LX.bp  = bp_target;
    LX.bpD = bp_target;
    LX.Wk  = bp_target * LX.Wk;
    LY = fbtt(L, LX);

    % 3) build normalized flux and plot in subplot i
    ax = nexttile;
    axis(ax,'equal'); hold(ax,'on');

    FNfull = nan(size(L.rrx));
    mask = L.zzx <= LY.zA;
    FNfull(mask) = (LY.Fx(mask)-LY.FA) ./ (LY.FB-LY.FA);

    % contours of flux and iso-psi from equilibrium solver
    contour(ax, L.rrx, L.zzx, FNfull, linspace(0,1,11), 'LineWidth',2, 'ShowText','on');
    contour(ax, LY.rA*LYe.RR(:,1:end/2), LY.zA + LY.rA*LYe.ZZ(:,1:end/2), ...
            [0;LYe.psiN].*ones(size(LYe.RR(:,1:end/2))), linspace(0,1,11), ...
            'LineWidth',2, 'ShowText','on');
    xlabel('$R$ [m]', 'Interpreter', 'latex', 'Fontsize',12)
    if i==1
       ylabel('$Z$ [m]', 'Interpreter', 'latex', 'Fontsize',12) 
    end
    title(ax, sprintf('$\\beta_{p,\\mathrm{LIU}} = %.1f$, $l_{i,\\mathrm{LIU}} = %.1f$, $l_{i,\\mathrm{FBT}} = %.1f$', ...
    LYe.bp_liu, LYe.li_liu, LY.li), 'Interpreter', 'latex');
end
sgtitle(sprintf('$\\psi_N$'), 'Interpreter', 'latex', 'Fontsize',16);


figure;hold on;
plot(Le.r_q, LXe.q_vec, '.')
plot(L.pQ.^2, 1./LY.iqQ, 'x')
xlabel('$\hat r, \psi_N$', 'Interpreter', 'latex', 'Fontsize',12)
ylabel('$q$', 'Interpreter', 'latex', 'Fontsize',12)


figure;hold on;
plot(LYe.r_plt, LYe.delta, '.')
plot(Le.r_q, LYe.delta_ana, 'LineWidth',2)
plot(Le.r_q, LYe.delta_ana + LXe.eps_val * interp1(LYe.r_fine, LYe.delta1_ana, Le.r_q, 'spline'), '--', 'LineWidth',2)
xlabel('$\hat r$', 'Interpreter', 'latex', 'Fontsize',12)
ylabel('$\Delta$', 'Interpreter', 'latex', 'Fontsize',12)
grid on;
legend({'$\Delta_{num}$','$\Delta_0$' ,'$\Delta_0 + \epsilon \Delta_1$'}, 'Interpreter','latex','FontSize',12,'Box','off', 'Location', 'Northwest')


