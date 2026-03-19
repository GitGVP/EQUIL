
beta0=0.427160000000000;

%beta0 = linspace(0., -1, 15)*4;

% Match inputs
q0 = 0.7795296521664806;
q1 = 2.9639395354123486;
s0 = 2 * (q1/q0 -1);
Sbc = [-0.35 0.06, 0];
eps_val = 0.1;
[L, LX] = equilSol('debug',4, 'Nb', 1, 'q0', q0,'s0' , s0, ...
    'beta',beta0, 'om_pts', 300, 'do_SFL',true);

LX.eps_val = eps_val;LX.Sbc=Sbc;


LY =equilY(L, LX);

% increase beta
for betas = [L.P.beta]
    L.P.beta = betas;
    LX = equilX(L);
    LX.eps_val = eps_val;LX.Sbc=Sbc;
    beta_interp = L.P.beta * (1 - LY.psiN_q);
    ord_b = 0:2:14;
    cb = (bsxfun(@power, L.r_q(:), ord_b) \ beta_interp(:)).';     % 1xK
    beta_fit  = @(rr) reshape(    bsxfun(@power, rr(:), ord_b)   * cb.',          size(rr));
    betap_fit = @(rr) reshape(    bsxfun(@power, rr(:), ord_b-1) * (ord_b.*cb).', size(rr));

    LX.kinetic_profiles.beta  = beta_fit;
    LX.kinetic_profiles.betap = betap_fit;
    LX.x = LY.x;
    LY = equilY(L,LX);
end

filename = "equi.h5";
to_venus(LX, LY, filename)




%% Example of running EQUIL and creating the HDF5 file used in VENUS-MHD
%beta0_vmec=0.425
beta1=linspace(0.005, 0.0903, 4);
beta2=linspace(0.1358, 0.5, 6);
beta0=[beta1 beta2];
%beta0 = linspace(0., -1, 15)*4;

% Match inputs
q0 = 0.7795296521664806;
q1 = 2.9639395354123486;
s0 = 2 * (q1/q0 -1);

Sbc = [-0.35 0.06, 0];
eps_val = 0.26;
[L, LX] = equilSol('debug',4, 'Nb', 1, 'q0', q0,'s0' , s0, ...
    'beta',0.3, 'om_pts', 300, 'do_SFL',true);

LX.eps_val = eps_val;LX.Sbc=Sbc;


LY =equilY(L, LX);
    
% increase beta
for i = 1:numel(beta0)
    L.P.beta = beta0(i);
    LX = equilX(L);
    LX.eps_val = eps_val;LX.Sbc=Sbc;
    beta_interp = L.P.beta * (1 - LY.psiN_q);
    ord_b = 0:2:14;
    cb = (bsxfun(@power, L.r_q(:), ord_b) \ beta_interp(:)).';     % 1xK
    beta_fit  = @(rr) reshape(    bsxfun(@power, rr(:), ord_b)   * cb.',          size(rr));
    betap_fit = @(rr) reshape(    bsxfun(@power, rr(:), ord_b-1) * (ord_b.*cb).', size(rr));

    LX.kinetic_profiles.beta  = beta_fit;
    LX.kinetic_profiles.betap = betap_fit;
    LX.x = LY.x;
    LY = equilY(L,LX);
    filename = sprintf("eqbeta%d.h5", i);
    to_venus(LX, LY, filename)
end


%%

% Example of plotting the flux surfaces
figure;axis equal;hold on;
contour(LY.RR, LY.ZZ, ...
        [0;LY.psiN;1].*ones(size(LY.RR)), linspace(0,1,11), ...
        'LineWidth',2, 'ShowText',1);
for i=floor(linspace(1,numel(LY.omega_plt),11)); plot(LY.RR(:,i), LY.ZZ(:,i),'--','Color', [0,0,0,0.5],'LineWidth',0.5); end
xlabel('$R/R_0$', 'Interpreter','latex', 'FontSize',14)
ylabel('$Z/R_0$', 'Interpreter','latex', 'FontSize',14)
title('$\psi_N$', 'Interpreter','latex', 'FontSize',14)




%% Scan epsilon
% First rerun just to compute the Next to Leading Order elongation
LX.x = LY.x;L.P.do_shift_NLO = true;L.P.hot_restart=true;
LY = equilY(L, LX);

n_equilibria = 10;
scan_eps_from_equilibrium(L,LX,LY,n_equilibria)
