
%beta0=0.425;
beta0=3;
%beta0=0.35432;
%beta0=0.000005;

%beta0 = linspace(0., -1, 15)*4;

% Match inputs
q0 = 0.7795296521664806;
q1 = 1.999;
s0 = 2 * (q1/q0 -1);
Sbc = [0 0 0];
eps_val = 0.1;
[L, LX] = equilSol('debug',4, 'Nb', 1, 'q0', q0,'s0' , s0, ...
    'beta',beta0, 'om_pts', 300, 'do_SFL',true);

LX.eps_val = eps_val;LX.Sbc=Sbc;



LY =equilY(L, LX);

r=LY.r_fine;
[~, I] = min(sqrt((LX.qfun(r)-1).^2));
[~, J] = min(sqrt((LX.qfun(r)-2).^2));
[X2_i, X2_ip] = solve_Xi_eq(r(1:I), 0, 1, LX.qfun, LX.qpfun);
[X2_e, X2_ep] = solve_Xi_eq(r(I:J), 1, 0, LX.qfun, LX.qpfun);


%%
% increase beta
for betas = [L.P.beta]
    L.P.beta = beta0;
    LX = equilX(L);
    LX.eps_val = eps_val;LX.Sbc=Sbc;
%     beta_interp = L.P.beta * (1 - LY.psiN_q);
%     ord_b = 0:2:14;
%     cb = (bsxfun(@power, L.r_q(:), ord_b) \ beta_interp(:)).';     % 1xK
%     beta_fit  = @(rr) reshape(    bsxfun(@power, rr(:), ord_b)   * cb.',          size(rr));
%     betap_fit = @(rr) reshape(    bsxfun(@power, rr(:), ord_b-1) * (ord_b.*cb).', size(rr));

%     LX.kinetic_profiles.beta  = beta_fit;
%     LX.kinetic_profiles.betap = betap_fit;
    LX.x = LY.x;
    LY = equilY(L,LX);
end
%%
filename = "equileps01.h5";
to_venus(LX, LY, filename, X2_i, X2_ip, X2_e, X2_ep)



%% Example of running EQUIL and creating the HDF5 file used in VENUS-MHD

beta0=0.9;
values_eps=[0.01 0.03 0.05 0.08 0.1 0.15 0.2 0.25 0.3 0.35];

% Match inputs
q0 = 0.7795296521664806;
q1 = 2.9639395354123486;
s0 = 2 * (q1/q0 -1);

%Sbc = [-0.35 0.06 0];
Sbc = [0 0 0];
eps_val = 0.3;

[L, LX] = equilSol('debug',4, 'Nb', 1, 'q0', q0,'s0' , s0, ...
    'beta',beta0, 'om_pts', 300, 'do_SFL',true);

LX.eps_val = eps_val;LX.Sbc=Sbc;


LY =equilY(L, LX);

r=LY.r_fine;
[~, I] = min(sqrt((LX.qfun(r)-1).^2));
[~, J] = min(sqrt((LX.qfun(r)-2).^2));
[X2_i, X2_ip] = solve_Xi_eq(r(1:I), 0, 1, LX.qfun, LX.qpfun);
[X2_e, X2_ep] = solve_Xi_eq(r(I:J), 1, 0, LX.qfun, LX.qpfun);


%%
% increase beta
for i = 1:numel(values_eps)
    L.P.beta = beta0;
    LX = equilX(L);
    LX.eps_val = values_eps(i);
    LX.Sbc=Sbc;
%     beta_interp = L.P.beta * (1 - LY.psiN_q);
%     ord_b = 0:2:14;
%     cb = (bsxfun(@power, L.r_q(:), ord_b) \ beta_interp(:)).';     % 1xK
%     beta_fit  = @(rr) reshape(    bsxfun(@power, rr(:), ord_b)   * cb.',          size(rr));
%     betap_fit = @(rr) reshape(    bsxfun(@power, rr(:), ord_b-1) * (ord_b.*cb).', size(rr));
% 
%     LX.kinetic_profiles.beta  = beta_fit;
%     LX.kinetic_profiles.betap = betap_fit;
%     qfct = @(r) 4*q0*r.^2 ./ (1-(1-r.^2).^4);
%     LX.qfun = qfct;
%     LX.qpfun = @(r) 2*r.*(3*r.^4 - 8*r.^2 + 6)./(r.^6 - 4*r.^4 + 6*r.^2 - 4).^2;
%     LX.q_vec = qfct(L.r_q);
%     LX.qp_vec = LX.qpfun(L.r_q);
    LX.x = LY.x;
    LY = equilY(L,LX);
    filename = sprintf("eqeps%d.h5", i);
    to_venus(LX, LY, filename, X2_i, X2_ip, X2_e, X2_ep)
end


%%

% Example of plotting the flux surfaces
figure;axis equal;hold on;
contour(LY.RR, LY.ZZ, ...
        LY.psiN.*ones(size(LY.RR)), linspace(0,1,11), ...
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

%%
%beta0=0.425;
beta0=0.8;
%beta0=0.35432;
%beta0=0.000005;

%beta0 = linspace(0., -1, 15)*4;

% Match inputs
q0 = 0.7795296521664806;
q1 = 2.9639395354123486;
s0 = 2 * (q1/q0 -1);
Sbc = [-0.35 0.06, 0];
eps_val = 0.01;
[L, LX] = equilSol('debug',4, 'Nb', 1, 'q0', q0,'s0' , s0, ...
    'beta',beta0, 'om_pts', 300, 'do_SFL',true);

LX.eps_val = eps_val;LX.Sbc=Sbc;



LY =equilY(L, LX);

r=LY.r_fine;
[~, I] = min(sqrt((LX.qfun(r)-1).^2));
[~, J] = min(sqrt((LX.qfun(r)-2).^2));
[X2_i, X2_ip] = solve_Xi_eq(r(1:I), 0, 1, LX.qfun, LX.qpfun);
[X2_e, X2_ep] = solve_Xi_eq(r(I:J), 1, 0, LX.qfun, LX.qpfun);


%%
% increase beta
for betas = [L.P.beta]
    L.P.beta = beta0;
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
%%
filename = "equileps.h5";
to_venus(LX, LY, filename, X2_i, X2_ip, X2_e, X2_ep)

%% Example of running EQUIL and creating the HDF5 file used in VENUS-MHD
%beta0_vmec=0.425
% beta1=linspace(0.005, 0.0903, 4);
% beta2=linspace(0.1358, 0.5, 6);
% beta0=[beta1 beta2];
%beta0 = linspace(0., -1, 15)*4;
beta0=0.6;
values_eps=logspace(-2,0,10);

% Match inputs
q0 = 0.7795296521664806;
q1 = 2.9639395354123486;
s0 = 2 * (q1/q0 -1);

%Sbc = [-0.35 0.06 0];
Sbc = [0 0 0];
eps_val = 0.1;
[L, LX] = equilSol('debug',4, 'Nb', 1, 'q0', q0,'s0' , s0, ...
    'beta',0.3, 'om_pts', 300, 'do_SFL',true);

LX.eps_val = eps_val;LX.Sbc=Sbc;


LY =equilY(L, LX);

r=LY.r_fine;
[~, I] = min(sqrt((LX.qfun(r)-1).^2));
[~, J] = min(sqrt((LX.qfun(r)-2).^2));
[X2_i, X2_ip] = solve_Xi_eq(r(1:I), 0, 1, LX.qfun, LX.qpfun);
[X2_e, X2_ep] = solve_Xi_eq(r(I:J), 1, 0, LX.qfun, LX.qpfun);


    
% increase beta
for i = 1:numel(beta0)
    L.P.beta = beta0(i);
    LX = equilX(L);
    LX.eps_val = eps_val;LX.Sbc=Sbc;
    %beta_interp = L.P.beta * (1 - LY.psiN_q);
%     ord_b = 0:2:14;
%     cb = (bsxfun(@power, L.r_q(:), ord_b) \ beta_interp(:)).';     % 1xK
%     beta_fit  = @(rr) reshape(    bsxfun(@power, rr(:), ord_b)   * cb.',          size(rr));
%     betap_fit = @(rr) reshape(    bsxfun(@power, rr(:), ord_b-1) * (ord_b.*cb).', size(rr));
% 
%     LX.kinetic_profiles.beta  = beta_fit;
%     LX.kinetic_profiles.betap = betap_fit;
    LX.x = LY.x;
    LY = equilY(L,LX);
    filename = sprintf("eqbeta%d.h5", i);
    to_venus(LX, LY, filename, X2_i, X2_ip, X2_e, X2_ep)
end

