
beta0=[1.2, 1.4, 1.5, 1.8, 2.0];
values_eps = [0.08, 0.135, 0.19, 0.245, 0.3];

% Match inputs
q0 = 0.7795296521664806;
q1 = 1.999;
s0 = 2 * (q1/q0 -1);
S2bc=[0.55, 0.49, 0.43, 0.37, 0.31]; %dq0=0.1
% S2bc=[0.55, 0.49, 0.43, 0.37, 0.31]; %dq0=0.06
Sbc = [S2bc(1) 0 0];
eps_val = values_eps(1);
[L, LX] = equilSol('debug',4, 'Nb', 1, 'q0', q0,'s0' , s0, ...
    'beta', 0, 'om_pts', 300, 'do_SFL',true, 'residuals_fun', @residuals_iso_static, 'jacobian_fun', @jacobian_iso_static);


deltaq0 = [0.03, 0.06, 0.1, 0.15, 0.2, 0.25];
dq0 = deltaq0(2);
qfct = @(r) (4*(-1+dq0))./(-4 + 6*r.^2 -4*r.^4 + r.^6);
LX.qfun = qfct;
LX.qpfun = @(r) (-8*r.*(6 - 8*r.^2 + 3*r.^4)*(-1 + dq0))./(-4 + 6*r.^2 - 4*r.^4 +r.^6).^2;

LX.eps_val = eps_val;LX.Sbc=Sbc;
LY =equilY(L, LX);

%%
L.P.beta = beta0(1);
% L.P.beta = 0.7;
LX = equilX(L);

qfct = @(r) (4*(-1+dq0))./(-4 + 6*r.^2 -4*r.^4 + r.^6);
LX.qfun = qfct;
LX.qpfun = @(r) (-8*r.*(6 - 8*r.^2 + 3*r.^4)*(-1 + dq0))./(-4 + 6*r.^2 - 4*r.^4 +r.^6).^2;


Sbc = [S2bc(1) 0 0];
LX.eps_val = eps_val; LX.Sbc=Sbc; LX.x = LY.x;
LY = equilY(L, LX);

%%
r=L.r_q;
[~, I] = min(sqrt((LX.qfun(r)-1).^2));
[~, J] = min(sqrt((LX.qfun(r)-2).^2));
V2fun = @(r) interp1(L.r_q, LY.V2, r, 'spline');
[X2_i, X2_ip] = solve_Xi_eq(r(1:I), 0, 1, LX.qfun, LX.qpfun, 2, false, V2fun);
if J == numel(r)
    neumann = false;
else
    neumann = true;
end
[X2_e, X2_ep] = solve_Xi_eq(r(I:J), 1, 0, LX.qfun, LX.qpfun, 2, neumann);

xijump = -2 * (LY.deltap_ana(I) * LX.qpfun(r(I)) + 2 * LX.kinetic_profiles.betap(r(I))) / r(I);

[xi2in, xi2out, xi2inp] = deal(zeros(size(L.r_q)));
xi2in(1:I) = xijump / (X2_ep(1) - X2_ip(end)) * X2_i;
xi2out(I:J) =  xijump / (X2_ep(1) - X2_ip(end)) * X2_e;



xi2inp(1:I) = xijump / (X2_ep(1) - X2_ip(end)) * X2_ip;
betapp = -2*L.P.beta;
% qpp = L.P.q0 * L.P.s0;
qpp = (8*(-1 + dq0)*(24 + 12*r.^2 - 156*r.^4 + 208*r.^6 - 108*r.^8 + 21*r.^10))./(-4 + 6*r.^2 - 4*r.^4 + r.^6).^3;


d2coeff = (-1 + 2./LX.qfun(L.r_q)).^2.*L.r_q.^2;
d1coeff = -(((2 - LX.qfun(L.r_q)).*L.r_q.*(-3.*2.*LX.qfun(L.r_q) + ...
3.*LX.qfun(L.r_q).^2 + 2.*2.*LX.qpfun(L.r_q).*L.r_q))./LX.qfun(L.r_q).^3);
d0coeff = (1 - 4).*(-1 + 2./LX.qfun(L.r_q)).^2;
xi2inpp = -(d1coeff .* xi2inp + d0coeff .* xi2inp + LY.V2) ./ d2coeff;

dW2 = trapz(L.r_q,(-((pi.*L.r_q.*(-2.*LY.deltap_ana.*LX.qpfun(L.r_q).^2.*L.r_q.^2.*xi2in + ...
2.*LX.qfun(L.r_q).*L.r_q.^2.*(LY.deltap_ana.*qpp.*xi2in + LX.qpfun(L.r_q).*(-2.*xi2in + ...
3.*LY.deltap_ana.*xi2inp)) + 4.*LX.qfun(L.r_q).^5.*(betapp.*L.r_q.*xi2in + ...
LX.kinetic_profiles.betap(L.r_q).*(xi2in - 2.*LX.qpfun(L.r_q).*L.r_q.*xi2in + 4.*L.r_q.*xi2inp)) - ...
LX.qfun(L.r_q).^6.*(2.*betapp.*L.r_q.*xi2in + LX.kinetic_profiles.betap(L.r_q).*(xi2in + ...
5.*L.r_q.*xi2inp)) + LX.qfun(L.r_q).^3.*(-6.*L.r_q.*(3.*xi2in + ...
2.*L.r_q.*xi2inp) + LX.qpfun(L.r_q).*L.r_q.*(2.*L.r_q.*xi2in + LY.deltap_ana.*(-3.*xi2in ...
+ 5.*L.r_q.*xi2inp)) + LY.deltap_ana.*(-9.*xi2in + 2.*qpp.*L.r_q.^2.*xi2in + ...
9.*L.r_q.*xi2inp - 3.*L.r_q.^2.*xi2inpp)) + ...
LX.qfun(L.r_q).^2.*(2.*LY.deltap_ana.*LX.qpfun(L.r_q).^2.*L.r_q.^2.*xi2in + ...
4.*L.r_q.*(3.*xi2in + 2.*L.r_q.*xi2inp) + LX.qpfun(L.r_q).*L.r_q.*(2.*L.r_q.*xi2in ...
+ LY.deltap_ana.*(3.*xi2in - 11.*L.r_q.*xi2inp)) + LY.deltap_ana.*(6.*xi2in - ...
4.*qpp.*L.r_q.^2.*xi2in - 6.*L.r_q.*xi2inp + 2.*L.r_q.^2.*xi2inpp)) + ...
LX.qfun(L.r_q).^4.*(LX.kinetic_profiles.betap(L.r_q).*(-5.*xi2in + 12.*LX.qpfun(L.r_q).*L.r_q.*xi2in - ...
11.*L.r_q.*xi2inp) + L.r_q.*((6 + betapp).*xi2in + 4.*L.r_q.*xi2inp) ...
+ LY.deltap_ana.*(3.*xi2in + L.r_q.*(-3.*xi2inp + ...
L.r_q.*xi2inpp)))))./LX.qfun(L.r_q).^4)));

growth_rate = - (LY.dW0 + dW2) / (r(I)^3 * LX.qpfun(r(I)) * sqrt(3));

[X2_i, X2_ip] = solve_Xi_eq(r(1:I), 0, 1, LX.qfun, LX.qpfun, 2);
if J == numel(r)
    neumann = false;
else
    neumann = true;
end
[X2_e, X2_ep] = solve_Xi_eq(r(I:J), 1, 0, LX.qfun, LX.qpfun, 2, neumann);


%%
figure;hold on;
plot(r, xi2in, '.')
plot(r, xi2out, '.')


%%

k=2;
if k ==2
    Sk = LY.S2_ana;
    Skp = interp1(LY.r_fine, LY.S2p_fine, L.r_q, 'spline');
elseif k == 3
    Sk = LY.S3_ana;
    Skp = interp1(LY.r_fine, LY.S3p_fine, L.r_q, 'spline');
end

% qpp = L.P.q0 * L.P.s0;
qpp = (8*(-1 + dq0)*(24 + 12*r.^2 - 156*r.^4 + 208*r.^6 - 108*r.^8 + 21*r.^10))./(-4 + 6*r.^2 - 4*r.^4 + r.^6).^3;

Vkp1 = -(((1 + k).*(1 + k - LX.qfun(L.r_q)).*((1 + k).*LX.qpfun(L.r_q).^2.*L.r_q.^2.*Skp ...
+ LX.qfun(L.r_q).*((-1 + k.^2).*LX.qpfun(L.r_q).*Sk + (LX.qpfun(L.r_q).^2 - (1 + ...
k).*qpp).*L.r_q.^2.*Skp) + LX.qfun(L.r_q).^2.*((-1 + k.^2).*LX.qpfun(L.r_q).*Sk + ...
L.r_q.*(-((2 + k).*LX.qpfun(L.r_q)) + qpp.*L.r_q).*Skp)))./(k.*LX.qfun(L.r_q).^4));


r=L.r_q;
[~, I] = min(sqrt((LX.qfun(r)-1).^2));
[~, J] = min(sqrt((LX.qfun(r)-(k+1)).^2));
Vkp1fun = @(r) interp1(L.r_q, Vkp1, r, 'spline');
[Xkp1_i, Xkp1_ip] = solve_Xi_eq(r(1:I), 0, 1, LX.qfun, LX.qpfun, k+1, false, Vkp1fun);
if J == numel(r)
    neumann = false;
else
    neumann = true;
end
[Xkp1_e, Xkp1_ep] = solve_Xi_eq(r(I:J), 1, 0, LX.qfun, LX.qpfun, k+1, neumann);

xijump =  (1+k)/k * LX.qpfun(r(I)) * Skp(I);

[xikp1in, xikp1out, xikp1inp] = deal(zeros(size(L.r_q)));
xikp1in(1:I) = xijump / (Xkp1_ep(1) - Xkp1_ip(end)) * Xkp1_i;
xikp1out(I:J) =  xijump / (Xkp1_ep(1) - Xkp1_ip(end)) * Xkp1_e;

%%
figure;hold on;
plot(r, xikp1in, '.')
plot(r, xikp1out, '.')

%% Y0_beta2_1.h5
%eps_beta2S2bc09_11.h5
%S2bc_beta2eps01_1.h5
%beta_S2bc0S3bc0eps01qnew_3
%dq0_S2bc0S3bc0eps01beta09_4

filename = "eps_S2bc0S3bc0beta12dq006_5.h5";
to_venus(LX, LY, filename, X2_i, X2_ip, X2_e, X2_ep, xi2in, xi2out, growth_rate, xikp1in, xikp1out)

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




%% Example of running EQUIL and creating the HDF5 file used in VENUS-MHD

beta0=0.9;
values_eps=[0.01 0.03 0.05 0.08 0.1 0.15 0.2 0.25 0.3 0.35];

% Match inputs
q0 = 0.7795296521664806;
q1 = 1.999;
s0 = 2 * (q1/q0 -1);

%Sbc = [-0.35 0.06 0];
Sbc = [0 0 0];
eps_val = 0.3;

[L, LX] = equilSol('debug',4, 'Nb', 1, 'q0', q0,'s0' , s0, ...
    'beta',beta0, 'om_pts', 300, 'do_SFL',true);

LX.eps_val = eps_val;LX.Sbc=Sbc;


LY =equilY(L, LX);
%%
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

