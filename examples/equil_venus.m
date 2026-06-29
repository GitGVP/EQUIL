%% q Felix

values_eps = [0.08, 0.135, 0.19, 0.245, 0.3];
S2bc=[0.55, 0.49, 0.43, 0.37, 0.31]; %dq0=0.1
% S2bc=[0.61, 0.54, 0.46, 0.39, 0.33]; %dq0=0.06
deltaq0 = [0.03, 0.06, 0.1, 0.15, 0.2, 0.25];
beta0=[0.1, 0.2, 0.3, 0.4, 0.5];

Sbc = [0 0 0];
eps_val = 0.3;
[L, LX] = equilSol('debug',4, 'Nb', 1, ...
    'beta', beta0(1), 'om_pts', 300, 'do_SFL',true, 'residuals_fun', @residuals_iso_static, 'jacobian_fun', @jacobian_iso_static, 'NLtol', 1e-7);

qfct = @(r) (4*(-1+dq0))./(-4 + 6*r.^2 -4*r.^4 + r.^6);
LX.qfun = qfct;
LX.qpfun = @(r) (-8*r.*(6 - 8*r.^2 + 3*r.^4)*(-1 + dq0))./(-4 + 6*r.^2 - 4*r.^4 +r.^6).^2;
%%
L.P.beta = beta0(1);
% L.P.beta = 0.7;
LX = equilX(L);

dq0 = 0.25;
lamb=4;
r1=0.52;
qfct = @(r) 1 - dq0*(1 - (r./r1).^lamb);
LX.qfun = qfct;
LX.qpfun = @(r) dq0*lamb*(r./r1).^(lamb-1)/r1;


LX.eps_val = eps_val;LX.Sbc=Sbc; LX.x = LY.x;
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
rs = L.r_q(I);

qpp = dq0 * lamb * (lamb-1) * (r./r1) .^(lamb-2) / r1^2;
betapp = -2*L.P.beta;

V2 = (L.r_q < rs) .*((-2.*(-2 + LX.qfun(L.r_q)).*(4.*LX.kinetic_profiles.betap(L.r_q).*LX.qfun(L.r_q).^3.*LX.qpfun(L.r_q).*L.r_q ...
+ 2.*LY.deltap_ana.*LX.qpfun(L.r_q).^2.*L.r_q.^2 + LX.qfun(L.r_q).*(LX.qpfun(L.r_q) + LY.deltap_ana.*LX.qpfun(L.r_q).^2 - ...
2.*LY.deltap_ana.*qpp).*L.r_q.^2 + LX.qfun(L.r_q).^5.*(LX.kinetic_profiles.betap(L.r_q) - ...
betapp.*L.r_q) + LX.qfun(L.r_q).^2.*L.r_q.*(LY.deltap_ana.*qpp.*L.r_q + ...
LX.qpfun(L.r_q).*(-3.*LY.deltap_ana + L.r_q)) - 2.*LX.qfun(L.r_q).^4.*(-(betapp.*L.r_q) + ...
LX.kinetic_profiles.betap(L.r_q).*(1 + 2.*LX.qpfun(L.r_q).*L.r_q))))./LX.qfun(L.r_q).^4);

V2fun = @(r) interp1(L.r_q, V2, r, 'spline');
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

% qpp = L.P.q0 * L.P.s0;
% qpp = (8*(-1 + dq0)*(24 + 12*r.^2 - 156*r.^4 + 208*r.^6 - 108*r.^8 + 21*r.^10))./(-4 + 6*r.^2 - 4*r.^4 + r.^6).^3;


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

filename = "beta_eps03lamb4r1052dq025_1.h5";
to_venus(LX, LY, filename, X2_i, X2_ip, X2_e, X2_ep, xi2in, xi2out, growth_rate, xikp1in, xikp1out)


%% q Graves 2007
beta0=[0.1, 0.2, 0.3, 0.4, 0.5];

Sbc = [0 0 0];
eps_val = 0.1;
[L, LX] = equilSol('debug',4, 'Nb', 1, ...
    'beta', beta0(1), 'om_pts', 300, 'do_SFL',true, 'residuals_fun', @residuals_iso_static, 'jacobian_fun', @jacobian_iso_static, 'NLtol', 1e-7);

dq0 = 0.1;
lamb=2;
r1=0.18;
qfct = @(r) 1 - dq0*(1 - (r./r1).^lamb);
LX.qfun = qfct;
LX.qpfun = @(r) dq0*lamb*(r./r1).^(lamb-1)/r1;

LX.eps_val = eps_val;LX.Sbc=Sbc;
LY =equilY(L, LX);

k=2;
if k ==2
    Sk = LY.S2_ana;
    Skp = interp1(LY.r_fine, LY.S2p_fine, L.r_q, 'spline');
elseif k == 3
    Sk = LY.S3_ana;
    Skp = interp1(LY.r_fine, LY.S3p_fine, L.r_q, 'spline');
end
r=L.r_q;
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

for i = 1:numel(beta0)
    L.P.beta = beta0(i);
    LX = equilX(L);
    LX.eps_val = eps_val;
    LX.Sbc=Sbc;
    LX.x = LY.x;
    LY = equilY(L,LX);

    r=L.r_q;
    [~, I] = min(sqrt((LX.qfun(r)-1).^2));
    [~, J] = min(sqrt((LX.qfun(r)-2).^2));
    rs = L.r_q(I);
    
    qpp = dq0 * lamb * (lamb-1) * (r./r1) .^(lamb-2) / r1^2;
    betapp = -2*L.P.beta;
    
    V2 = (L.r_q < rs) .*((-2.*(-2 + LX.qfun(L.r_q)).*(4.*LX.kinetic_profiles.betap(L.r_q).*LX.qfun(L.r_q).^3.*LX.qpfun(L.r_q).*L.r_q ...
    + 2.*LY.deltap_ana.*LX.qpfun(L.r_q).^2.*L.r_q.^2 + LX.qfun(L.r_q).*(LX.qpfun(L.r_q) + LY.deltap_ana.*LX.qpfun(L.r_q).^2 - ...
    2.*LY.deltap_ana.*qpp).*L.r_q.^2 + LX.qfun(L.r_q).^5.*(LX.kinetic_profiles.betap(L.r_q) - ...
    betapp.*L.r_q) + LX.qfun(L.r_q).^2.*L.r_q.*(LY.deltap_ana.*qpp.*L.r_q + ...
    LX.qpfun(L.r_q).*(-3.*LY.deltap_ana + L.r_q)) - 2.*LX.qfun(L.r_q).^4.*(-(betapp.*L.r_q) + ...
    LX.kinetic_profiles.betap(L.r_q).*(1 + 2.*LX.qpfun(L.r_q).*L.r_q))))./LX.qfun(L.r_q).^4);
    
    V2fun = @(r) interp1(L.r_q, V2, r, 'spline');
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

    filename = "beta_eps01lamb4r1042dq01_5.h5";
    to_venus(LX, LY, filename, X2_i, X2_ip, X2_e, X2_ep, xi2in, xi2out, growth_rate, xikp1in, xikp1out)
end


%% q 

values_eps = [0.01, 0.05, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4];
% S2bc=[0.05, 0.085, 0.12, 0.16, 0.193]; %for eps=0.2
% S2bc=[0.058, 0.088, 0.115, 0.142, 0.163]; %for eps=0.3
S2bc=[0.002 0.018 0.027 0.032 0.043 0.052 0.055 0.058 0.0585 0.059]; %for eps=0.3
beta0=[0.5, 0.75, 1.0, 1.25, 1.5];

q0 = 0.7795296521664806;
q1 = 1.8;
s0 = 2 * (q1/q0 -1);

Sbc = [0 0 0];
eps_val = values_eps(9);
[L, LX] = equilSol('debug',4, 'Nb', 1, 'q0', q0,'s0' , s0, 'beta', beta0(1), 'om_pts', 300, ...
    'do_SFL',true, 'residuals_fun', @residuals_iso_static, 'jacobian_fun', @jacobian_iso_static, 'NLtol', 1e-7);

LX.eps_val = eps_val; LX.Sbc=Sbc;
LY =equilY(L, LX);


%%
L.P.beta = beta0(1);
% L.P.beta = 0.7;
LX = equilX(L);

Sbc = [S2bc(1) 0 0];
LX.eps_val = eps_val; LX.Sbc=Sbc; LX.x = LY.x;
LY = equilY(L, LX);

%%
r=L.r_q;
[~, I] = min(sqrt((LX.qfun(r)-1).^2));
[~, J] = min(sqrt((LX.qfun(r)-2).^2));
rs = L.r_q(I);

qpp = L.P.q0 * L.P.s0;
betapp = -2*L.P.beta;

V2 = (L.r_q < rs) .*((-2.*(-2 + LX.qfun(L.r_q)).*(4.*LX.kinetic_profiles.betap(L.r_q).*LX.qfun(L.r_q).^3.*LX.qpfun(L.r_q).*L.r_q ...
+ 2.*LY.deltap_ana.*LX.qpfun(L.r_q).^2.*L.r_q.^2 + LX.qfun(L.r_q).*(LX.qpfun(L.r_q) + LY.deltap_ana.*LX.qpfun(L.r_q).^2 - ...
2.*LY.deltap_ana.*qpp).*L.r_q.^2 + LX.qfun(L.r_q).^5.*(LX.kinetic_profiles.betap(L.r_q) - ...
betapp.*L.r_q) + LX.qfun(L.r_q).^2.*L.r_q.*(LY.deltap_ana.*qpp.*L.r_q + ...
LX.qpfun(L.r_q).*(-3.*LY.deltap_ana + L.r_q)) - 2.*LX.qfun(L.r_q).^4.*(-(betapp.*L.r_q) + ...
LX.kinetic_profiles.betap(L.r_q).*(1 + 2.*LX.qpfun(L.r_q).*L.r_q))))./LX.qfun(L.r_q).^4);

V2fun = @(r) interp1(L.r_q, V2, r, 'spline');
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

d2coeff = (-1 + 2./LX.qfun(L.r_q)).^2.*L.r_q.^2;
d1coeff = -(((2 - LX.qfun(L.r_q)).*L.r_q.*(-3.*2.*LX.qfun(L.r_q) + ...
3.*LX.qfun(L.r_q).^2 + 2.*2.*LX.qpfun(L.r_q).*L.r_q))./LX.qfun(L.r_q).^3);
d0coeff = (1 - 4).*(-1 + 2./LX.qfun(L.r_q)).^2;
xi2inpp = -(d1coeff .* xi2inp + d0coeff .* xi2inp + V2) ./ d2coeff;

dW0 = - trapz(L.r_q, (L.r_q < rs) .*(pi.*L.r_q.*(LX.kinetic_profiles.betap(L.r_q).^2.*LX.qfun(L.r_q).^6 + LY.deltap_ana.*LX.qpfun(L.r_q).*(1 + ...
LY.deltap_ana.*LX.qpfun(L.r_q)).*L.r_q.^2 + LX.qfun(L.r_q).^3.*L.r_q.*(-((1 + ...
betapp).*LY.deltap_ana) + LX.kinetic_profiles.betap(L.r_q).*(2 - 5.*LY.deltap_ana.*LX.qpfun(L.r_q)) + L.r_q) + ...
LX.qfun(L.r_q).^5.*(-(LX.kinetic_profiles.betap(L.r_q).*L.r_q) + betapp.*LY.deltap_ana.*L.r_q) - ...
LX.qfun(L.r_q).^4.*(LX.kinetic_profiles.betap(L.r_q).^2 + betapp.*LY.deltap_ana.*L.r_q - ...
4.*LX.kinetic_profiles.betap(L.r_q).*LY.deltap_ana.*LX.qpfun(L.r_q).*L.r_q) - LX.qfun(L.r_q).*L.r_q.*(-L.r_q + ...
LY.deltap_ana.*(1 + LX.qpfun(L.r_q).*L.r_q) + LY.deltap_ana.^2.*(2.*LX.qpfun(L.r_q) + LX.qpfun(L.r_q).^2.*L.r_q - ...
qpp.*L.r_q)) - LX.qfun(L.r_q).^2.*L.r_q.*(LX.kinetic_profiles.betap(L.r_q) - (2 + ...
betapp).*LY.deltap_ana + 2.*L.r_q + LY.deltap_ana.^2.*(-2.*LX.qpfun(L.r_q) + ...
qpp.*L.r_q))))./LX.qfun(L.r_q).^3);

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

growth_rate = - (dW0 + dW2) / (r(I)^3 * LX.qpfun(r(I)) * sqrt(3));

[X2_i, X2_ip] = solve_Xi_eq(r(1:I), 0, 1, LX.qfun, LX.qpfun, 2);
if J == numel(r)
    neumann = false;
else
    neumann = true;
end
[X2_e, X2_ep] = solve_Xi_eq(r(I:J), 1, 0, LX.qfun, LX.qpfun, 2, neumann);

k=2;
if k ==2
    Sk = LY.S2_ana;
    Skp = interp1(LY.r_fine, LY.S2p_fine, L.r_q, 'spline');
elseif k == 3
    Sk = LY.S3_ana;
    Skp = interp1(LY.r_fine, LY.S3p_fine, L.r_q, 'spline');
end

qpp = L.P.q0 * L.P.s0;

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
filename = "eps_beta05S2rs0_9.h5";
to_venus(LX, LY, filename, X2_i, X2_ip, X2_e, X2_ep, xi2in, xi2out, growth_rate, xikp1in, xikp1out)


%% 

% beta0=[0.5, 0.75, 0.9, 1.00, 1.1, 1.25, 1.3, 1.4, 1.5, 1.75]; %for eps=0.28
beta0=[1.0, 1.5, 2.0, 2.5, 3.0, 3.25, 3.5, 4.0, 4.5, 5.0]; %for eps=0.1
% values_eps=[0.01 0.03 0.05 0.08 0.1 0.15 0.2 0.25 0.3 0.35];

% Match inputs
q0 = 0.7795296521664806;
q1 = 1.8;
s0 = 2 * (q1/q0 -1);

%Sbc = [-0.35 0.06 0];
Sbc = [0 0 0];
eps_val = 0.1;

[L, LX] = equilSol('debug',4, 'Nb', 1, 'q0', q0,'s0' , s0, 'beta', beta0(end), 'om_pts', 300, ...
    'do_SFL',true, 'residuals_fun', @residuals_iso_static, 'jacobian_fun', @jacobian_iso_static, 'NLtol', 1e-6);

LX.eps_val = eps_val;LX.Sbc=Sbc;


LY =equilY(L, LX);

k=2;
    if k ==2
        Sk = LY.S2_ana;
        Skp = interp1(LY.r_fine, LY.S2p_fine, L.r_q, 'spline');
    elseif k == 3
        Sk = LY.S3_ana;
        Skp = interp1(LY.r_fine, LY.S3p_fine, L.r_q, 'spline');
    end
    
    qpp = L.P.q0 * L.P.s0;
    % qpp = (8*(-1 + dq0)*(24 + 12*r.^2 - 156*r.^4 + 208*r.^6 - 108*r.^8 + 21*r.^10))./(-4 + 6*r.^2 - 4*r.^4 + r.^6).^3;
    
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

% increase beta
for i = numel(beta0):-1:1
    L.P.beta = beta0(i);
    LX = equilX(L);
    LX.eps_val = eps_val;
    LX.Sbc=Sbc;
    LX.x = LY.x;
    LY = equilY(L,LX);

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
    qpp = L.P.q0 * L.P.s0;
    % qpp = (8*(-1 + dq0)*(24 + 12*r.^2 - 156*r.^4 + 208*r.^6 - 108*r.^8 + 21*r.^10))./(-4 + 6*r.^2 - 4*r.^4 + r.^6).^3;
    
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
    
    [~, I] = min(sqrt((LX.qfun(L.r_q)-1).^2));
    [~, J] = min(sqrt((LX.qfun(L.r_q)-2).^2));
    
    %qrs = LX.qfun(L.r_q(I));
    qrs= 1;
    rs = L.r_q(I);
    srs = rs * LX.qpfun(L.r_q(I));
    betap = -2 * qrs / rs^4 * trapz(L.r_q(1:I), LX.kinetic_profiles.betap(L.r_q(1:I)) .* L.r_q(1:I).^2);
    
    dWc = -rs^2/qrs^2 * ( betap/qrs^2 + (1/rs)^4 * trapz(L.r_q(1:I), (3/qrs + 1./LX.qfun(L.r_q(1:I)) .* (1./LX.qfun(L.r_q(1:I))-1/qrs))));
    
    sigma = (1/rs)^4 *  trapz(L.r_q(1:I), L.r_q(1:I).^3 .* ( 1./LX.qfun(L.r_q(1:I)).^2 - 1) );
    
    % Compute b and c
    
    betapp = -2*L.P.beta;
    qpp = (8*(-1 + dq0)*(24 + 12*L.r_q.^2 - 156*L.r_q.^4 + 208*L.r_q.^6 - 108*L.r_q.^8 + 21*L.r_q.^10))./(-4 + 6*L.r_q.^2 - 4*L.r_q.^4 + L.r_q.^6).^3;
    
    V2 = (L.r_q < rs) .*((-2.*(-2 + LX.qfun(L.r_q)).*(4.*LX.kinetic_profiles.betap(L.r_q).*LX.qfun(L.r_q).^3.*LX.qpfun(L.r_q).*L.r_q ...
    + 2.*LY.deltap_ana.*LX.qpfun(L.r_q).^2.*L.r_q.^2 + LX.qfun(L.r_q).*(LX.qpfun(L.r_q) + LY.deltap_ana.*LX.qpfun(L.r_q).^2 - ...
    2.*LY.deltap_ana.*qpp).*L.r_q.^2 + LX.qfun(L.r_q).^5.*(LX.kinetic_profiles.betap(L.r_q) - ...
    betapp.*L.r_q) + LX.qfun(L.r_q).^2.*L.r_q.*(LY.deltap_ana.*qpp.*L.r_q + ...
    LX.qpfun(L.r_q).*(-3.*LY.deltap_ana + L.r_q)) - 2.*LX.qfun(L.r_q).^4.*(-(betapp.*L.r_q) + ...
    LX.kinetic_profiles.betap(L.r_q).*(1 + 2.*LX.qpfun(L.r_q).*L.r_q))))./LX.qfun(L.r_q).^4);
    
    V2fun = @(r) interp1(L.r_q, V2, r, 'spline');
    [X2_i, X2_ip] = solve_Xi_eq(L.r_q(1:I), 0, 1, LX.qfun, LX.qpfun, 2, false, V2fun);
    if J == numel(L.r_q)
        neumann = false;
    else
        neumann = true;
    end
    [X2_e, X2_ep] = solve_Xi_eq(L.r_q(I:J-1), 1, 0, LX.qfun, LX.qpfun, 2, neumann);
    
    b = rs * X2_ip(end) / X2_i(end);
    c = rs * X2_ep(1) / X2_e(1);
    bb = -1/2 + rs * X2_ip(end) / X2_i(end) / 4;
    cc = 1/2 + rs * X2_ep(1) / X2_e(1) / 4;
    
    dWt = rs^2/qrs^4 * ( (32 * (b-c) * sigma + 9 * (b-1)*(1-c))/(64*(b-c)) -...
        (betap + sigma) * (3*(b-1)*(c+3))/ (8*(b-c)) - ...
        (betap + sigma)^2 * (c+3)*(b+3)/(4*(b-c)));
    
    % factors eps_val^2 missing
    growth_rate_bussac = -pi /sqrt(3) /srs * dWt;

   
    filename = sprintf("beta_S2bc0S3bc0eps01qold%d.h5", i);
    to_venus(LX, LY, filename, X2_i, X2_ip, X2_e, X2_ep, xi2in, xi2out, growth_rate, growth_rate_bussac, xikp1in, xikp1out)
end

%%
% beta0=[0.5, 0.75, 0.9, 1.00, 1.1, 1.25, 1.3, 1.4, 1.5, 1.75]; %for eps=0.28
%beta0=[1.0, 1.5, 2.0, 2.5, 3.0, 3.25, 3.5, 4.0, 4.5, 5.0]; %for eps=0.1
% values_eps=[0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4];
values_eps = [0.01, 0.05, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4];

% Match inputs
q0 = 0.7795296521664806;
q1 = 1.8;
s0 = 2 * (q1/q0 -1);

%Sbc = [-0.35 0.06 0];
Sbc = [0 0 0];
eps_val = values_eps(end);
[L, LX] = equilSol('debug',4, 'Nb', 1, 'q0', q0,'s0' , s0, 'beta', 0.5, 'om_pts', 300, ...
    'do_SFL',true, 'residuals_fun', @residuals_iso_static, 'jacobian_fun', @jacobian_iso_static, 'NLtol', 1e-6);

LX.eps_val = eps_val;LX.Sbc=Sbc;


LY =equilY(L, LX);

k=2;
    if k ==2
        Sk = LY.S2_ana;
        Skp = interp1(LY.r_fine, LY.S2p_fine, L.r_q, 'spline');
    elseif k == 3
        Sk = LY.S3_ana;
        Skp = interp1(LY.r_fine, LY.S3p_fine, L.r_q, 'spline');
    end
    
    qpp = L.P.q0 * L.P.s0;
    % qpp = (8*(-1 + dq0)*(24 + 12*r.^2 - 156*r.^4 + 208*r.^6 - 108*r.^8 + 21*r.^10))./(-4 + 6*r.^2 - 4*r.^4 + r.^6).^3;
    
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

% increase beta
for i = numel(values_eps):-1:1
    L.P.beta = 0.5;
    LX = equilX(L);
    LX.eps_val = values_eps(i);
    LX.Sbc=Sbc;
    LX.x = LY.x;
    LY = equilY(L,LX);

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
    qpp = L.P.q0 * L.P.s0;
    % qpp = (8*(-1 + dq0)*(24 + 12*r.^2 - 156*r.^4 + 208*r.^6 - 108*r.^8 + 21*r.^10))./(-4 + 6*r.^2 - 4*r.^4 + r.^6).^3;
    
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
    
   
    filename = sprintf("eps_beta05S2bc0_%d.h5", i);
    to_venus(LX, LY, filename, X2_i, X2_ip, X2_e, X2_ep, xi2in, xi2out, growth_rate, xikp1in, xikp1out)
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


%% Example of running EQUIL and creating the HDF5 file used in VENUS-MHD
%beta0_vmec=0.425
% beta1=linspace(0.005, 0.0903, 4);
% beta2=linspace(0.1358, 0.5, 6);
% beta0=[beta1 beta2];
%beta0 = linspace(0., -1, 15)*4;
beta0=0.06;
values_eps=logspace(-2,0,10);

% Match inputs
q0 = 0.7795296521664806;
q1 = 2.9639395354123486;
s0 = 2 * (q1/q0 -1);

%Sbc = [-0.35 0.06 0];
Sbc = [0 0 0];
eps_val = 0.1;
[L, LX] = equilSol('debug',4, 'Nb', 1, 'q0', q0,'s0' , s0, ...
    'beta',0.06, 'om_pts', 300, 'do_SFL',true);

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

%% Bussac growth rate
% recompute from scratch.

beta0=[1.2, 1.4, 1.5, 1.8, 2.0];
deltaq0 = [0.03, 0.06, 0.1, 0.15, 0.2, 0.25];

values_eps = [0.08, 0.135, 0.19, 0.245, 0.3];
S2bc=[0.55, 0.49, 0.43, 0.37, 0.31]; %dq0=0.1
%Sbc = [S2bc(1) 0 0];

gr_vec = zeros(6,1);
% dq0 = 0.1, beta0 = 1.2, 
for ii =1:6
    dq0 = deltaq0(ii);
    Sbc = zeros(3,1);
    eps_val = values_eps(end);
    [L, LX] = equilSol('debug',4, 'Nb', 1, ...
        'beta', beta0(1), 'om_pts', 300, 'do_SFL',true, 'residuals_fun', @residuals_iso_static, 'jacobian_fun', @jacobian_iso_static);
    
    qfct = @(r) (4*(-1+dq0))./(-4 + 6*r.^2 -4*r.^4 + r.^6);
    LX.qfun = qfct;
    LX.qpfun = @(r) (-8*r.*(6 - 8*r.^2 + 3*r.^4)*(-1 + dq0))./(-4 + 6*r.^2 - 4*r.^4 +r.^6).^2;
    
    LX.eps_val = eps_val;LX.Sbc=Sbc;
    LY =equilY(L, LX);
    
    
    
    % Bussac growth rate computation
    % might want to enforce qrs = 1.
    % in fact, might want to really make this whole thing more accurate with
    % respect to the proper position of the q = 1 surface.
    
    [~, I] = min(sqrt((LX.qfun(L.r_q)-1).^2));
    [~, J] = min(sqrt((LX.qfun(L.r_q)-2).^2));
    
    %qrs = LX.qfun(L.r_q(I));
    qrs= 1;
    rs = L.r_q(I);
    srs = rs * LX.qpfun(L.r_q(I));
    betap = -2 * qrs / rs^4 * trapz(L.r_q(1:I), LX.kinetic_profiles.betap(L.r_q(1:I)) .* L.r_q(1:I).^2);
%     
%     dWc = -rs^2/qrs^2 * ( betap/qrs^2 + (1/rs)^4 * trapz(L.r_q(1:I), (3/qrs + 1./LX.qfun(L.r_q(1:I)) .* (1./LX.qfun(L.r_q(1:I))-1/qrs))));
    
    sigma = (1/rs)^4 *  trapz(L.r_q(1:I), L.r_q(1:I).^3 .* ( qrs^2./LX.qfun(L.r_q(1:I)).^2 - 1) );
    
    % Compute b and c
    
    betapp = -2*L.P.beta;
    qpp = (8*(-1 + dq0)*(24 + 12*L.r_q.^2 - 156*L.r_q.^4 + 208*L.r_q.^6 - 108*L.r_q.^8 + 21*L.r_q.^10))./(-4 + 6*L.r_q.^2 - 4*L.r_q.^4 + L.r_q.^6).^3;
    
    V2 = (L.r_q < rs) .*((-2.*(-2 + LX.qfun(L.r_q)).*(4.*LX.kinetic_profiles.betap(L.r_q).*LX.qfun(L.r_q).^3.*LX.qpfun(L.r_q).*L.r_q ...
    + 2.*LY.deltap_ana.*LX.qpfun(L.r_q).^2.*L.r_q.^2 + LX.qfun(L.r_q).*(LX.qpfun(L.r_q) + LY.deltap_ana.*LX.qpfun(L.r_q).^2 - ...
    2.*LY.deltap_ana.*qpp).*L.r_q.^2 + LX.qfun(L.r_q).^5.*(LX.kinetic_profiles.betap(L.r_q) - ...
    betapp.*L.r_q) + LX.qfun(L.r_q).^2.*L.r_q.*(LY.deltap_ana.*qpp.*L.r_q + ...
    LX.qpfun(L.r_q).*(-3.*LY.deltap_ana + L.r_q)) - 2.*LX.qfun(L.r_q).^4.*(-(betapp.*L.r_q) + ...
    LX.kinetic_profiles.betap(L.r_q).*(1 + 2.*LX.qpfun(L.r_q).*L.r_q))))./LX.qfun(L.r_q).^4);
    
    V2fun = @(r) interp1(L.r_q, V2, r, 'spline');
    [X2_i, X2_ip] = solve_Xi_eq(L.r_q(1:I), 0, 1, LX.qfun, LX.qpfun, 2, false, V2fun);
    if J == numel(L.r_q)
        neumann = false;
    else
        neumann = true;
    end
    [X2_e, X2_ep] = solve_Xi_eq(L.r_q(I:J-1), 1, 0, LX.qfun, LX.qpfun, 2, neumann);
    
    xijump = -2 * (LY.deltap_ana(I) * LX.qpfun(L.r_q(I)) + 2 * LX.kinetic_profiles.betap(L.r_q(I))) / L.r_q(I);
    xi2p_in =  xijump / (X2_ep(1) - X2_ip(end)) * X2_ip(end);
    xi2p_out = xijump / (X2_ep(1) - X2_ip(end)) * X2_ep(1);
    xi2_in = xijump / (X2_ep(1) - X2_ip(end)) * X2_i(end);
    xi2_out = xijump / (X2_ep(1) - X2_ip(end)) * X2_e(1);
    
    b = rs * xi2p_in / xi2_in;
    c = rs * xi2p_out / xi2_out;
    
    %rs^2/qrs^4 *
    dWt = rs^2/qrs^4 * ( (32 * (b-c) * sigma + 9 * (b-1)*(1-c))/(64*(b-c)) -...
        (betap + sigma) * (3*(b-1)*(c+3))/ (8*(b-c)) - ...
        (betap + sigma)^2 *(c+3)*(b+3)/(4*(b-c)));
    
    % factors eps_val^2 missing
    growth_rate = -pi /sqrt(3) /srs * dWt;
    gr_vec(ii) = growth_rate;
end

%%
figure;
for ii=1:6
    plot(values_eps, gr_vec(ii) * values_eps.^2, 'o-', ...
    'DisplayName', sprintf('$\\Delta q_0 = %.2f, \\hat\\gamma/\\epsilon^2 = %.2f$', deltaq0(ii),gr_vec(ii))); 
    hold on;
end
xlabel('$\epsilon$', 'Interpreter', 'latex', 'Fontsize',14)
ylabel('$\hat \gamma^\mathrm{Bussac}$', 'Interpreter', 'latex', 'Fontsize',14)
legend('Interpreter', 'latex', 'Fontsize',14, 'Box','off')
grid on;
title('$\beta_p(\hat r_s)=1.2$', 'Interpreter', 'latex', 'Fontsize',14)




%% Bussac growth rate
% recompute from scratch.

beta0=[1.2, 1.4, 1.5, 1.8, 2.0];
deltaq0 = [0.03, 0.06, 0.1, 0.15, 0.2, 0.25];

values_eps = [0.08, 0.135, 0.19, 0.245, 0.3];
S2bc=[0.55, 0.49, 0.43, 0.37, 0.31]; %dq0=0.1
%Sbc = [S2bc(1) 0 0];

gr_vec = zeros(6,1);
gr_vec_b = zeros(6,1);
% dq0 = 0.1, beta0 = 1.2, 
for ii =1:6
%ii=2;
    dq0 = deltaq0(ii);
    Sbc = zeros(3,1);
    eps_val = values_eps(5);
    [L, LX] = equilSol('debug',4, 'Nb', 1, ...
        'beta', beta0(1), 'om_pts', 300, 'do_SFL',true, 'residuals_fun', @residuals_iso_static, 'jacobian_fun', @jacobian_iso_static, 'NLtol', 1e-7);
    
    qfct = @(r) (4*(-1+dq0))./(-4 + 6*r.^2 -4*r.^4 + r.^6);
    LX.qfun = qfct;
    LX.qpfun = @(r) (-8*r.*(6 - 8*r.^2 + 3*r.^4)*(-1 + dq0))./(-4 + 6*r.^2 - 4*r.^4 +r.^6).^2;
    
    LX.eps_val = eps_val;LX.Sbc=Sbc;
    LY =equilY(L, LX);
    
    
    
    % Bussac growth rate computation
    % might want to enforce qrs = 1.
    % in fact, might want to really make this whole thing more accurate with
    % respect to the proper position of the q = 1 surface.
    
    [~, I] = min(sqrt((LX.qfun(L.r_q)-1).^2));
    [~, J] = min(sqrt((LX.qfun(L.r_q)-2).^2));
    
    %qrs = LX.qfun(L.r_q(I));
    qrs= 1;
    rs = L.r_q(I);
    srs = rs * LX.qpfun(L.r_q(I));
    betap = -2 * qrs / rs^4 * trapz(L.r_q(1:I), LX.kinetic_profiles.betap(L.r_q(1:I)) .* L.r_q(1:I).^2);
    
    dWc = -rs^2/qrs^2 * ( betap/qrs^2 + (1/rs)^4 * trapz(L.r_q(1:I), (3/qrs + 1./LX.qfun(L.r_q(1:I)) .* (1./LX.qfun(L.r_q(1:I))-1/qrs))));
    
    sigma = (1/rs)^4 *  trapz(L.r_q(1:I), L.r_q(1:I).^3 .* ( 1./LX.qfun(L.r_q(1:I)).^2 - 1) );
    
    % Compute b and c
    
    betapp = -2*L.P.beta;
    qpp = (8*(-1 + dq0)*(24 + 12*L.r_q.^2 - 156*L.r_q.^4 + 208*L.r_q.^6 - 108*L.r_q.^8 + 21*L.r_q.^10))./(-4 + 6*L.r_q.^2 - 4*L.r_q.^4 + L.r_q.^6).^3;
    
    V2 = (L.r_q < rs) .*((-2.*(-2 + LX.qfun(L.r_q)).*(4.*LX.kinetic_profiles.betap(L.r_q).*LX.qfun(L.r_q).^3.*LX.qpfun(L.r_q).*L.r_q ...
    + 2.*LY.deltap_ana.*LX.qpfun(L.r_q).^2.*L.r_q.^2 + LX.qfun(L.r_q).*(LX.qpfun(L.r_q) + LY.deltap_ana.*LX.qpfun(L.r_q).^2 - ...
    2.*LY.deltap_ana.*qpp).*L.r_q.^2 + LX.qfun(L.r_q).^5.*(LX.kinetic_profiles.betap(L.r_q) - ...
    betapp.*L.r_q) + LX.qfun(L.r_q).^2.*L.r_q.*(LY.deltap_ana.*qpp.*L.r_q + ...
    LX.qpfun(L.r_q).*(-3.*LY.deltap_ana + L.r_q)) - 2.*LX.qfun(L.r_q).^4.*(-(betapp.*L.r_q) + ...
    LX.kinetic_profiles.betap(L.r_q).*(1 + 2.*LX.qpfun(L.r_q).*L.r_q))))./LX.qfun(L.r_q).^4);
    
    V2fun = @(r) interp1(L.r_q, V2, r, 'spline');
    [X2_i, X2_ip] = solve_Xi_eq(L.r_q(1:I), 0, 1, LX.qfun, LX.qpfun, 2, false, V2fun);
    if J == numel(L.r_q)
        neumann = false;
    else
        neumann = true;
    end
    [X2_e, X2_ep] = solve_Xi_eq(L.r_q(I:J-1), 1, 0, LX.qfun, LX.qpfun, 2, neumann);
    
    
    
    b = rs * X2_ip(end) / X2_i(end);
    c = rs * X2_ep(1) / X2_e(1);
    bb = -1/2 + rs * X2_ip(end) / X2_i(end) / 4;
    cc = 1/2 + rs * X2_ep(1) / X2_e(1) / 4;
    
    dWt = rs^2/qrs^4 * ( (32 * (b-c) * sigma + 9 * (b-1)*(1-c))/(64*(b-c)) -...
        (betap + sigma) * (3*(b-1)*(c+3))/ (8*(b-c)) - ...
        (betap + sigma)^2 * (c+3)*(b+3)/(4*(b-c)));

    dWtb = rs^2/qrs^4 * ( 8*sigma*(1+bb-cc) +9*bb*(1-cc)-24*bb*cc*(betap+sigma) - 16 * cc * (1+bb) * (betap+sigma)^2) / (16 * (1+bb-cc));
    %dWt =  rs^2/qrs^4 * ( ) / (16 * (1+b-c))
    
    % factors eps_val^2 missing
    growth_rate = -pi /sqrt(3) /srs * dWt
    growth_rate_bussac = -pi /sqrt(3) /srs * dWtb
    gr_vec(ii) = growth_rate;
    gr_vec_b(ii) = growth_rate_bussac;
end

%%

figure;
for ii=1:6
    plot(values_eps, gr_vec(ii) * values_eps.^2, 'o-', ...
    'DisplayName', sprintf('$\\Delta q_0 = %.2f, \\hat\\gamma/\\epsilon^2 = %.2f$', deltaq0(ii),gr_vec(ii))); 
    hold on;
end
xlabel('$\epsilon$', 'Interpreter', 'latex', 'Fontsize',14)
ylabel('$\hat \gamma^\mathrm{Bussac}$', 'Interpreter', 'latex', 'Fontsize',14)
legend('Interpreter', 'latex', 'Fontsize',14, 'Box','off')
grid on;
title('$\beta_p(\hat r_s)=1.2$', 'Interpreter', 'latex', 'Fontsize',14)