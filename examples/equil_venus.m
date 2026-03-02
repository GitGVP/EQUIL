%% Example of running EQUIL and creating the HDF5 file used in VENUS-MHD


% Match inputs
q0 = 0.7795296521664806;
q1 = 2.9639395354123486;
s0 = 2 * (q1/q0 -1);
[L, LX] = equilSol('debug',4, 'Nb', 1, 'q0', q0,'s0' , s0, ...
    'beta',1, 'om_pts', 299);

% steeper q profile
LX.qfun = @(r) L.P.q0*(1+0.5*L.P.s0*r.^4);
LX.qpfun = @(r) L.P.q0*(0.5*4*L.P.s0*r.^3);
LX.q_vec = LX.qfun(L.r_q);LX.qp_vec = LX.qpfun(L.r_q);
LX.eps_val = 0.1;

LX.Sbc=zeros(numel(LX.Sbc),1);
LX.S2bc=0; LX.S3bc=0;

LY =equilY(L, LX);

filename = 'equilibrium.h5';
to_venus(LX, LY, filename)


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
