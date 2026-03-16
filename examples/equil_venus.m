%% Example of running EQUIL and creating the HDF5 file used in VENUS-MHD


% Match inputs
q0 = 0.7795296521664806;
q1 = 2.9639395354123486;
s0 = 2 * (q1/q0 -1);

eps_val = 0.0975; Sbc = [-1.2 0.4 0];
[L, LX] = equilSol('debug',4, 'Nb', 1, 'q0', q0,'s0' , s0, ...
    'beta',1, 'om_pts', 299, 'do_SFL',true);

LX.eps_val = eps_val;LX.Sbc=Sbc;


LY =equilY(L, LX);

% increase beta
for betas = [3 4.6 4.6 4.6]
    L.P.beta = betas;
    LX = equilX(L);
    LX.eps_val = eps_val;LX.Sbc=Sbc;
    beta_interp = L.P.beta * (1 - LY.psiN);
    ord_b = 0:2:14;
    cb = (bsxfun(@power, L.r_q(:), ord_b) \ beta_interp(:)).';     % 1xK
    beta_fit  = @(rr) reshape(    bsxfun(@power, rr(:), ord_b)   * cb.',          size(rr));
    betap_fit = @(rr) reshape(    bsxfun(@power, rr(:), ord_b-1) * (ord_b.*cb).', size(rr));

    LX.kinetic_profiles.beta  = @(rr) reshape(interp1(L.r_q, beta_fit(L.r_q),  rr(:), 'pchip', 'extrap'), size(rr));
    LX.kinetic_profiles.betap = @(rr) reshape(interp1(L.r_q, betap_fit(L.r_q), rr(:), 'pchip', 'extrap'), size(rr));
    LX.x = LY.x;
    LY = equilY(L,LX);
end
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
