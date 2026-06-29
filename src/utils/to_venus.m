function to_venus(LX, LY, filename, X2_i, X2_ip, X2_e, X2_ep, xi2in, xi2out, growth_rate2, xikp1in, xikp1out)
%
 
    if exist(filename, 'file')
        delete(filename)
    end
    
    Nr = numel(LY.r_plt); Ntheta = numel(LY.omega_plt)-1;
    %R0 = 3.1923984520629225; B0 =  6.9617498446849172;
    %R0 = 1; B0 = 1;
    R0 = 3.1923984520629225; B0 =  0.9617498446849172;
    P0 = LX.eps_val^2*B0^2/4./pi/1.0E-07;
    Lref = 1/LX.eps_val^2; Bref =1;
    %Lref = 100; Bref =50;
    %Lref = 1; Bref =1;
    rho0=1;

    % added constants 
    T0 = 0.01999995793231724;
    
    
    %% GEOMETRY

    
    h5create(filename, '/geometry/R', [Nr, Ntheta]);
    h5write(filename, '/geometry/R', R0*LY.RR_sfl(:,1:end-1));
    
    h5create(filename, '/geometry/Z',[Nr, Ntheta]);
    h5write(filename, '/geometry/Z', R0*LY.ZZ_sfl(:,1:end-1));
    
    %% NORMALISATION
    h5create(filename, '/normalisation/R0', 1);
    h5write(filename, '/normalisation/R0', R0);
    
    h5create(filename, '/normalisation/B0', 1);
    h5write(filename, '/normalisation/B0', B0);
    
    h5create(filename, '/normalisation/Bref', 1);
    h5write(filename, '/normalisation/Bref', Bref);
    
    h5create(filename, '/normalisation/P0', 1);
    h5write(filename, '/normalisation/P0', P0);
    
    h5create(filename, '/normalisation/M02', 1);
    h5write(filename, '/normalisation/M02', 0);
    
    h5create(filename, '/normalisation/Lref', 1);
    h5write(filename, '/normalisation/Lref', Lref);
    
    %% PROFILES
    F = -R0*B0*(1 + LX.eps_val.^2.*LY.t2);
    g = R0*F./LX.qfun(LY.r_plt)*LX.eps_val^2;

%     q = [LX.qfun(LY.r_plt(2)); LX.qfun(LY.r_plt(2:end))];
%     g = R0*F./q*LX.eps_val^2;

    %g =  -2* LY.psi(end) * ones(Nr,1); 
    h5create(filename, '/profiles/F',Nr);
    h5write(filename, '/profiles/F', F);
    
    h5create(filename, '/profiles/g', Nr);
    h5write(filename, '/profiles/g', g);
    
    % try to match VENUS radial coordinate sqrt(psi)
    h5create(filename, '/profiles/h', Nr);
    h5write(filename, '/profiles/h', LY.r_plt);
    
    h5create(filename, '/profiles/q', Nr);
%     h5write(filename, '/profiles/q', q);
    h5write(filename, '/profiles/q', LX.qfun(LY.r_plt));

    h5create(filename, '/profiles/P', Nr);
    h5write(filename, '/profiles/P', LX.kinetic_profiles.beta(LY.r_plt)*LX.eps_val^2/4./pi/1.0E-07*B0^2)
    
    h5create(filename, '/profiles/Pperp', [Nr, Ntheta]);
    h5write(filename, '/profiles/Pperp', LY.betaperp_sfl(:,1:end-1)*LX.eps_val^2/4./pi/1.0E-07*B0^2)

    h5create(filename, '/profiles/Ppar', [Nr, Ntheta]);
    h5write(filename, '/profiles/Ppar', LY.betapar_sfl(:,1:end-1)*LX.eps_val^2/4./pi/1.0E-07*B0^2)

    h5create(filename, '/profiles/Prot', [Nr, Ntheta]);
    h5write(filename, '/profiles/Prot', LX.kinetic_profiles.beta(LY.r_plt).*ones(Nr, Ntheta)*LX.eps_val^2);
    
    h5create(filename, '/profiles/rho', Nr);
    %h5write(filename, '/profiles/rho', 1-LY.psiN);
    h5write(filename, '/profiles/rho', rho0*ones(Nr, 1));

    h5create(filename, '/profiles/rhorot', [Nr, Ntheta]);
    h5write(filename, '/profiles/rhorot', rho0*ones(Nr, Ntheta)); % careful
    
    h5create(filename, '/profiles/T', Nr);
    h5write(filename, '/profiles/T', T0*ones(Nr,1));
    
    h5create(filename, '/profiles/U', Nr);
    h5write(filename, '/profiles/U', zeros(Nr,1));
    
    h5create(filename, '/profiles/Omega', Nr);
    h5write(filename, '/profiles/Omega', zeros(Nr,1));
    
    h5create(filename, '/profiles/Uthi', Nr);
    h5write(filename, '/profiles/Uthi', zeros(Nr,1));
    
    h5create(filename, '/profiles/s', Nr);
    h5write(filename, '/profiles/s', LY.r_plt); 
    %h5write(filename, '/profiles/s', sqrt(LY.psiN)); %LY.r_plt

    %% Calculations for post processing
    r=LY.r_fine;
    dbetapardr = mean(LY.dbetapardr,2).';
    d2betapardrdB = mean(LY.d2betapardrdB,2).';

    beta_poloidal = LX.qfun(r).^2 ./ r .^ 4 .* cumtrapz(r_fine, r_fine.^2 .* (-2 * dbetapardr + d2betapardrdB));
    [~, I] = min(sqrt((LX.qfun(r)-1).^2));
    beta_rs = beta_poloidal(I);
    
    li = 2*LX.qfun(r).^2 ./ r .^ 4 .* cumtrapz(r, r.^3 ./ LX.qfun(r).^2);
    
%     deltap = interp1(LY.r_plt, LY.deltap, r, 'spline');
%     deltapp = interp1(LY.r_plt, LY.deltapp, r, 'spline');
    deltap=LY.deltap_fine;
    deltapp = 1 - 2 * LX.qfun(LY.r_fine).^2 .* LX.kinetic_profiles.betap(LY.r_fine) ./ LY.r_fine - LY.deltap_fine .* (3 ./LY.r_fine - 2 * LX.qpfun(LY.r_fine)./LX.qfun(LY.r_fine));

    rs = LY.r_fine(I);
    s_shear=rs*LX.qpfun(rs);

    zeta  = LX.eps_val * (deltap + r/2);
    zetap = LX.eps_val * (deltapp + 1/2);
    zetas = LX.eps_val * (deltap(I) + rs/2);

    b = rs/4 * (X2_ip(end)/X2_i(end) - 1/rs);
    c = rs/4 * (X2_ep(1) /X2_e(1)  + 3/rs);

    s_hat = li(I)/2 - 1/4;
    l2  = r .* (1./LX.qfun(r) - 1/2);
    nu  = LX.eps_val*rs * ((-3/4 + c*(3/4 + s_hat + beta_rs)) / (1 + b - c));

    intUTC_1 = rs^2/4 * (3*(zetas + nu)*(zetas - LX.eps_val*rs) + nu*zetas*(1 + 4*b));

    to_intUTC_2 = r .* zetap.^2 .* l2.^2 ...
            + 3* l2.^2 .* zeta.^2 ./ r ...
            + r.^3 .* (1/2*(zetap + 3./r.*zeta - LX.eps_val) - LX.eps_val./LX.qfun(r)).^2;

    mask = r <= rs;
    intUTC_2 = trapz(r(mask), to_intUTC_2(mask));
    intUTC   = intUTC_1 + intUTC_2;

    W_1 = -3/4*rs^2 * (zetas - LX.eps_val*rs/2)*(zetas - 3*LX.eps_val*rs/2) ...
      - s_hat * LX.eps_val^2 * rs^4 / 2;

    to_W_2 = r .* zetap.^2 .* l2.^2 ...
       + 3* l2.^2 .* zeta.^2 ./ r ...
       + r.^3 .* (1/2*(zetap + 3*zeta./r) - LX.eps_val*(1./LX.qfun(r) + 1/2)).^2;

    W_2 = trapz(r(mask), to_W_2(mask));
    W   = W_1 - W_2;

    intrU       = (W + intUTC);
    growth_rate = pi / (rs^2 * s_shear * sqrt(3)) * intrU;
%     fprintf('bp(rs) = %f, bp^2-13/144 = %f, (gamma/omega_A) = %.1e, (gamma/omega_A)^2 = %.1e, (gamma/omega_A)new = %.1e \n', beta_rs, beta_rs^2-13/144, growth_rate, growth_rate^2, growth_rate2*LX.eps_val^2)
    fprintf('bp(rs) = %f, rs = %f, (gamma/omega_A) = %.1e, (gamma/omega_A)^2 = %.1e, (gamma/omega_A)new = %.1e \n', beta_rs, r(I), growth_rate, growth_rate^2, growth_rate2*LX.eps_val^2)

    %% POST PROCESSING
    
    
    h5create(filename, '/postprocessing/eps', 1);
    h5write(filename, '/postprocessing/eps', LX.eps_val);

    h5create(filename, '/postprocessing/S2bc', 1);
    h5write(filename, '/postprocessing/S2bc', LX.Sbc(1));

    h5create(filename, '/postprocessing/S3bc', 1);
    h5write(filename, '/postprocessing/S3bc', LX.Sbc(2));

    h5create(filename, '/postprocessing/beta_poloidal', numel(beta_poloidal));
    h5write(filename, '/postprocessing/beta_poloidal', beta_poloidal);
    
    h5create(filename, '/postprocessing/betap_rs', 1);
    h5write(filename, '/postprocessing/betap_rs', beta_rs);

    h5create(filename, '/postprocessing/kappa', 1);
    h5write(filename, '/postprocessing/kappa', LY.kappa(end));

    h5create(filename, '/postprocessing/deltatrig', 1);
    h5write(filename, '/postprocessing/deltatrig', LY.deltatrig(end));

    h5create(filename, '/postprocessing/betap', Nr);
    h5write(filename, '/postprocessing/betap', LX.kinetic_profiles.betap(LY.r_plt));
    
    h5create(filename, '/postprocessing/psiN', numel(LY.psiN));
    h5write(filename, '/postprocessing/psiN', LY.psiN);

    h5create(filename, '/postprocessing/deltap', numel(LY.r_fine));
    h5write(filename, '/postprocessing/deltap', deltap);
    
    h5create(filename, '/postprocessing/deltapp', numel(LY.r_fine));
    h5write(filename, '/postprocessing/deltapp', deltapp);
    
    h5create(filename, '/postprocessing/li', numel(li));
    h5write(filename, '/postprocessing/li', li);
    
    h5create(filename, '/postprocessing/X2_i', numel(X2_i));
    h5write(filename, '/postprocessing/X2_i', X2_i);
    
    h5create(filename, '/postprocessing/X2_ip', numel(X2_i));
    h5write(filename, '/postprocessing/X2_ip', X2_ip);
    
    h5create(filename, '/postprocessing/X2_e', numel(X2_e));
    h5write(filename, '/postprocessing/X2_e', X2_e);
    
    h5create(filename, '/postprocessing/X2_ep', numel(X2_e));
    h5write(filename, '/postprocessing/X2_ep', X2_ep);

    h5create(filename, '/postprocessing/new_xi2in', numel(xi2in));
    h5write(filename, '/postprocessing/new_xi2in', xi2in);

    h5create(filename, '/postprocessing/new_xi2out', numel(xi2out));
    h5write(filename, '/postprocessing/new_xi2out', xi2out);

    h5create(filename, '/postprocessing/xikin', numel(xikp1in));
    h5write(filename, '/postprocessing/xikin', xikp1in);

    h5create(filename, '/postprocessing/xikout', numel(xikp1out));
    h5write(filename, '/postprocessing/xikout', xikp1out);
    
    h5create(filename, '/postprocessing/r_fine', numel(LY.r_fine));
    h5write(filename, '/postprocessing/r_fine', LY.r_fine);

    h5create(filename, '/postprocessing/gamma_ana', 1);
    h5write(filename, '/postprocessing/gamma_ana', growth_rate);

    h5create(filename, '/postprocessing/newgamma_ana', 1);
    h5write(filename, '/postprocessing/newgamma_ana', growth_rate2*LX.eps_val^2);

%     h5create(filename, '/postprocessing/gamma_bussac', 1);
%     h5write(filename, '/postprocessing/gamma_bussac', growthrate_bussac*LX.eps_val^2);

%     h5create(filename, '/postprocessing/L1', numel(LY.L1));
%     h5write(filename, '/postprocessing/L1', real(LY.L1));

%     h5create(filename, '/postprocessing/N0', numel(LY.N0));
%     h5write(filename, '/postprocessing/N0', LY.N0);
% 
%     h5create(filename, '/postprocessing/Nm1', numel(LY.Nm1));
%     h5write(filename, '/postprocessing/Nm1', real(LY.Nm1));
% 
%     h5create(filename, '/postprocessing/Nimag', numel(LY.Nm1));
%     h5write(filename, '/postprocessing/Nimag', imag(LY.Nm1));

%     h5create(filename, '/postprocessing/Mreal', numel(LY.Mm1));
%     h5write(filename, '/postprocessing/Mreal', real(LY.Mm1));
% 
%     h5create(filename, '/postprocessing/Mm1', numel(LY.Mm1));
%     h5write(filename, '/postprocessing/Mm1', imag(LY.Mm1));
    
    h5create(filename, '/postprocessing/JoverBm1', numel(LY.JoverBm1));
    h5write(filename, '/postprocessing/JoverBm1', real(LY.JoverBm1));

%     h5create(filename, '/postprocessing/oneoverB1', numel(LY.oneoverB1));
%     h5write(filename, '/postprocessing/oneoverB1', real(LY.oneoverB1));
% 
%     h5create(filename, '/postprocessing/JoverB1', numel(LY.JoverB1));
%     h5write(filename, '/postprocessing/JoverB1', real(LY.JoverB1));

    h5create(filename, '/postprocessing/Y0_ana', numel(LY.Y0_ana));
    h5write(filename, '/postprocessing/Y0_ana', real(LY.Y0_ana));

    h5create(filename, '/postprocessing/Y0', numel(LY.Y0));
    h5write(filename, '/postprocessing/Y0', real(LY.Y0));

    h5create(filename, '/postprocessing/Y1', numel(-LY.l1p));
    h5write(filename, '/postprocessing/Y1', real(-LY.l1p));

    h5create(filename, '/postprocessing/Z1', numel(LY.l1));
    h5write(filename, '/postprocessing/Z1', imag(LY.l1));
    
%     h5create(filename, '/postprocessing/g22', [Nr, Ntheta]);
%     h5write(filename, '/postprocessing/g22', LY.gtt(:,1:end-1));
% 
%     h5create(filename, '/postprocessing/g12', [Nr, Ntheta]);
%     h5write(filename, '/postprocessing/g12', LY.grt(:,1:end-1));
% 
%     h5create(filename, '/postprocessing/g11', [Nr, Ntheta]);
%     h5write(filename, '/postprocessing/g11', LY.grr(:,1:end-1));

%     h5create(filename, '/postprocessing/Ja', [Nr, Ntheta]);
%     h5write(filename, '/postprocessing/Ja', LY.J_SFL(:,1:end-1)./LY.r_plt*R0^2);
   
    %% Quantities not used as input in VENUS_MHD

    %h5create(filename, '/profiles/dqds', Nr);
    %h5write(filename, '/profiles/dqds', LX.qpfun(LY.r_plt));
    %h5create(filename, '/profiles/dPds', Nr);
    %h5write(filename, '/profiles/dPds', Pax/LX.kinetic_profiles.beta(0)*LX.kinetic_profiles.betap(LY.r_plt));
    %h5create(filename, '/normalisation/a', 1);
    %h5write(filename, '/normalisation/a', R0 * LX.eps_val);
end