function to_venus(LX, LY, filename)
%
% Prepares inputs for the stability problem to be solved by VENUS-MHD
% Inputs:
%   LX: EQUIL input structure
%   LY: EQUIL output structure
%   filename: string of the HDF5 file to be created
% Things TODO/check:
% 1. Check that derivatives like dqds, dPds that are output by venus-mhd
% are the same as the ones from EQUIL
% 2. Check conversion factors for pressure
% 3. Check meaning of g and compute it from EQUIL
% 4. Check meaning of T, rho etc. which are necessary for VENUS-MHD
% 5. Check the angular variable (i.e. straight field line coordinates conversion)
% 6. Check sizes (currently L.P.om_pts =299 to have Ntheta=300)
    
    if exist(filename, 'file')
        delete(filename)
    end
    
    Nr = numel(LY.r_plt); Ntheta = numel(LY.omega_plt)-1;
    R0 = 3.1923984520629225; B0 =  0.9367286560046748;
    P0 = LX.eps_val^2*B0^2/4./pi/1.0E-07;
    Lref =1; Bref =0.9617498446849172;
    P0_vmec = 0.425/31830.921665582493;
    eps_vmec=0.32;

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
    
    % epsilon^2 * B0^2 / mu0
    h5create(filename, '/normalisation/P0', 1);
    h5write(filename, '/normalisation/P0', P0);
    
    h5create(filename, '/normalisation/M02', 1);
    h5write(filename, '/normalisation/M02', 0);
    
    h5create(filename, '/normalisation/Lref', 1);
    h5write(filename, '/normalisation/Lref', Lref);
    
    %% PROFILES
    F = -R0*B0*(1 + LX.eps_val.^2.*LY.t2);
    g = F./LX.qfun(LY.r_plt)/R0;
    %g =  -2* LY.psi(end) * ones(Nr,1); 
    h5create(filename, '/profiles/F',Nr);
    h5write(filename, '/profiles/F', F);
    
    h5create(filename, '/profiles/g', Nr);
    h5write(filename, '/profiles/g', g);
    
    % try to match VENUS radial coordinate sqrt(psi)
    h5create(filename, '/profiles/h', Nr);
    h5write(filename, '/profiles/h', LY.r_plt); 
    %h5write(filename, '/profiles/h', sqrt([0;LY.psiN(2)/2;LY.psiN(2:end-1);(1+LY.psiN(end-1))/2;1])); %LY.r_plt
    
    h5create(filename, '/profiles/q', Nr);
    h5write(filename, '/profiles/q', LX.qfun(LY.r_plt));

    h5create(filename, '/profiles/P', Nr);
    h5write(filename, '/profiles/P', LX.kinetic_profiles.beta(LY.r_plt)/LX.eps_val^2/P0_vmec*eps_vmec^2)%/LX.eps_val^2);
    
    h5create(filename, '/profiles/Prot', [Nr, Ntheta]);
    h5write(filename, '/profiles/Prot', LX.kinetic_profiles.beta(LY.r_plt).*ones(Nr, Ntheta)/LX.eps_val^2);
    
    h5create(filename, '/profiles/rho', Nr);
    %h5write(filename, '/profiles/rho', 1-([0;LY.psiN(2)/2;LY.psiN(2:end-1);(1+LY.psiN(end-1))/2;1])); % careful
    h5write(filename, '/profiles/rho', 1-LY.psiN);

    h5create(filename, '/profiles/rhorot', [Nr, Ntheta]);
    h5write(filename, '/profiles/rhorot', 1*ones(Nr, Ntheta)); % careful
    
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
    %h5write(filename, '/profiles/s', sqrt([0;LY.psiN(2)/2;LY.psiN(2:end-1);(1+LY.psiN(end-1))/2;1])); %LY.r_plt

    %% Calculations for post processing
    r=LY.r_fine;
    beta_poloidal = LX.qfun(r).^2 ./ r .^ 4 .* cumtrapz(r, r.^2 .* (-2 * LX.kinetic_profiles.betap(r)));
    [~, I] = min(sqrt((LX.qfun(r)-1).^2));
    [~, J] = min(sqrt((LX.qfun(r)-2).^2));
    beta_rs = beta_poloidal(I);
    
    li = 2*LX.qfun(r).^2 ./ r .^ 4 .* cumtrapz(r, r.^3 ./ LX.qfun(r).^2);
    
    [X2_i, X2_ip] = solve_Xi_eq(r(1:I), 0, 1, LX.qfun, LX.qpfun);
    [X2_e, X2_ep] = solve_Xi_eq(r(I:J), 1, 0, LX.qfun, LX.qpfun);
    
    deltap = interp1(LY.r_plt, LY.deltap, r, 'spline');
    deltapp = interp1(LY.r_plt, LY.deltapp, r, 'spline');
    %deltap=LY.deltap_fine;
    %deltapp = 1 - 2 * LX.qfun(LY.r_fine).^2 .* LX.kinetic_profiles.betap(LY.r_fine) ./ LY.r_fine - LY.deltap_fine .* (3 ./LY.r_fine - 2 * LX.qpfun(LY.r_fine)./LX.qfun(LY.r_fine));

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

    intrU       = LX.eps_val^2 * R0^2 * (W + intUTC);
    growth_rate = pi / (rs^2 * s_shear * sqrt(3)) * intrU;


    %% POST PROCESSING
    
    
    h5create(filename, '/postprocessing/eps', 1);
    h5write(filename, '/postprocessing/eps', LX.eps_val);
    
    h5create(filename, '/postprocessing/betap_rs', 1);
    h5write(filename, '/postprocessing/betap_rs', beta_rs);
    
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
    
    h5create(filename, '/postprocessing/r_fine', numel(LY.r_fine));
    h5write(filename, '/postprocessing/r_fine', LY.r_fine);

    h5create(filename, '/postprocessing/gamma_ana', numel(growth_rate));
    h5write(filename, '/postprocessing/gamma_ana', growth_rate);
    
   
    %% Quantities not used as input in VENUS_MHD

    %h5create(filename, '/profiles/dqds', Nr);
    %h5write(filename, '/profiles/dqds', LX.qpfun(LY.r_plt));
    %h5create(filename, '/profiles/dPds', Nr);
    %h5write(filename, '/profiles/dPds', Pax/LX.kinetic_profiles.beta(0)*LX.kinetic_profiles.betap(LY.r_plt));
    %h5create(filename, '/normalisation/a', 1);
    %h5write(filename, '/normalisation/a', R0 * LX.eps_val);
end