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
    R0 = 10.334454903971825; B0 =  0.9617498446849172;
    P0 = LX.eps_val^2*B0^2/4./pi/1.0E-07;
    Lref =1; Bref =0.9617498446849172;

    % added constants 
    T0 = 0.01999995793231724; Pax = 0.03999991586463448;
    g0 = -0.6; % TODO ? related to flux at the edge
    
    
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
    %f = LX.eps_val*LY.r_plt.*(1-LX.eps_val^2*LY.t2)./LX.qfun(LY.r_plt);
    %g = LX.eps_val*(1-LX.eps_val^2*LY.t2)./LX.qfun(LY.r_plt);
    g=F./LX.qfun(LY.r_plt)/R0;

    h5create(filename, '/profiles/F',Nr);
    h5write(filename, '/profiles/F', F);
    
    h5create(filename, '/profiles/g', Nr);
    h5write(filename, '/profiles/g', g);
    %h5write(filename, '/profiles/g', g0*ones(Nr,1));
    
    h5create(filename, '/profiles/h', Nr);
    h5write(filename, '/profiles/h', LY.r_plt);
    
    h5create(filename, '/profiles/q', Nr);
    h5write(filename, '/profiles/q', LX.qfun(LY.r_plt));

    h5create(filename, '/profiles/P', Nr);
    h5write(filename, '/profiles/P', LX.kinetic_profiles.beta(LY.r_plt)/LX.eps_val^2);
    
    h5create(filename, '/profiles/Prot', [Nr, Ntheta]);
    h5write(filename, '/profiles/Prot', LX.kinetic_profiles.beta(LY.r_plt).*ones(Nr, Ntheta)/LX.eps_val^2);
    
    h5create(filename, '/profiles/rho', Nr);
    h5write(filename, '/profiles/rho', 1*ones(Nr,1)); % careful
    
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

    %% POST PROCESSING
    beta_poloidal = LX.qfun(LY.r_fine).^2 ./ LY.r_fine .^ 4 .* cumtrapz(LY.r_fine, LY.r_fine.^2 .* (-2 * LX.kinetic_profiles.betap(LY.r_fine)));
    [~, I] = min(sqrt((LX.qfun(LY.r_fine)-1).^2));
    [~, J] = min(sqrt((LX.qfun(LY.r_fine)-2).^2));
    beta_rs = beta_poloidal(I);

    li = 2*LX.qfun(LY.r_fine).^2 ./ LY.r_fine .^ 4 .* cumtrapz(LY.r_fine, LY.r_fine.^3 ./ LX.qfun(LY.r_fine).^2);

    [X2_i, X2_ip] = solve_Xi_eq(LY.r_fine(1:I), 0, 1, LX.qfun, LX.qpfun);
    [X2_e, X2_ep] = solve_Xi_eq(LY.r_fine(I:J), 1, 0, LX.qfun, LX.qpfun);

    deltap = interp1(LY.r_plt, LY.deltap, LY.r_fine, 'spline');
    deltapp = interp1(LY.r_plt, LY.deltapp, LY.r_fine, 'spline');

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
    
   
    %% Quantities not used as input in VENUS_MHD

    %h5create(filename, '/profiles/dqds', Nr);
    %h5write(filename, '/profiles/dqds', LX.qpfun(LY.r_plt));
    %h5create(filename, '/profiles/dPds', Nr);
    %h5write(filename, '/profiles/dPds', Pax/LX.kinetic_profiles.beta(0)*LX.kinetic_profiles.betap(LY.r_plt));
    %h5create(filename, '/normalisation/a', 1);
    %h5write(filename, '/normalisation/a', R0 * LX.eps_val);
end