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
    
    Nr = numel(LY.r_plt); Ntheta = numel(LY.omega_plt);
    R0 = 10.334454903971825;B0 =  0.9617498446849172;
    P0 = 31830.921665582508;
    Lref =1; Bref =1;

    % added constants 
    T0 = 0.01999995793231724; Pax = 0.03999991586463448;
    g0 = -0.6; % TODO ? related to flux at the edge
    
    %% GEOMETRY
    h5create(filename, '/geometry/R', [Nr, Ntheta]);
    h5write(filename, '/geometry/R', R0*LY.RR);
    
    h5create(filename, '/geometry/Z',[Nr, Ntheta]);
    h5write(filename, '/geometry/Z', R0*LY.ZZ);
    
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
    h5create(filename, '/profiles/F',Nr);
    h5write(filename, '/profiles/F', F);
    
    h5create(filename, '/profiles/g', Nr);
    h5write(filename, '/profiles/g', g0*ones(Nr,1));
    
    h5create(filename, '/profiles/h', Nr);
    h5write(filename, '/profiles/h', LY.r_plt);
    
    h5create(filename, '/profiles/q', Nr);
    h5write(filename, '/profiles/q', LX.qfun(LY.r_plt));

    h5create(filename, '/profiles/P', Nr);
    h5write(filename, '/profiles/P', Pax/LX.kinetic_profiles.beta(0)*LX.kinetic_profiles.beta(LY.r_plt));
    
    h5create(filename, '/profiles/Prot', [Nr, Ntheta]);
    h5write(filename, '/profiles/Prot', Pax/LX.kinetic_profiles.beta(0)*LX.kinetic_profiles.beta(LY.r_plt).*ones(size(LY.RR)));
    
    h5create(filename, '/profiles/rho', Nr);
    h5write(filename, '/profiles/rho', 1-LY.r_plt.^2); % careful
    
    h5create(filename, '/profiles/rhorot', [Nr, Ntheta]);
    h5write(filename, '/profiles/rhorot', 1-LY.RR.^2); % careful

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
    %% Quantities not used as input in VENUS_MHD

    %h5create(filename, '/profiles/dqds', Nr);
    %h5write(filename, '/profiles/dqds', LX.qpfun(LY.r_plt));
    %h5create(filename, '/profiles/dPds', Nr);
    %h5write(filename, '/profiles/dPds', Pax/LX.kinetic_profiles.beta(0)*LX.kinetic_profiles.betap(LY.r_plt));
    %h5create(filename, '/normalisation/a', 1);
    %h5write(filename, '/normalisation/a', R0 * LX.eps_val);
end