function to_venus_cleaner(LX, LY, filename)
%
    if exist(filename, 'file')
        delete(filename)
    end
    Nr = numel(LY.r_plt); Ntheta = numel(LY.omega_plt)-1;
    R0 = 3.1923984520629225; B0 =  0.9617498446849172;
    P0 = LX.eps_val^2*B0^2/4./pi/1.0E-07;
    Lref = 1/LX.eps_val^2; Bref =1;
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

    h5create(filename, '/profiles/F',Nr);
    h5write(filename, '/profiles/F', F);
    
    h5create(filename, '/profiles/g', Nr);
    h5write(filename, '/profiles/g', g);
    
    h5create(filename, '/profiles/h', Nr);
    h5write(filename, '/profiles/h', LY.r_plt);
    
    h5create(filename, '/profiles/q', Nr);
    h5write(filename, '/profiles/q', LX.qfun(LY.r_plt));

    h5create(filename, '/profiles/P', Nr);
    h5write(filename, '/profiles/P', LX.kinetic_profiles.beta(LY.r_plt)*LX.eps_val^2/4./pi/1.0E-07*B0^2)

    h5create(filename, '/profiles/Prot', [Nr, Ntheta]);
    h5write(filename, '/profiles/Prot', LX.kinetic_profiles.beta(LY.r_plt).*ones(Nr, Ntheta)*LX.eps_val^2);
    
    h5create(filename, '/profiles/rho', Nr);
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
end