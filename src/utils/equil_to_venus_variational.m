function info = equil_to_venus_variational(L,LX,LY,filename,varargin)
%EQUIL_TO_VENUS_VARIATIONAL Write the dimensionless VENUS-MHD interface.
%   The stability interface contains native SFL geometry, tensors, and the
%   small set of their derivatives used by MishkaA/B and
%   MishkaAnisotropicA/B. VENUS interpolates these fields but does not
%   differentiate sampled geometry.
%
%   Optional name/value pair: 'gamma', adiabatic index (default 5/3).
%
%   HDF5 implementation:
%   The output file is opened exactly once. All groups, datasets and root
%   attributes are created through the same low-level HDF5 file handle.

    parser = inputParser;
    addParameter(parser,'gamma',5/3,@(x) isscalar(x) && isfinite(x) && x > 0);
    parse(parser,varargin{:});

    required_sfl = {'RR_sfl','ZZ_sfl','Rs_sfl','Rtheta_sfl', ...
        'Zs_sfl','Ztheta_sfl','g11_sfl','g12_sfl','g22_sfl', ...
        'Ja_sfl','B2_sfl','dR2ds_sfl','dR2dtheta_sfl', ...
        'dJads_sfl','dJadtheta_sfl','dg12dtheta_sfl', ...
        'dg22dtheta_sfl', ...
        'dB2ds_sfl','dB2dtheta_sfl','Pi_parallel_sfl', ...
        'Pi_perp_sfl','Dparallel_sfl','Dperp_sfl', ...
        'kappas_sfl','kappatheta_sfl','kappaphi_sfl', ...
        'current_s_sfl','current_theta_sfl','current_phi_sfl'};

    if ~L.P.do_SFL || ~all(isfield(LY,required_sfl))
        error(['VENUS export requires straight-field-line coordinates. ', ...
            'Run equilVariationalSol with ''do_SFL'',true.']);
    end

    %% Profiles and closure quantities

    s = LY.r_plt(:);

    q = LX.qfun(s);
    q = q(:);

    dqds = LX.qpfun(s);
    dqds = dqds(:);

    epsilon = LX.eps_val;

    T = LY.a0 + epsilon^2*LY.t2(:);
    dTds = epsilon^2*LY.t2p(:);

    one_minus_sigma_cyl = LY.one_minus_sigma_cyl(:);
    sigma_cyl_r = LY.sigma_cyl_r(:);

    % In anisotropy, R^2 B^phi=T/(1-sigma) is not a flux function.
    % The one-dimensional F below is the auxiliary profile which gives the
    % native SFL contravariant field through g=epsilon^2*F/q and
    % B^theta=g/Ja. It reduces to the usual toroidal F in isotropy.
    F = -T./one_minus_sigma_cyl;

    dFds = -dTds./one_minus_sigma_cyl ...
        - T.*sigma_cyl_r./one_minus_sigma_cyl.^2;

    g = epsilon^2*F./q;
    dgds = epsilon^2*(dFds.*q - F.*dqds)./q.^2;

    P = epsilon^2*LX.kinetic_profiles.beta(s);
    P = P(:);

    dPds = epsilon^2*LX.kinetic_profiles.betap(s);
    dPds = dPds(:);

    ntheta = size(LY.RR_sfl,2)-1;
    theta = (0:ntheta-1)'*(2*pi/ntheta);

    Dparallel_total = LY.Dparallel_sfl(:,1:ntheta);
    Dperp_total     = LY.Dperp_sfl(:,1:ntheta);

    is_anisotropic = max(abs(LY.Pi_parallel(:)-LY.Pi_perp(:))) > ...
        1e-10*max(1,max(abs(LY.Pi_parallel(:))));

    pressure_model = func2str(L.P.equation_of_state);

    has_split_closure = ~is_anisotropic || ...
        strcmp(pressure_model,'thermal_and_deBlank');

    if is_anisotropic && abs(L.P.mach20) > 1e-13
        error(['The anisotropic VENUS model is static and neglects the ', ...
            'explicit R pressure derivative; export requires mach20=0.']);
    end

    if has_split_closure
        Dparallel_A = Dparallel_total - dPds;
        Dperp_A     = Dperp_total     - dPds;

        Pi_parallel_A = LY.Pi_parallel_sfl(:,1:ntheta) - P;
        Pi_perp_A     = LY.Pi_perp_sfl(:,1:ntheta) - P;
    end

    %% Magnetic-field consistency check

    Btheta_venus = -LY.Btheta_sfl(:,1:ntheta);
    Bphi_venus   = -LY.Bphi_sfl(:,1:ntheta);

    Btheta_from_profiles = g./LY.Ja_sfl(:,1:ntheta);
    Bphi_from_profiles = q.*Btheta_from_profiles;

    relative_Btheta_error = max(abs( ...
        Btheta_from_profiles(:)-Btheta_venus(:))) ...
        / max(1,max(abs(Btheta_venus(:))));

    relative_Bphi_error = max(abs( ...
        Bphi_from_profiles(:)-Bphi_venus(:))) ...
        / max(1,max(abs(Bphi_venus(:))));

    assert(relative_Btheta_error < 1e-7 && relative_Bphi_error < 1e-7, ...
        ['The exported profiles do not reconstruct the native SFL ', ...
         'magnetic field (Btheta %.3e, Bphi %.3e).'], ...
        relative_Btheta_error,relative_Bphi_error);


    %% Create HDF5 file ONCE

    % H5F_ACC_TRUNC replaces an existing file atomically at creation level.
    file_id = H5F.create( ...
        filename, ...
        'H5F_ACC_TRUNC', ...
        'H5P_DEFAULT', ...
        'H5P_DEFAULT');

    % Guarantee closure even if any dataset/attribute write throws.
    file_cleanup = onCleanup(@() H5F.close(file_id));


    %% Create groups

    groups = { ...
        '/coordinates', ...
        '/geometry', ...
        '/tensors', ...
        '/magnetic', ...
        '/anisotropy', ...
        '/curvature', ...
        '/current', ...
        '/profiles', ...
        '/postprocessing'};

    for k = 1:numel(groups)
        create_group(file_id,groups{k});
    end


    %% Coordinates

    write_dataset(file_id,'/coordinates/s',s);
    write_dataset(file_id,'/coordinates/theta',theta);


    %% Geometry

    write_dataset(file_id,'/geometry/R', ...
        LY.RR_sfl(:,1:ntheta));

    write_dataset(file_id,'/geometry/Z', ...
        LY.ZZ_sfl(:,1:ntheta));

    write_dataset(file_id,'/geometry/Rs', ...
        LY.Rs_sfl(:,1:ntheta));

    write_dataset(file_id,'/geometry/Rtheta', ...
        LY.Rtheta_sfl(:,1:ntheta));

    write_dataset(file_id,'/geometry/Zs', ...
        LY.Zs_sfl(:,1:ntheta));

    write_dataset(file_id,'/geometry/Ztheta', ...
        LY.Ztheta_sfl(:,1:ntheta));


    %% Tensors

    write_dataset(file_id,'/tensors/g11', ...
        LY.g11_sfl(:,1:ntheta));

    write_dataset(file_id,'/tensors/g12', ...
        LY.g12_sfl(:,1:ntheta));

    write_dataset(file_id,'/tensors/g22', ...
        LY.g22_sfl(:,1:ntheta));

    write_dataset(file_id,'/tensors/Ja', ...
        LY.Ja_sfl(:,1:ntheta));

    write_dataset(file_id,'/tensors/B2', ...
        LY.B2_sfl(:,1:ntheta));

    write_dataset(file_id,'/tensors/dR2ds', ...
        LY.dR2ds_sfl(:,1:ntheta));

    write_dataset(file_id,'/tensors/dR2dtheta', ...
        LY.dR2dtheta_sfl(:,1:ntheta));

    write_dataset(file_id,'/tensors/dJads', ...
        LY.dJads_sfl(:,1:ntheta));

    write_dataset(file_id,'/tensors/dJadtheta', ...
        LY.dJadtheta_sfl(:,1:ntheta));

    write_dataset(file_id,'/tensors/dg12dtheta', ...
        LY.dg12dtheta_sfl(:,1:ntheta));

    write_dataset(file_id,'/tensors/dg22dtheta', ...
        LY.dg22dtheta_sfl(:,1:ntheta));

    write_dataset(file_id,'/tensors/dB2ds', ...
        LY.dB2ds_sfl(:,1:ntheta));

    write_dataset(file_id,'/tensors/dB2dtheta', ...
        LY.dB2dtheta_sfl(:,1:ntheta));


    %% Magnetic field

    write_dataset(file_id,'/magnetic/Bs', ...
        LY.Bs_sfl(:,1:ntheta));

    % VENUS uses the opposite global field orientation to EQUIL, consistent
    % with the negative F and g profiles above. Curvature is unchanged by
    % this simultaneous reversal of Btheta and Bphi.
    write_dataset(file_id,'/magnetic/Btheta',Btheta_venus);
    write_dataset(file_id,'/magnetic/Bphi',Bphi_venus);


    %% Anisotropy

    write_dataset(file_id,'/anisotropy/B', ...
        LY.B_sfl(:,1:ntheta));

    write_dataset(file_id,'/anisotropy/sigma', ...
        LY.sigma_sfl(:,1:ntheta));

    write_dataset(file_id,'/anisotropy/Pi_parallel', ...
        LY.Pi_parallel_sfl(:,1:ntheta));

    write_dataset(file_id,'/anisotropy/Pi_perp', ...
        LY.Pi_perp_sfl(:,1:ntheta));

    if has_split_closure
        write_dataset(file_id,'/anisotropy/Pi_parallel_A', ...
            Pi_parallel_A);

        write_dataset(file_id,'/anisotropy/Pi_perp_A', ...
            Pi_perp_A);

        write_dataset(file_id,'/anisotropy/Dparallel_A', ...
            Dparallel_A);

        write_dataset(file_id,'/anisotropy/Dperp_A', ...
            Dperp_A);
    end

    % Native total-pressure closure coefficients for the unsplit
    % perpendicular anisotropic model.
    write_dataset(file_id,'/anisotropy/Dparallel', ...
        Dparallel_total);

    write_dataset(file_id,'/anisotropy/Dperp', ...
        Dperp_total);


    %% Curvature

    write_dataset(file_id,'/curvature/kappas', ...
        LY.kappas_sfl(:,1:ntheta));

    write_dataset(file_id,'/curvature/kappatheta', ...
        LY.kappatheta_sfl(:,1:ntheta));

    write_dataset(file_id,'/curvature/kappaphi', ...
        LY.kappaphi_sfl(:,1:ntheta));


    %% Current

    % Current reverses with the global B orientation used by VENUS.
    write_dataset(file_id,'/current/js', ...
        -LY.current_s_sfl(:,1:ntheta));

    write_dataset(file_id,'/current/jtheta', ...
        -LY.current_theta_sfl(:,1:ntheta));

    write_dataset(file_id,'/current/jphi', ...
        -LY.current_phi_sfl(:,1:ntheta));


    %% Profiles

    write_dataset(file_id,'/profiles/h',s);
    write_dataset(file_id,'/profiles/q',q);
    write_dataset(file_id,'/profiles/dqds',dqds);

    write_dataset(file_id,'/profiles/F',F);
    write_dataset(file_id,'/profiles/dFds',dFds);

    write_dataset(file_id,'/profiles/g',g);
    write_dataset(file_id,'/profiles/dgds',dgds);

    write_dataset(file_id,'/profiles/P',P);
    write_dataset(file_id,'/profiles/dPds',dPds);

    write_dataset(file_id,'/profiles/rho', ...
        ones(size(s)));

    write_dataset(file_id,'/profiles/gamma', ...
        parser.Results.gamma*ones(size(s)));

    write_dataset(file_id,'/postprocessing/eps',epsilon);
    write_dataset(file_id,'/postprocessing/bp',LY.bp);
    write_dataset(file_id,'/postprocessing/bp_bussac',LY.bp_bussac);
    write_dataset(file_id,'/postprocessing/tau',LY.tau);
    if isfield(LY,'sfl_force_parallel_error_raw')
        write_dataset(file_id,'/postprocessing/sfl_force_parallel_error_raw', ...
            LY.sfl_force_parallel_error_raw);
        write_dataset(file_id, ...
            '/postprocessing/sfl_curvature_parallel_error_raw', ...
            LY.sfl_curvature_parallel_error_raw);
        write_dataset(file_id,'/postprocessing/sfl_force_closure_error', ...
            LY.sfl_force_closure_error);
    end


    %% Root attributes

    write_attribute(file_id,'schema_name', ...
        'equil_variational_venus');

    write_attribute(file_id,'schema_version', ...
        int32(7));

    write_attribute(file_id,'normalisation', ...
        'dimensionless EQUIL');

    write_attribute(file_id,'coordinates', ...
        's=r/a; theta straight-field-line; phi geometric');

    write_attribute(file_id,'magnetic_field', ...
        ['g=(d psi/ds)/s; B^theta=g*s/J; ', ...
         'F=q*g/epsilon^2 (R^2 B^phi only in isotropy)']);

    write_attribute(file_id,'native_sfl_fields', ...
        'geometry and Mishka derivatives evaluated analytically by EQUIL');

    write_attribute(file_id,'pressure', ...
        'P=mu0*p/B0^2');

    write_attribute(file_id,'pressure_model', ...
        pressure_model);

    write_attribute(file_id,'split_pressure_closure_available', ...
        int32(has_split_closure));

    write_attribute(file_id,'anisotropic_derivatives', ...
        ['Dparallel and Dperp are total-pressure partial derivatives; ', ...
         'Dparallel_A and Dperp_A are species-A partial derivatives. ', ...
         'All are with respect to s at fixed R and B']);

    write_attribute(file_id,'equilibrium_current', ...
        ['contravariant SFL components; perpendicular current from the ', ...
         'anisotropic force balance and jtheta fixed by Ampere''s law']);

    write_attribute(file_id,'equil_radial_discretization', ...
        L.P.radial_discretization);

    write_attribute(file_id,'equil_m', ...
        int32(L.P.m));

    write_attribute(file_id,'equil_nq', ...
        int32(L.P.nq));

    write_attribute(file_id,'equil_spline_p', ...
        int32(L.P.spline_p));

    write_attribute(file_id,'equil_om_pts', ...
        int32(L.P.om_pts));

    write_attribute(file_id,'equil_Ns', ...
        int32(L.P.Ns));

    write_attribute(file_id,'equil_NLtol', ...
        double(L.P.NLtol));


    %% Close file

    % Clearing the onCleanup object closes file_id here.
    clear file_cleanup


    %% Return information

    info.filename = filename;
    info.radial_points = numel(s);
    info.poloidal_points = ntheta;
    info.theta_periodicity_error = ...
        max(abs(LY.theta_SFL_periodicity_error));

    fprintf( ...
        'Wrote dimensionless SFL equilibrium %s (%d x %d).\n', ...
        filename,numel(s),ntheta);
end


% ========================================================================
% HDF5 helpers
% ========================================================================

function create_group(file_id,path)

    group_id = H5G.create( ...
        file_id, ...
        path, ...
        'H5P_DEFAULT', ...
        'H5P_DEFAULT', ...
        'H5P_DEFAULT');

    H5G.close(group_id);
end


function write_dataset(file_id,path,value)

    value = double(value);

    % MATLAB stores dimensions in column-major order whereas the HDF5
    % low-level interface describes dimensions in C order. Reversing the
    % dimensions here preserves the shape seen by h5read/h5py.
    if isscalar(value)
        space_id = H5S.create('H5S_SCALAR');
    else
        dims = fliplr(size(value));
        space_id = H5S.create_simple(numel(dims),dims,[]);
    end

    type_id = H5T.copy('H5T_NATIVE_DOUBLE');

    try
        dataset_id = H5D.create( ...
            file_id, ...
            path, ...
            type_id, ...
            space_id, ...
            'H5P_DEFAULT');

        try
            H5D.write( ...
                dataset_id, ...
                'H5ML_DEFAULT', ...
                'H5S_ALL', ...
                'H5S_ALL', ...
                'H5P_DEFAULT', ...
                value);
        catch ME
            H5D.close(dataset_id);
            rethrow(ME);
        end

        H5D.close(dataset_id);

    catch ME
        H5S.close(space_id);
        H5T.close(type_id);
        rethrow(ME);
    end

    H5S.close(space_id);
    H5T.close(type_id);
end


function write_attribute(file_id,name,value)

    if ischar(value) || isstring(value)

        value = char(value);

        space_id = H5S.create('H5S_SCALAR');
        type_id = H5T.copy('H5T_C_S1');

        % Fixed-length string attribute. NULL padding gives the same
        % effective representation expected by h5readatt/h5py.
        H5T.set_size(type_id,numel(value));
        H5T.set_strpad(type_id,'H5T_STR_NULLPAD');

        try
            attr_id = H5A.create( ...
                file_id, ...
                name, ...
                type_id, ...
                space_id, ...
                'H5P_DEFAULT');

            try
                H5A.write(attr_id,type_id,value);
            catch ME
                H5A.close(attr_id);
                rethrow(ME);
            end

            H5A.close(attr_id);

        catch ME
            H5S.close(space_id);
            H5T.close(type_id);
            rethrow(ME);
        end

        H5S.close(space_id);
        H5T.close(type_id);

    elseif isa(value,'int32')

        space_id = H5S.create('H5S_SCALAR');
        type_id = H5T.copy('H5T_NATIVE_INT32');

        try
            attr_id = H5A.create( ...
                file_id, ...
                name, ...
                type_id, ...
                space_id, ...
                'H5P_DEFAULT');

            try
                H5A.write(attr_id,type_id,value);
            catch ME
                H5A.close(attr_id);
                rethrow(ME);
            end

            H5A.close(attr_id);

        catch ME
            H5S.close(space_id);
            H5T.close(type_id);
            rethrow(ME);
        end

        H5S.close(space_id);
        H5T.close(type_id);

    else

        value = double(value);

        space_id = H5S.create('H5S_SCALAR');
        type_id = H5T.copy('H5T_NATIVE_DOUBLE');

        try
            attr_id = H5A.create( ...
                file_id, ...
                name, ...
                type_id, ...
                space_id, ...
                'H5P_DEFAULT');

            try
                H5A.write(attr_id,type_id,value);
            catch ME
                H5A.close(attr_id);
                rethrow(ME);
            end

            H5A.close(attr_id);

        catch ME
            H5S.close(space_id);
            H5T.close(type_id);
            rethrow(ME);
        end

        H5S.close(space_id);
        H5T.close(type_id);
    end
end
