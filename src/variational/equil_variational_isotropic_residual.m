function [residual,cache] = equil_variational_isotropic_residual( ...
        L,LX,x,Bguess)
%EQUIL_VARIATIONAL_ISOTROPIC_RESIDUAL Reduced static-isotropic residual.
%   For the unmodified isotropic equation of state, Bphi=T/R is inserted
%   into the action before variation.  Thus no local B equation is solved
%   or differentiated.  The constant a0/R vacuum contribution is removed
%   analytically, as in the general assembly, to avoid cancellation at
%   small inverse aspect ratio.

    if ~is_native_isotropic(L,LX)
        [residual,cache] = equil_variational_residual(L,LX,x,Bguess);
        return
    end

    state = equil_variational_state(L,LX,x);
    epsilon = state.epsilon;
    pressure = isotropic_pressure(LX,state.r,size(state.R));
    edge_pressure = isotropic_pressure(LX,state.edge.r,size(state.edge.R));
    local = isotropic_local_fields(state,pressure);
    edge = isotropic_local_fields(state.edge,edge_pressure);

    residual = zeros(L.total_dofs,1);
    w = L.w_r;

    % Flux variation of the reduced action.
    flux_derivative_kernel = ...
        epsilon*state.psir.*state.goo./state.J;
    flux_value_kernel = (state.J/epsilon)./state.psir.*( ...
        pressure.Pr+state.T.*state.Tr./state.R.^2);
    row = block_rows(L,1);
    residual(row) = L.B1{1}'*(w.*mean(flux_derivative_kernel,2)) ...
        -L.B0{1}'*(w.*mean(flux_value_kernel,2));
    edge_flux = epsilon*state.edge.psir.*state.edge.goo./state.edge.J;
    residual(row) = residual(row) ...
        -L.B0_edge{1}'*mean(edge_flux,2);
    residual(row) = residual(row)/epsilon^2;

    % Coefficients in the first geometrical variation.  They are the
    % reduced isotropic action written in the same weak, vacuum-subtracted
    % form as equil_variational_assemble.
    [coeff_J,coeff_R,coeff_goo] = ...
        isotropic_coefficients(state,pressure);
    [edge_coeff_J,~,~] = ...
        isotropic_coefficients(state.edge,edge_pressure);

    delta_value = -coeff_J.*state.JoverR-coeff_R;
    delta_derivative = -coeff_J.*state.R.*state.Zw;
    row = block_rows(L,2);
    residual(row) = L.B0{2}'*(w.*mean(delta_value,2)) ...
        +L.B1{2}'*(w.*mean(delta_derivative,2));
    edge_delta_derivative = ...
        -edge_coeff_J.*state.edge.R.*state.edge.Zw;
    residual(row) = residual(row) ...
        -L.B0_edge{2}'*mean(edge_delta_derivative,2);
    residual(row) = residual(row)/epsilon^4;

    omega = state.omega;
    for is = 1:L.P.Ns
        n = is;
        cn = cos(n*omega);
        sn = sin(n*omega);
        deltaJ_value = state.JoverR.*cn+state.R.*( ...
            -n*state.Rr.*cn+n*state.Zr.*sn);
        deltaJ_derivative = state.R.*(state.Zw.*cn+state.Rw.*sn);
        delta_goo_value = -2*n*(state.Rw.*sn+state.Zw.*cn);
        shape_value = coeff_J.*deltaJ_value+coeff_R.*cn ...
            +coeff_goo.*delta_goo_value;
        shape_derivative = coeff_J.*deltaJ_derivative;
        profile = 3+is;
        row = block_rows(L,profile);
        residual(row) = L.B0{profile}'*(w.*mean(shape_value,2)) ...
            +L.B1{profile}'*(w.*mean(shape_derivative,2));
        residual(row) = residual(row)/epsilon^4;
    end

    for ia = 1:L.P.Na
        n = L.P.A_modes(ia)-1;
        cn = cos(n*omega);
        sn = sin(n*omega);
        deltaJ_value = state.JoverR.*sn+state.R.*( ...
            -n*state.Rr.*sn-n*state.Zr.*cn);
        deltaJ_derivative = state.R.*(state.Zw.*sn-state.Rw.*cn);
        delta_goo_value = 2*n*(state.Rw.*cn-state.Zw.*sn);
        shape_value = coeff_J.*deltaJ_value+coeff_R.*sn ...
            +coeff_goo.*delta_goo_value;
        shape_derivative = coeff_J.*deltaJ_derivative;
        profile = 3+L.P.Ns+ia;
        row = block_rows(L,profile);
        residual(row) = L.B0{profile}'*(w.*mean(shape_value,2)) ...
            +L.B1{profile}'*(w.*mean(shape_derivative,2));
        residual(row) = residual(row)/epsilon^4;
    end

    if L.P.vertical_shift
        profile = 4+L.P.Ns+L.P.Na;
        row = block_rows(L,profile);
        deltaJ_derivative = -state.R.*state.Rw;
        residual(row) = L.B1{profile}'*( ...
            w.*mean(coeff_J.*deltaJ_derivative,2))/epsilon^4;
    end

    gauge = mean(state.J./state.R.^2,2)/epsilon^2-state.r/state.a0;
    row = block_rows(L,3);
    residual(row) = L.B0{3}'*(w.*gauge)/epsilon;

    fields = struct;
    fields.pressure = pressure;
    fields.edge = edge;
    fields.action_density = ...
        0.5*state.Bp2-pressure.Pi-0.5*state.T.^2./state.R.^2;
    fields.gauge = gauge;
    fields.flux_derivative_kernel = flux_derivative_kernel;
    fields.flux_value_kernel = flux_value_kernel;
    fields.coeff_J = coeff_J;
    fields.coeff_R = coeff_R;
    fields.coeff_goo = coeff_goo;

    cache.state = state;
    cache.local = local;
    cache.fields = fields;
    cache.residual = residual;
    cache.isotropic_reduced = true;
end

function tf = is_native_isotropic(L,LX)
    tf = strcmp(func2str(L.P.equation_of_state),'isotropic') ...
        && (~isfield(LX,'pressure_map') || isempty(LX.pressure_map));
end

function pressure = isotropic_pressure(LX,r,target_size)
    pressure.Pi = expand_radial(LX.eps_val^2*LX.kinetic_profiles.beta(r), ...
        target_size);
    pressure.Pr = expand_radial(LX.eps_val^2*LX.kinetic_profiles.betap(r), ...
        target_size);
    pressure.PR = zeros(target_size);
    pressure.PB = zeros(target_size);
    pressure.PBB = zeros(target_size);
    pressure.PrB = zeros(target_size);
    pressure.PRB = zeros(target_size);
    pressure.PrR = zeros(target_size);
    pressure.PRR = zeros(target_size);
    pressure.Pperp = pressure.Pi;
end

function value = expand_radial(value,target_size)
    if isscalar(value)
        value = value+zeros(target_size);
    elseif ~isequal(size(value),target_size)
        value = value(:)+zeros(target_size);
    end
end

function local = isotropic_local_fields(state,pressure)
    local.Bphi = state.T./state.R;
    local.B = sqrt(state.Bp2+local.Bphi.^2);
    local.sigma = zeros(size(state.R));
    local.sigma_B = zeros(size(state.R));
    local.G = zeros(size(state.R));
    local.GB = local.B.*state.R.^2./state.T.^2;
    local.pressure = pressure;
    local.iterations = 0;
    local.max_abs_G = 0;
end

function [coeff_J,coeff_R,coeff_goo] = ...
        isotropic_coefficients(state,pressure)
    vacuum_subtracted_T2 = state.T.^2-state.a0^2;
    coeff_J = -0.5*state.Bp2-pressure.Pi ...
        -0.5*vacuum_subtracted_T2./state.R.^2;
    coeff_R = state.J.*vacuum_subtracted_T2./state.R.^3;
    coeff_goo = (state.epsilon*state.psir).^2./(2*state.J);
end

function rows = block_rows(L,profile)
    first = L.profile_starts(profile);
    rows = first:first+L.profile_lengths(profile)-1;
end
