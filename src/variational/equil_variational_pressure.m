function pressure = equil_variational_pressure(L, LX, r, R, B)
%EQUIL_VARIATIONAL_PRESSURE Evaluate normalized Pi_parallel and derivatives.
%   Existing EQUIL equations of state return beta-scale maps.  Their
%   pressure and all derivatives are multiplied by epsilon^2 here so that
%   Pi_parallel=mu0*P_parallel/B0^2 is used everywhere in the action.

    if isfield(LX, 'pressure_map') && ~isempty(LX.pressure_map)
        pressure = LX.pressure_map(r, R, B);
        required = {'Pi','Pr','PR','PB','PBB'};
        for k = 1:numel(required)
            if ~isfield(pressure, required{k})
                error('pressure_map did not return field %s.', required{k});
            end
        end
        pressure = expand_pressure_fields(pressure, size(R));
        if ~isfield(pressure, 'Pperp')
            pressure.Pperp = pressure.Pi-B.*pressure.PB;
        end
        return
    end

    eos = L.P.equation_of_state;
    [~, beta_r, beta_B, beta_R, beta_BB, beta_rB, beta_RB, ...
     beta_rR, beta_RR, ~, ~, ~, ~, ~, ~, beta_parallel, beta_perp] = ...
        eos(LX.kinetic_profiles, r, R, B);

    scale = LX.eps_val^2;
    pressure.Pi = scale*beta_parallel;
    pressure.Pr = scale*beta_r;
    pressure.PR = scale*beta_R;
    pressure.PB = scale*beta_B;
    pressure.PBB = scale*beta_BB;
    pressure.PrB = scale*beta_rB;
    pressure.PRB = scale*beta_RB;
    pressure.PrR = scale*beta_rR;
    pressure.PRR = scale*beta_RR;
    pressure.Pperp = scale*beta_perp;
    pressure = expand_pressure_fields(pressure, size(R));
    pressure.Pperp = pressure.Pi-B.*pressure.PB;
end

function pressure = expand_pressure_fields(pressure, target_size)
    names = fieldnames(pressure);
    for k = 1:numel(names)
        value = pressure.(names{k});
        if isscalar(value)
            pressure.(names{k}) = value + zeros(target_size);
        elseif isequal(size(value), target_size)
            continue
        elseif numel(target_size) == 2 && size(value,1) == target_size(1) ...
                && size(value,2) == 1
            pressure.(names{k}) = value + zeros(target_size);
        else
            error('Pressure field %s has size [%s], expected [%s].', ...
                  names{k}, num2str(size(value)), num2str(target_size));
        end
    end
end
