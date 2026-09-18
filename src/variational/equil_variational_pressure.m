function pressure = equil_variational_pressure(L, LX, r, R, B)
%EQUIL_VARIATIONAL_PRESSURE Evaluate normalized Pi_parallel and derivatives.
%   Existing EQUIL equations of state return beta-scale maps.  Their
%   pressure and all derivatives are multiplied by epsilon^2 here so that
%   Pi_parallel=mu0*P_parallel/B0^2 is used everywhere in the action.

    if isfield(LX, 'pressure_map') && ~isempty(LX.pressure_map)
        pressure = LX.pressure_map(r, R, B);
        required = {'Pi','Pr','PR','PB','PBB','PrB','PRB','PrR','PRR'};
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
    try
        pressure = eos(LX.kinetic_profiles,r,R,B,'variational');
    catch exception
        if strcmp(exception.identifier,'MATLAB:TooManyInputs')
            error('equilVariational:PressureDerivativeInterface', ...
                ['The equation of state must accept the request ', ...
                 '''variational'' and return the named pressure ', ...
                 'derivative interface through second order.']);
        end
        rethrow(exception)
    end
    required = {'Pi','Pr','PR','PB','PBB','PrB','PRB','PrR','PRR'};
    if ~isstruct(pressure)
        error('equilVariational:PressureDerivativeInterface', ...
            'The variational equation-of-state output must be a struct.');
    end
    for k = 1:numel(required)
        if ~isfield(pressure,required{k})
            error('equilVariational:PressureDerivativeInterface', ...
                'The variational pressure output is missing field %s.', ...
                required{k});
        end
    end
    if ~isfield(pressure,'Pperp')
        pressure.Pperp = pressure.Pi-B.*pressure.PB;
    end
    names = fieldnames(pressure);
    scale = LX.eps_val^2;
    for k = 1:numel(names)
        pressure.(names{k}) = scale*pressure.(names{k});
    end
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
