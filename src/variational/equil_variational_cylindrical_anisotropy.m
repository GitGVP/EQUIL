function reference = equil_variational_cylindrical_anisotropy(L,LX,r)
%EQUIL_VARIATIONAL_CYLINDRICAL_ANISOTROPY Prescribed sigma(r,B=1).

    native_pressure = ...
        ~isfield(LX,'pressure_map') || isempty(LX.pressure_map);
    eos_name = func2str(L.P.equation_of_state);
    if native_pressure && any(strcmp(eos_name, ...
            {'isotropic','isotropic_rotating'}))
        reference.sigma = zeros(size(r));
        reference.one_minus_sigma = ones(size(r));
        reference.sigmar = zeros(size(r));
        return
    end

    unit_field = ones(size(r));
    pressure = equil_variational_pressure( ...
        L,LX,r,unit_field,unit_field);
    reference.sigma = pressure.PB;
    reference.one_minus_sigma = 1-reference.sigma;
    if any(~isfinite(reference.one_minus_sigma(:))) ...
            || any(reference.one_minus_sigma(:) ...
                <= L.P.min_one_minus_sigma)
        error(['The cylindrical-reference anisotropy reaches ', ...
            '1-sigma_cyl <= prescribed minimum.']);
    end
    if isfield(pressure,'PrB')
        reference.sigmar = pressure.PrB;
    else
        reference.sigmar = [];
    end
end
