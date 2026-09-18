function varargout = thermal_and_deBlank(kinetic_profiles,r,RR,BB,varargin)
    thermal_profiles = struct( ...
        'beta',  kinetic_profiles.beta, ...
        'betap', kinetic_profiles.betap);

    anisotropic_profiles = kinetic_profiles;
    anisotropic_profiles.beta  = kinetic_profiles.beta_A;
    anisotropic_profiles.betap = kinetic_profiles.betap_A;

    if nargin > 4
        pressure_th = isotropic(thermal_profiles,r,RR,BB,'variational');
        pressure_A = de_Blank( ...
            anisotropic_profiles,r,RR,BB,'variational');
        names = fieldnames(pressure_th);
        for k = 1:numel(names)
            pressure_th.(names{k}) = pressure_th.(names{k}) ...
                +pressure_A.(names{k});
        end
        varargout{1} = pressure_th;
        return
    end

    nout = 17;
    out_th = cell(1,nout);
    out_A  = cell(1,nout);

    [out_th{:}] = isotropic( ...
        thermal_profiles,r,RR,BB);

    [out_A{:}] = de_Blank( ...
        anisotropic_profiles,r,RR,BB);

    out = cell(1,nout);
    for k = 1:nout
        out{k} = out_th{k} + out_A{k};
    end

    varargout = out(1:nargout);
end
