function LX = equil_variational_set_beta(LX, beta_value)
%EQUIL_VARIATIONAL_SET_BETA Rescale the prescribed radial pressure profile.

    if ~isfield(LX,'beta_shape') || ~isfield(LX,'betap_shape')
        error('Automatic beta continuation requires beta_shape profiles.');
    end
    target_beta = LX.beta_value;
    beta_shape = LX.beta_shape;
    betap_shape = LX.betap_shape;
    LX.kinetic_profiles.beta = @(r) beta_value*beta_shape(r);
    LX.kinetic_profiles.betap = @(r) beta_value*betap_shape(r);

    % In the runaway closure gamma(r) is the coefficient of B in the
    % parallel pressure, not the dimensionless anisotropy ratio itself:
    % Pi_parallel = epsilon^2*(beta+B*gamma).  It therefore has to vanish
    % and ramp with beta during the automatic pressure continuation.  The
    % old continuation left gamma at its target value even in the beta=0
    % seed, so the nominal vacuum seed still carried anisotropic pressure.
    if isfield(LX.kinetic_profiles,'gamma')
        if target_beta == 0
            if beta_value ~= 0
                error(['Cannot rescale a runaway gamma profile from a ', ...
                       'zero target beta.']);
            end
            fraction = 0;
        else
            fraction = beta_value/target_beta;
        end
        target_gamma = LX.kinetic_profiles.gamma;
        target_gammap = LX.kinetic_profiles.gammap;
        LX.kinetic_profiles.gamma = @(r) fraction*target_gamma(r);
        LX.kinetic_profiles.gammap = @(r) fraction*target_gammap(r);
    end
    LX.beta_value = beta_value;
end
