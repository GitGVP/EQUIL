function LX = equil_variational_set_beta(LX, beta_value)
%EQUIL_VARIATIONAL_SET_BETA Rescale the prescribed radial pressure profile.

    if ~isfield(LX,'beta_shape') || ~isfield(LX,'betap_shape')
        error('Automatic beta continuation requires beta_shape profiles.');
    end
    LX.beta_value = beta_value;
    beta_shape = LX.beta_shape;
    betap_shape = LX.betap_shape;
    LX.kinetic_profiles.beta = @(r) beta_value*beta_shape(r);
    LX.kinetic_profiles.betap = @(r) beta_value*betap_shape(r);
end
