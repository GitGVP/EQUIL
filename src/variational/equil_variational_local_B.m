function local = equil_variational_local_B(L, LX, state, Bguess)
%EQUIL_VARIATIONAL_LOCAL_B Eliminate B by a vectorized scalar Newton solve.

    % For the two isotropic closures sigma is identically zero, so the
    % local equation has the exact solution below.  Apart from being
    % cheaper, this makes it explicit that an isotropic run cannot fail
    % because of a numerical local-B residual.
    native_pressure = ~isfield(LX,'pressure_map') || isempty(LX.pressure_map);
    eos_name = func2str(L.P.equation_of_state);
    if native_pressure && any(strcmp(eos_name, ...
            {'isotropic','isotropic_rotating'}))
        B = sqrt(state.Bp2+(state.T./state.R).^2);
        constraint = equil_variational_B_constraint(L,LX,state,B);
        local = constraint;
        local.B = B;
        local.iterations = 0;
        local.max_abs_G = max(abs(constraint.G(:)));
        return
    end

    lower = sqrt(max(state.Bp2, 0));
    if nargin < 4 || isempty(Bguess) || ~isequal(size(Bguess), size(state.Bp2))
        B = sqrt(state.Bp2 + (state.T./state.R).^2);
    else
        B = Bguess;
    end
    B = max(B, lower.*(1+1e-10) + 1e-12);

    converged = false;
    for iteration = 1:L.P.local_B_maxit
        constraint = equil_variational_B_constraint(L, LX, state, B);
        if max(abs(constraint.G(:))) < L.P.local_B_tol
            converged = true;
            break
        end
        if any(abs(constraint.GB(:)) < L.P.min_GB)
            error('Local B condensation encountered |G_B| below the limit.');
        end
        step = constraint.G./constraint.GB;
        current_norm = max(abs(constraint.G(:)));
        alpha = 1;
        accepted = false;
        for trial = 1:12
            candidate = B-alpha*step;
            valid = all(isfinite(candidate(:))) ...
                && all(candidate(:) > lower(:).*(1+1e-10)+1e-12);
            if valid
                try
                    candidate_constraint = ...
                        equil_variational_B_constraint(L,LX,state,candidate);
                    candidate_norm = max(abs(candidate_constraint.G(:)));
                    if isfinite(candidate_norm) && candidate_norm < current_norm
                        accepted = true;
                        break
                    end
                catch
                    % Reduce the scalar-Newton step.
                end
            end
            alpha = alpha/2;
        end
        if ~accepted
            error('Local B Newton backtrack failed; max|G| = %.3e.', ...
                  current_norm);
        end
        B = candidate;
    end

    constraint = equil_variational_B_constraint(L, LX, state, B);
    max_abs_G = max(abs(constraint.G(:)));
    if ~converged && max_abs_G >= 10*L.P.local_B_tol
        error('Local B solve did not converge; max|G| = %.3e.', max_abs_G);
    end
    if any(1-constraint.sigma(:) <= L.P.min_one_minus_sigma)
        error('Local B solve reached 1-sigma <= prescribed minimum.');
    end
    if any(constraint.GB(:) <= L.P.min_GB)
        error('Local B solution reached the mirror boundary G_B <= limit.');
    end

    local = constraint;
    local.B = B;
    local.iterations = iteration;
    local.max_abs_G = max_abs_G;
end
