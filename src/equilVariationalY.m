function LY = equilVariationalY(L, LX)
%EQUILVARIATIONALY Solve the condensed variational finite-element system.

    has_initial_state = isfield(LX, 'x') && ~isempty(LX.x);
    can_ramp_beta = isfield(LX,'beta_value') ...
        && (~isfield(LX,'pressure_map') || isempty(LX.pressure_map));
    if ~has_initial_state && L.P.beta_continuation && can_ramp_beta ...
            && LX.beta_value ~= 0
        LY = adaptive_beta_continuation(L,LX);
        return
    end

    if ~has_initial_state
        x = equil_variational_x0(L, LX);
    else
        x = LX.x(:);
    end
    if numel(x) ~= L.total_dofs
        error('Initial state has %d entries; expected %d.', ...
              numel(x), L.total_dofs);
    end

    Bguess = [];
    if isfield(LX,'local_B_guess') && ~isempty(LX.local_B_guess)
        Bguess = LX.local_B_guess;
    end
    if L.P.do_ana
        LY = equil_ana(L,LX);
    else
        LY = struct;
    end
    hard_limit = 3*L.P.nk;
    LY.res_norms = zeros(hard_limit,1);
    LY.isconverged = false;
    LY.convergence_reason = 'maximum iterations';
    effective_tol = L.P.NLtol;
    continuation_trial = isfield(LX,'continuation_trial') ...
        && LX.continuation_trial;
    iteration = 0;
    iteration_limit = L.P.nk;

    while iteration < iteration_limit
        iteration = iteration+1;
        [residual, cache] = equil_variational_residual(L, LX, x, Bguess);
        res_norm = norm(residual);
        LY.res_norms(iteration) = res_norm;

        if L.P.debug > 0
            fprintf('Variational iter %d: |R| = %.4e, max|G| = %.3e\n', ...
                    iteration, res_norm, cache.local.max_abs_G);
        end
        if res_norm < effective_tol
            LY.isconverged = true;
            LY.convergence_reason = 'residual';
            Bguess = cache.local.B;
            break
        end
        if continuation_trial && iteration >= 8
            recent = LY.res_norms(iteration-7:iteration);
            if res_norm > 1e4*effective_tol ...
                    && res_norm > 0.97*recent(1)
                LY.convergence_reason = 'stagnated continuation attempt';
                break
            end
        end

        jacobian = equil_variational_jacobian(L, LX, x, cache);
        update = jacobian \ residual;
        if res_norm < max(effective_tol,L.P.NLstepTol) ...
                && norm(update) < L.P.NLstepTol
            LY.isconverged = true;
            LY.convergence_reason = 'roundoff-scale update';
            Bguess = cache.local.B;
            break
        end
        alpha = L.P.damping;
        accepted = false;

        % A short residual backtrack is sufficient for this proof of
        % principle and also rejects inverted geometries/local B failures.
        for trial = 1:12
            x_trial = x - alpha*update;
            try
                [res_trial, trial_cache] = ...
                    equil_variational_residual(L, LX, x_trial, cache.local.B);
                if norm(res_trial) < res_norm
                    accepted = true;
                    break
                end
            catch
                % Try the next, smaller step.
            end
            alpha = alpha/2;
        end
        if ~accepted
            LY.convergence_reason = 'line search failed';
            break
        end

        x = x_trial;
        Bguess = trial_cache.local.B;
        if L.P.debug > 0
            fprintf('  |dx| = %.3e, alpha = %.3g\n', norm(update), alpha);
        end
        if iteration == iteration_limit && iteration_limit < hard_limit
            extra = predicted_extra_iterations( ...
                LY.res_norms(1:iteration),effective_tol,hard_limit-iteration);
            iteration_limit = iteration+extra;
        end
    end

    LY.res_norms = LY.res_norms(1:iteration);
    LY.x = x;
    [LY.residual, final_cache] = ...
        equil_variational_residual(L, LX, x, Bguess);
    LY.local_B_residual = final_cache.local.max_abs_G;
    LY.local_B_quadrature = final_cache.local.B;
    LY.effective_NLtol = effective_tol;
    nprofiles = 3+L.P.Ns+L.P.Nh;
    LY.residual_block_norms = zeros(nprofiles,1);
    LY.residual_block_names = cell(nprofiles,1);
    LY.residual_block_names(1:3) = {'t2','Delta','P'};
    for profile = 1:nprofiles
        rows = L.profile_starts(profile) ...
            +(0:L.profile_lengths(profile)-1);
        LY.residual_block_norms(profile) = norm(LY.residual(rows));
        if profile > 3 && profile <= 3+L.P.Ns
            LY.residual_block_names{profile} = sprintf('S%d',profile-2);
        elseif profile > 3+L.P.Ns
            LY.residual_block_names{profile} = ...
                sprintf('V%d',profile-(2+L.P.Ns));
        end
    end
    if norm(LY.residual) < effective_tol
        LY.isconverged = true;
        LY.convergence_reason = 'residual';
    end
    LY = equilVariationalPP(L, LX, LY, final_cache);
end

function LY = adaptive_beta_continuation(L,LX)
% Seed the nonlinear branch at low beta, then bisect only failed jumps.

    target_beta = LX.beta_value;
    stage_input = equil_variational_set_beta(LX,0);
    stage_input.x = equil_variational_x0(L,stage_input);

    accepted_betas = zeros(0,1);
    stage_norms = cell(0,1);
    stage_iterations = zeros(0,1);
    attempted_betas = zeros(0,1);
    attempted_epsilons = zeros(0,1);
    attempt_success = false(0,1);
    attempt_messages = cell(0,1);

    [ok,stage_solution,message] = attempt_beta_stage( ...
        L,stage_input,0);
    attempted_betas(end+1,1) = 0;
    attempted_epsilons(end+1,1) = stage_input.eps_val;
    attempt_success(end+1,1) = ok;
    attempt_messages{end+1,1} = message;
    if ~ok && stage_input.eps_val > 1e-3
        seed_input = stage_input;
        seed_input.eps_val = 1e-3;
        seed_input.x = equil_variational_x0(L,seed_input);
        [seed_ok,seed_solution,seed_message] = attempt_beta_stage( ...
            L,seed_input,0);
        attempted_betas(end+1,1) = 0;
        attempted_epsilons(end+1,1) = seed_input.eps_val;
        attempt_success(end+1,1) = seed_ok;
        attempt_messages{end+1,1} = seed_message;
        if seed_ok
            stage_input.x = seed_solution.x;
            stage_input.local_B_guess = seed_solution.local_B_quadrature;
            [ok,stage_solution,message] = attempt_beta_stage( ...
                L,stage_input,0);
            attempted_betas(end+1,1) = 0;
            attempted_epsilons(end+1,1) = stage_input.eps_val;
            attempt_success(end+1,1) = ok;
            attempt_messages{end+1,1} = message;
        else
            message = seed_message;
        end
    end
    if ~ok
        error('Variational beta=0 seed failed: %s',message);
    end
    [accepted_betas,stage_norms,stage_iterations] = ...
        append_accepted_stage(accepted_betas,stage_norms, ...
                              stage_iterations,0,stage_solution);

    current_beta = 0;
    failed_attempts = 0;
    low_beta = min(abs(target_beta),1);
    trial_beta = sign(target_beta)*low_beta;
    while current_beta ~= target_beta
        stage_input = equil_variational_set_beta(LX,trial_beta);
        stage_input.x = stage_solution.x;
        stage_input.local_B_guess = stage_solution.local_B_quadrature;
        [ok,trial_solution,message] = attempt_beta_stage( ...
            L,stage_input,trial_beta);
        attempted_betas(end+1,1) = trial_beta; %#ok<AGROW>
        attempted_epsilons(end+1,1) = stage_input.eps_val; %#ok<AGROW>
        attempt_success(end+1,1) = ok; %#ok<AGROW>
        attempt_messages{end+1,1} = message; %#ok<AGROW>

        if ok
            current_beta = trial_beta;
            stage_solution = trial_solution;
            [accepted_betas,stage_norms,stage_iterations] = ...
                append_accepted_stage(accepted_betas,stage_norms, ...
                    stage_iterations,current_beta,stage_solution);
            trial_beta = target_beta;
        else
            failed_attempts = failed_attempts+1;
            beta_step = 0.5*(trial_beta-current_beta);
            if failed_attempts > 8 || abs(beta_step) < 1e-4
                error(['Adaptive variational continuation could not cross ', ...
                    'beta %.6g -> %.6g after %d failed attempts. Last ', ...
                    'failure: %s'],current_beta,target_beta, ...
                    failed_attempts,message);
            end
            trial_beta = current_beta+beta_step;
        end
    end

    LY = stage_solution;
    LY.beta_continuation = accepted_betas.';
    LY.continuation_res_norms = stage_norms;
    LY.continuation_iterations = stage_iterations.';
    LY.beta_continuation_attempts = attempted_betas.';
    LY.beta_continuation_attempt_epsilons = attempted_epsilons.';
    LY.beta_continuation_attempt_success = attempt_success.';
    LY.beta_continuation_attempt_messages = attempt_messages;
    LY = add_analytical_output(L,LX,LY);
end

function [ok,stage_solution,message] = ...
        attempt_beta_stage(L,stage_input,beta_value)
    if L.P.debug > 0
        fprintf('Variational continuation attempt beta = %.6g at epsilon = %.4g\n', ...
                beta_value,stage_input.eps_val);
    end
    ok = false;
    stage_solution = [];
    message = '';
    try
        Ltrial = L;
        Ltrial.P.do_ana = false;
        Ltrial.P.do_shift_NLO = false;
        stage_input.continuation_trial = true;
        stage_solution = equilVariationalY(Ltrial,stage_input);
        ok = stage_solution.isconverged;
        if ~ok
            message = sprintf('|R| = %.3e (%s)', ...
                norm(stage_solution.residual), ...
                stage_solution.convergence_reason);
        end
    catch exception
        message = exception.message;
    end
    if L.P.debug > 0
        if ok
            fprintf('  accepted beta = %.6g in %d Newton iterations\n', ...
                    beta_value,numel(stage_solution.res_norms));
        else
            fprintf('  rejected beta = %.6g: %s\n',beta_value,message);
        end
    end
end

function extra = predicted_extra_iterations(history,tolerance,available)
    extra = 0;
    count = min(8,numel(history));
    values = history(end-count+1:end);
    slope = polyfit((1:count).',log(values),1);
    if slope(1) < 0 && all(diff(values) < 0)
        estimate = ceil(log(tolerance/values(end))/slope(1))+2;
        if estimate > 0
            extra = min(estimate,available);
        end
    elseif values(end) < 100*tolerance
        extra = min(10,available);
    end
end

function [betas,norms,iterations] = append_accepted_stage( ...
        betas,norms,iterations,beta_value,stage_solution)
    betas(end+1,1) = beta_value;
    norms{end+1,1} = stage_solution.res_norms;
    iterations(end+1,1) = numel(stage_solution.res_norms);
end

function LY = add_analytical_output(L,LX,LY)
% Analytical LO/NLO profiles are beta-stage independent and computed once.
    if ~L.P.do_ana
        return
    end
    analytical = equil_ana(L,LX);
    names = fieldnames(analytical);
    for k = 1:numel(names)
        LY.(names{k}) = analytical.(names{k});
    end
end
