function [L,X,Y,info] = meq_to_equil(meq_L,meq_Y,varargin)
%MEQ_TO_EQUIL Reproduce one retained MEQ/FBT equilibrium with EQUIL.
%   [L,X,Y,INFO] = MEQ_TO_EQUIL(MEQ_L,MEQ_Y) converts the last finite-q
%   domain, fits its pressure and q profiles, fixes its canonical boundary
%   harmonics, and iterates the EQUIL aspect-ratio scale.  The default
%   closure is isotropic.  Name/value options are intended for examples
%   and focused studies, not as a second solver interface.

    options = conversion_options(varargin{:});
    R0 = meq_Y.rA;
    Z0 = meq_Y.zA;
    assert(isscalar(R0) && isscalar(Z0), ...
        'meq_to_equil supports a single magnetic axis only.');

    fbt_psi = meq_L.pQ.^2;
    fbt_q = 1./meq_Y.iqQ;
    if isempty(options.q_axis) && isfield(meq_Y,'qA') ...
            && isscalar(meq_Y.qA) && isfinite(meq_Y.qA)
        options.q_axis = meq_Y.qA;
    end
    if ~isempty(options.q_axis)
        fbt_q(1) = options.q_axis;
    end
    finite_q = isfinite(fbt_q) & fbt_q > 0;
    assert(any(finite_q),'MEQ did not return a finite positive q profile.');
    last_finite = max(fbt_psi(finite_q));
    if ischar(options.psiN)
        psiN = last_finite;
    else
        psiN = options.psiN;
    end
    assert(psiN > 0 && psiN <= last_finite+10*eps, ...
        'Requested psiN exceeds the last finite-q MEQ surface.');

    fbt_R = [R0*ones(size(meq_Y.rq,1),1),meq_Y.rq];
    fbt_Z = [Z0*ones(size(meq_Y.zq,1),1),meq_Y.zq];
    boundary_R = interp1(fbt_psi,fbt_R.',psiN,'pchip').';
    boundary_Z = interp1(fbt_psi,fbt_Z.',psiN,'pchip').';
    [A_modes,classification] = choose_A_modes( ...
        meq_L.oq,boundary_R,boundary_Z,options);
    [Sraw,boundary, Araw] = get_Sbc( ...
        meq_L.oq,boundary_R,boundary_Z,false,options.Ns, ...
        numel(A_modes),A_modes);

    gamma_minus1 = real(boundary.negative_complex(1))/R0;
    gamma0 = (mean(boundary.R-1i*boundary.Z)-(R0-1i*Z0))/R0;
    Sscale = Sraw(:)/R0;
    Ascale = Araw(:)/R0;
    Zscale = (mean(boundary.Z)-Z0)/R0;
    boundary_radius = mean(hypot(boundary_R-R0,boundary_Z-Z0));
    vertical_shift = abs(Zscale)*R0 > ...
        options.vertical_shift_tolerance*boundary_radius;

    B0 = options.B0;
    if isempty(B0)
        B0 = meq_Y.TQ(1)/R0;
    end
    assert(isfinite(B0) && B0 > 0,'The MEQ axis field must be positive.');

    solver_arguments = {'beta',options.initial_beta, ...
        'equation_of_state',options.equation_of_state, ...
        'gamma0',options.initial_beta*options.gamma_hat, ...
        'Ns',options.Ns,'A_modes',A_modes, ...
        'vertical_shift',vertical_shift,'m',options.m,'nq',options.nq, ...
        'spline_p',options.spline_p,'om_pts',options.om_pts, ...
        'debug',options.debug};
    [L,X] = equilVariationalSol(solver_arguments{:});
    epsilon = gamma_minus1;
    X.eps_val = epsilon;
    X.Sbc = Sscale/epsilon^2;
    X.Abc = Ascale/epsilon^2;
    X.Zbc = Zscale/epsilon^2;
    X = set_closure_profiles(X,options, ...
        X.kinetic_profiles.beta,X.kinetic_profiles.betap);
    target_Sbc = X.Sbc;
    target_Abc = X.Abc;
    target_Zbc = X.Zbc;
    if isempty(A_modes)
        try
            Y = equilVariationalY(L,X);
        catch
            [X,Y] = ramp_boundaries( ...
                L,X,target_Sbc,target_Abc,target_Zbc);
        end
    else
        [X,Y] = ramp_boundaries( ...
            L,X,target_Sbc,target_Abc,target_Zbc);
    end
    if ~Y.isconverged
        [X,Y] = ramp_boundaries( ...
            L,X,target_Sbc,target_Abc,target_Zbc);
    end

    epsilon_history = zeros(options.niter,1);
    F_epsilon_history = zeros(options.niter,1);
    F_R_history = zeros(options.niter,1);
    profile_fraction = zeros(options.niter,1);
    for iteration = 1:options.niter
        mu0 = 4e-7*pi;
        pressure_scale = epsilon^2*B0^2/mu0;
        mapped_psi = psiN*Y.psiN_q;
        q_values = interp1( ...
            fbt_psi(finite_q),fbt_q(finite_q),mapped_psi,'pchip');
        [qfun,qpfun] = fit_even_profile_safe( ...
            L.r_q,q_values,2:2:options.q_fit_order,fbt_q(1), ...
            0.25*min(q_values),2*max(q_values));
        pressure_values = interp1( ...
            fbt_psi,meq_Y.PQ,mapped_psi,'pchip')/pressure_scale;
        pressure_bound = max(abs(pressure_values));
        [betafun,betapfun] = fit_even_profile_safe( ...
            L.r_q,pressure_values,0:2:options.pressure_fit_order,0, ...
            -1e-6*pressure_bound,2*pressure_bound);

        new_Sbc = Sscale/epsilon^2;
        new_Abc = Ascale/epsilon^2;
        new_Zbc = Zscale/epsilon^2;
        X.x = rescale_boundaries( ...
            L,Y.x,X.Sbc,new_Sbc,X.Abc,new_Abc,X.Zbc,new_Zbc);
        X.Sbc = new_Sbc;
        X.Abc = new_Abc;
        X.Zbc = new_Zbc;
        X.eps_val = epsilon;
        X.local_B_guess = Y.local_B_quadrature;

        old_q = X.qfun;
        old_qp = X.qpfun;
        old_beta = X.kinetic_profiles.beta;
        old_betap = X.kinetic_profiles.betap;
        accepted = false;
        last_failure = 'no converged trial';
        for fraction = 2.^-(0:7)
            trial = X;
            trial.qfun = blend(old_q,qfun,fraction);
            trial.qpfun = blend(old_qp,qpfun,fraction);
            beta_trial = blend(old_beta,betafun,fraction);
            betap_trial = blend(old_betap,betapfun,fraction);
            trial = set_closure_profiles( ...
                trial,options,beta_trial,betap_trial);
            try
                candidate = equilVariationalY(L,trial);
                if candidate.isconverged
                    X = trial;
                    Y = candidate;
                    accepted = true;
                    break
                end
                last_failure = sprintf('|R|=%.3e (%s)', ...
                    norm(candidate.residual),candidate.convergence_reason);
            catch profile_exception
                if isempty(profile_exception.stack)
                    last_failure = profile_exception.message;
                else
                    last_failure = sprintf('%s (%s:%d)', ...
                        profile_exception.message, ...
                        profile_exception.stack(1).name, ...
                        profile_exception.stack(1).line);
                end
            end
        end
        assert(accepted, ...
            ['MEQ-to-EQUIL profile matching failed at iteration %d. ', ...
             'Last trial: %s'],iteration,last_failure);
        profile_fraction(iteration) = fraction;

        F_epsilon_history(iteration) = ...
            epsilon*(1+epsilon^2*Y.P(end))-gamma_minus1;
        F_R_history(iteration) = ...
            -epsilon^2*Y.delta(end)-real(gamma0);
        epsilon_target = boundary_epsilon( ...
            gamma_minus1,Y.P(end),epsilon);
        if fraction < 1
            epsilon_next = epsilon;
        else
            maximum_step = options.max_epsilon_step*epsilon;
            epsilon_next = epsilon+max(-maximum_step, ...
                min(maximum_step,epsilon_target-epsilon));
        end
        epsilon_history(iteration) = epsilon;
        if options.debug > 0
            fprintf(['MEQ->EQUIL %2d: eps %.9f, F_eps %.3e, ', ...
                'F_R %.3e, profile %.3g\n'],iteration,epsilon, ...
                F_epsilon_history(iteration),F_R_history(iteration),fraction)
        end
        epsilon = epsilon_next;
    end

    info.success = Y.isconverged && profile_fraction(end) == 1 ...
        && abs(F_epsilon_history(end)) < options.scale_tolerance;
    info.R0 = R0;
    info.Z0 = Z0;
    info.B0 = B0;
    info.psiN = psiN;
    info.fbt_psiN = fbt_psi;
    info.fbt_q = fbt_q;
    info.fbt_R = fbt_R;
    info.fbt_Z = fbt_Z;
    info.boundary = boundary;
    info.boundary_R = boundary_R;
    info.boundary_Z = boundary_Z;
    info.gamma_minus1 = gamma_minus1;
    info.gamma0 = gamma0;
    info.F_epsilon = F_epsilon_history(end);
    info.F_R = F_R_history(end);
    info.F_epsilon_history = F_epsilon_history;
    info.F_R_history = F_R_history;
    info.epsilon_history = epsilon_history;
    info.profile_fraction = profile_fraction;
    info.A_modes = A_modes;
    info.auto_A = classification;
    info.vertical_shift = vertical_shift;
    info.boundary_RMS_over_R0 = curve_rms( ...
        boundary_R,boundary_Z,R0*Y.RR(end,:),Z0+R0*Y.ZZ(end,:))/R0;
    assert(info.success, ...
        ['MEQ-to-EQUIL conversion did not reach the scale tolerance: ', ...
         'F_epsilon=%.3e.'],info.F_epsilon)
    fprintf(['MEQ->EQUIL matched psi_N<=%.4f: F_epsilon %.2e, ', ...
        'F_R %.2e, boundary RMS/R0 %.2e.\n'], ...
        psiN,info.F_epsilon,info.F_R,info.boundary_RMS_over_R0)
end


function options = conversion_options(varargin)
    options.psiN = 'auto';
    options.Ns = 3;
    options.A_modes = 'auto';
    options.niter = 12;
    options.m = 8;
    options.nq = 6;
    options.spline_p = 4;
    options.om_pts = 96;
    options.q_fit_order = 10;
    options.pressure_fit_order = 14;
    options.max_epsilon_step = 0.05;
    options.scale_tolerance = 1e-5;
    options.vertical_shift_tolerance = 0.02;
    options.initial_beta = 1e-2;
    options.equation_of_state = @isotropic;
    options.gamma_hat = 0;
    options.B0 = [];
    options.q_axis = [];
    options.debug = 0;
    assert(mod(numel(varargin),2) == 0, ...
        'MEQ conversion options must be name/value pairs.');
    for k = 1:2:numel(varargin)
        name = varargin{k};
        assert(isfield(options,name),'Unknown MEQ conversion option %s.',name)
        options.(name) = varargin{k+1};
    end
end


function [modes,diagnostics] = choose_A_modes(oq,R,Z,options)
    if ~ischar(options.A_modes)
        modes = options.A_modes(:).';
        diagnostics.relative_norm = NaN;
        return
    end
    assert(strcmp(options.A_modes,'auto'), ...
        'A_modes must be a numeric vector or ''auto''.');
    reliable_max = floor(numel(oq)/2)-2;
    candidates = 2:min(6,reliable_max+1);
    [~,~,coefficients] = get_Sbc( ...
        oq,R,Z,false,options.Ns,numel(candidates),candidates);
    radius = mean(hypot(R-mean(R),Z-mean(Z)));
    relative_norm = norm(coefficients)/max(radius,eps);
    diagnostics.relative_norm = relative_norm;
    diagnostics.candidates = candidates;
    diagnostics.coefficients = coefficients;
    if relative_norm < 5e-4
        modes = [];
        return
    end
    [energy,order] = sort(abs(coefficients).^2,'descend');
    count = find(cumsum(energy)/sum(energy) >= 0.95,1,'first');
    modes = sort(candidates(order(1:min(count,4))));
end


function X = set_closure_profiles(X,options,beta,betap)
    X.kinetic_profiles.beta = beta;
    X.kinetic_profiles.betap = betap;
    if strcmp(func2str(options.equation_of_state),'runaways')
        X.kinetic_profiles.gamma = @(r) options.gamma_hat*beta(r);
        X.kinetic_profiles.gammap = @(r) options.gamma_hat*betap(r);
    elseif ~strcmp(func2str(options.equation_of_state),'isotropic')
        error('meq_to_equil currently supports isotropic and runaways.');
    end
end


function [profile,derivative] = fit_even_profile_safe( ...
        r,values,orders,offset,lower_bound,upper_bound)
    for order = fliplr(orders)
        retained = orders(orders <= order);
        design = bsxfun(@power,r(:),retained);
        coefficients = design\(values(:)-offset);
        candidate = @(rr) reshape(offset+ ...
            bsxfun(@power,rr(:),retained)*coefficients,size(rr));
        positive = retained > 0;
        positive_orders = retained(positive);
        positive_orders = positive_orders(:);
        derivative_orders = positive_orders-1;
        positive_coefficients = coefficients(positive);
        positive_coefficients = positive_coefficients(:);
        derivative_coefficients = ...
            (derivative_orders+1).*positive_coefficients;
        candidate_derivative = @(rr) reshape( ...
            bsxfun(@power,rr(:),derivative_orders.') ...
            *derivative_coefficients,size(rr));
        test_values = candidate(linspace(0,1,401).');
        if all(isfinite(test_values)) && min(test_values) >= lower_bound ...
                && max(test_values) <= upper_bound
            profile = candidate;
            derivative = candidate_derivative;
            return
        end
    end
    error('Unable to construct a bounded even-polynomial profile fit.');
end


function epsilon = boundary_epsilon(gamma_minus1,Pedge,current)
    if abs(Pedge) < 100*eps
        epsilon = gamma_minus1;
        return
    end
    candidates = roots([Pedge,0,1,-gamma_minus1]);
    candidates = real(candidates( ...
        abs(imag(candidates)) < 1e-10 & real(candidates) > 0));
    assert(~isempty(candidates),'The conformal boundary scale has no root.');
    [~,index] = min(abs(candidates-current));
    epsilon = candidates(index);
end


function output = blend(old,new,fraction)
    if fraction == 1
        output = new;
    elseif fraction == 0
        output = old;
    else
        output = @(r) (1-fraction)*old(r)+fraction*new(r);
    end
end


function x = rescale_boundaries( ...
        L,x,old_S,new_S,old_A,new_A,old_Z,new_Z)
    for is = 1:L.P.Ns
        profile = 3+is;
        rows = L.profile_starts(profile)+ ...
            (0:L.profile_lengths(profile)-1);
        if old_S(is) ~= 0
            x(rows) = x(rows)*(new_S(is)/old_S(is));
        elseif new_S(is) ~= 0
            target = new_S(is)*L.r_q.^is-L.Sbc0{is}*new_S(is);
            x(rows) = L.B0{profile}\target;
        end
    end
    for ia = 1:L.P.Na
        profile = 3+L.P.Ns+ia;
        rows = L.profile_starts(profile)+ ...
            (0:L.profile_lengths(profile)-1);
        if old_A(ia) ~= 0
            x(rows) = x(rows)*(new_A(ia)/old_A(ia));
        elseif new_A(ia) ~= 0
            power = L.P.A_leading_powers(ia);
            target = new_A(ia)*L.r_q.^power-L.Abc0{ia}*new_A(ia);
            x(rows) = L.B0{profile}\target;
        end
    end
    if L.P.vertical_shift
        profile = 4+L.P.Ns+L.P.Na;
        rows = L.profile_starts(profile)+ ...
            (0:L.profile_lengths(profile)-1);
        if old_Z ~= 0
            x(rows) = x(rows)*(new_Z/old_Z);
        elseif new_Z ~= 0
            target = new_Z*L.r_q.^2-L.Zbc0*new_Z;
            x(rows) = L.B0{profile}\target;
        end
    end
end


function [X,Y] = ramp_boundaries(L,X,target_S,target_A,target_Z)
    X.Sbc = zeros(size(target_S));
    X.Abc = zeros(size(target_A));
    X.Zbc = 0;
    X.x = [];
    X.local_B_guess = [];
    Y = equilVariationalY(L,X);
    for component = 1:3
        if (component == 2 && ~L.P.vertical_shift) ...
                || (component == 3 && L.P.Na == 0)
            continue
        end
        fraction = 0;
        if component == 1
            step = 0.1;
        else
            step = 0.05;
        end
        while fraction < 1
            trial_fraction = min(1,fraction+step);
            new_S = target_S;
            new_A = target_A;
            new_Z = target_Z;
            if component == 1
                new_S = trial_fraction*target_S;
                new_A(:) = 0;
                new_Z = 0;
            elseif component == 2
                new_A(:) = 0;
                new_Z = trial_fraction*target_Z;
            else
                new_A = trial_fraction*target_A;
            end
            trial = X;
            trial.x = rescale_boundaries( ...
                L,Y.x,X.Sbc,new_S,X.Abc,new_A,X.Zbc,new_Z);
            trial.Sbc = new_S;
            trial.Abc = new_A;
            trial.Zbc = new_Z;
            trial.local_B_guess = Y.local_B_quadrature;
            converged = false;
            try
                candidate = equilVariationalY(L,trial);
                converged = candidate.isconverged;
            catch
                % Retry a shorter fixed-boundary step.
            end
            if converged
                fraction = trial_fraction;
                X = trial;
                Y = candidate;
                step = min(1.4*step,0.15);
            else
                step = step/2;
                assert(step >= 1e-3, ...
                    'Boundary continuation failed at fraction %.3f.', ...
                    fraction);
            end
        end
        if component == 1 && ~L.P.vertical_shift
            X.Zbc = target_Z;
        end
    end
end


function rms = curve_rms(R1,Z1,R2,Z2)
    if hypot(R1(end)-R1(1),Z1(end)-Z1(1)) < 1e-12
        R1 = R1(1:end-1); Z1 = Z1(1:end-1);
    end
    if hypot(R2(end)-R2(1),Z2(end)-Z2(1)) < 1e-12
        R2 = R2(1:end-1); Z2 = Z2(1:end-1);
    end
    count = 256;
    R1 = interpft(R1(:),count); Z1 = interpft(Z1(:),count);
    R2 = interpft(R2(:),count); Z2 = interpft(Z2(:),count);
    distance2 = bsxfun(@minus,R1,R2.').^2+ ...
        bsxfun(@minus,Z1,Z2.').^2;
    values = [min(distance2,[],2);min(distance2,[],1).'];
    rms = sqrt(mean(values));
end
