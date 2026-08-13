function [Sbc, diagnostics, Abc] = get_Sbc( ...
        oq, rq, zq, debug, Nm, Na, A_modes, angle_mode_limit, ...
        angle_initial_coefficients)
%GET_SBC Recover EQUIL S_m and A_m boundary harmonics from an FBT contour.
% Inputs:
%  oq,rq,zq  - samples on the same (increasing) oq grid (column vectors)
%  debug      - print Newton and Fourier diagnostics (default false)
%  Nm         - number of returned positive harmonics (default preserves
%               the historical output length)
%  Na         - number of returned A_m sine-quadrature harmonics
%  A_modes    - selected physical mode numbers m>=2 (default 2:Na+1)
%  angle_mode_limit - harmonics used to reconstruct omega (default is the
%               conservative cutoff of the supplied samples). A smaller
%               explicit value is useful for sensitivity tests on a dense
%               contour reconstructed from an interpolated flux grid.
%  angle_initial_coefficients - optional [sine;cosine] starting vector.
%               Coefficients are padded/truncated to angle_mode_limit.
% Outputs:
%  Sbc         - real positive Fourier coefficients of R-iZ
%  diagnostics - parameterization and reconstruction accuracy information
%  Abc         - A_m coefficients, in A_modes order. In the fixed EQUIL
%                gauge, S_m-i*A_m is the positive Fourier coefficient at
%                angular frequency m-1 of R-i*Z.
    if nargin <4
        debug=false;
    end

    oq = oq(:);
    rq = rq(:);
    zq = zq(:);
    if numel(rq) ~= numel(oq) || numel(zq) ~= numel(oq)
        error('oq, rq and zq must contain the same number of samples.');
    end

    Noriginal = numel(oq);
    Nfine = 10*Noriginal;
    rq_fine = interpft(rq, Nfine);
    zq_fine = interpft(zq, Nfine);
    oq_fine = linspace(0,2*pi,Nfine+1);
    oq_fine = oq_fine(1:end-1);

    maxIter = 80;tol = 1e-11;
    Nt = numel(oq_fine);
    historical_Nm = Nt/4-2;
    if nargin < 5 || isempty(Nm)
        Nm = historical_Nm;
    end
    if Nm < 1 || Nm ~= floor(Nm)
        error('Nm must be a positive integer.');
    end
    if nargin < 6 || isempty(Na)
        Na = 0;
    end
    if Na < 0 || Na ~= floor(Na)
        error('Na must be a non-negative integer.');
    end
    if nargin < 7 || isempty(A_modes)
        A_modes = 2:(Na+1);
    else
        A_modes = A_modes(:).';
    end
    if numel(A_modes) ~= Na || any(A_modes < 2) ...
            || any(A_modes ~= round(A_modes)) ...
            || numel(unique(A_modes)) ~= numel(A_modes)
        error('A_modes must contain Na distinct integers >= 2.');
    end
    A_angular_modes = A_modes-1;

    % Upsampling makes interpolation smoother but does not create new MEQ
    % information. Restrict the angle-reconstruction solve to modes that
    % are identifiable on the original angular grid.
    reliable_mode_max = floor(Noriginal/2)-2;
    if Nm > reliable_mode_max
        warning('get_Sbc:SymmetricModeLimit', ...
            ['Nm=%d exceeds the conservative cutoff %d for %d input ', ...
            'points; returned higher S coefficients are interpolation-', ...
            'dependent.'],Nm,reliable_mode_max,Noriginal);
    end
    if any(A_angular_modes > reliable_mode_max)
        error('get_Sbc:AsymmetricModeLimit', ...
            ['Requested A_%d (angular frequency %d) exceeds the ', ...
            'conservative ', ...
            'cutoff %d for %d input points. Supply a denser independently ', ...
            'sampled contour, for example with get_fbt_boundary.'], ...
            max(A_modes),max(A_angular_modes),reliable_mode_max,Noriginal);
    end
    if nargin < 8 || isempty(angle_mode_limit)
        angle_mode_limit = reliable_mode_max;
    end
    if angle_mode_limit < 1 || angle_mode_limit ~= round(angle_mode_limit) ...
            || angle_mode_limit > reliable_mode_max
        error('angle_mode_limit must be an integer from 1 through %d.', ...
            reliable_mode_max);
    end
    nA = angle_mode_limit;

    omega = linspace(0,2*pi,Nt+1).';
    omega = omega(1:end-1);

    % Transform from L.oq to correct orientation and start of angle
    theta_std = mod(pi - oq_fine(:), 2*pi);
    [theta_std, idx] = sort(theta_std);
    oq2 = theta_std;
    rq2 = rq_fine(idx);
    zq2 = zq_fine(idx);

    % Periodic extension. The old code appended rq(1), which generally
    % belongs to the opposite geometrical angle after sorting and creates
    % an artificial jump at 2*pi.
    theta_ext = [oq2(:); oq2(1)+2*pi];
    R_ext = [rq2(:); rq2(1)];
    Z_ext = [zq2(:); zq2(1)];

    % Spectral derivatives are consistent with the periodic interpft
    % representation and avoid a second-order differentiation error.
    Rp = periodic_derivative(rq2(:));
    Zp = periodic_derivative(zq2(:));
    Rp_ext = [Rp; Rp(1)]; Zp_ext = [Zp; Zp(1)];

    angle_sine = zeros(nA,1);
    angle_cosine = zeros(nA,1);
    if nargin >= 9 && ~isempty(angle_initial_coefficients)
        initial = angle_initial_coefficients(:);
        initial_half = floor(numel(initial)/2);
        sine_count = min(nA,initial_half);
        angle_sine(1:sine_count) = initial(1:sine_count);
        if numel(initial) >= 2*initial_half
            cosine_count = min(nA,numel(initial)-initial_half);
            angle_cosine(1:cosine_count) = ...
                initial(initial_half+(1:cosine_count));
        end
    end
    kvec = (1:nA).';
    S = sin(kvec * omega.'); % nA x Nt
    C = cos(kvec * omega.'); % nA x Nt
    asymmetric_parameterization = Na > 0;

    for it = 1:maxIter
        [res,Jac] = angle_system( ...
            angle_sine,angle_cosine,omega,S,C, ...
            theta_ext,R_ext,Z_ext,Rp_ext,Zp_ext, ...
            asymmetric_parameterization);
        if debug; fprintf('iter = %i, |res| = %.4e\n', it, norm(res)); end
        if norm(res) < tol, break; end
        delta = Jac \ res;
        damping = 1;
        accepted = false;
        for line_iteration = 1:14
            trial_sine = angle_sine-damping*delta(1:nA);
            trial_cosine = angle_cosine;
            if asymmetric_parameterization
                trial_cosine = angle_cosine-damping*delta(nA+1:end);
            end
            theta_derivative = 1+((kvec.*trial_sine).'*C).';
            if asymmetric_parameterization
                theta_derivative = theta_derivative- ...
                    ((kvec.*trial_cosine).'*S).';
            end
            if min(theta_derivative) > 0.05
                trial_residual = angle_system( ...
                    trial_sine,trial_cosine,omega,S,C, ...
                    theta_ext,R_ext,Z_ext,Rp_ext,Zp_ext, ...
                    asymmetric_parameterization);
                if norm(trial_residual) < norm(res)
                    angle_sine = trial_sine;
                    angle_cosine = trial_cosine;
                    accepted = true;
                    break
                end
            end
            damping = damping/2;
        end
        if ~accepted
            break
        end
    end

    % Refresh the samples after the last accepted Newton step.
    [res,~,Rint,Zint] = angle_system( ...
        angle_sine,angle_cosine,omega,S,C, ...
        theta_ext,R_ext,Z_ext,Rp_ext,Zp_ext, ...
        asymmetric_parameterization);

    if norm(res) >= tol
        warning('get_Sbc:NoConvergence', ...
            'Angle reconstruction ended with |res| = %.3e.',norm(res));
    end

    w = Rint-1i*Zint;
    Sbc = zeros(1,Nm);
    positive_complex = zeros(1,Nm);
    for m = 1:Nm
        positive_complex(m) = mean(w.*exp(-1i*m*omega));
        Sbc(m) = real(positive_complex(m));
        if debug
            negative = mean(w.*exp(1i*(m+1)*omega));
            fprintf(['m = %i, Sbc = %.6e, imag = %.3e, ', ...
                '(R-iZ)_{-%i} = %.3e%+.3ei\n'],m,Sbc(m), ...
                imag(positive_complex(m)),m+1,real(negative),imag(negative));
        end
    end

    negative_complex = zeros(1,reliable_mode_max);
    positive_reliable = zeros(1,reliable_mode_max);
    for m = 1:reliable_mode_max
        negative_complex(m) = mean(w.*exp(1i*m*omega));
        positive_reliable(m) = mean(w.*exp(-1i*m*omega));
    end
    negative_modes = negative_complex(2:end);
    c0 = mean(w);
    cminus1 = negative_complex(1);
    Abc = zeros(1,Na);
    for ia = 1:Na
        angular_mode = A_angular_modes(ia);
        Abc(ia) = -imag(positive_reliable(angular_mode));
    end
    reliable_positive = positive_complex(1:min(Nm,reliable_mode_max));
    w_reconstructed = c0+real(cminus1)*exp(-1i*omega);
    for m = 1:numel(reliable_positive)
        w_reconstructed = w_reconstructed+ ...
            real(reliable_positive(m))*exp(1i*m*omega);
    end
    for ia = 1:Na
        angular_mode = A_angular_modes(ia);
        if angular_mode <= reliable_mode_max
            w_reconstructed = w_reconstructed+ ...
                1i*imag(positive_reliable(angular_mode))* ...
                exp(1i*angular_mode*omega);
        end
    end

    omitted_modes = setdiff(1:reliable_mode_max,A_angular_modes);
    omitted_A = imag(positive_reliable(omitted_modes));
    if isempty(omitted_A)
        max_omitted_A = 0;
    else
        max_omitted_A = max(abs(omitted_A));
    end

    diagnostics.iterations = it;
    diagnostics.parameterization_residual = norm(res);
    diagnostics.reliable_mode_max = reliable_mode_max;
    diagnostics.angle_mode_limit = angle_mode_limit;
    diagnostics.asymmetric_parameterization = ...
        asymmetric_parameterization;
    diagnostics.omega = omega;
    diagnostics.R = Rint;
    diagnostics.Z = Zint;
    diagnostics.angle_coefficients = [angle_sine;angle_cosine];
    diagnostics.angle_sine_coefficients = angle_sine;
    diagnostics.angle_cosine_coefficients = angle_cosine;
    diagnostics.positive_complex = positive_complex;
    diagnostics.positive_reliable = positive_reliable;
    diagnostics.Abc = Abc;
    diagnostics.A_modes = A_modes;
    diagnostics.A_angular_modes = A_angular_modes;
    diagnostics.negative_complex = negative_complex;
    diagnostics.negative_modes = negative_modes;
    diagnostics.max_negative_mode = max(abs(negative_modes));
    diagnostics.max_omitted_A = max_omitted_A;
    diagnostics.max_positive_imaginary = ...
        max(abs(imag(positive_reliable)));
    diagnostics.reconstruction_relative_error = ...
        norm(w-w_reconstructed)/norm(w-c0);
end


function [residual,Jacobian,R,Z] = angle_system( ...
        angle_sine,angle_cosine,omega,S,C, ...
        theta_ext,R_ext,Z_ext,Rp_ext,Zp_ext,is_asymmetric)
%ANGLE_SYSTEM Fourier gauge conditions and their analytic Jacobian.
    nmode = numel(angle_sine);
    theta_model = omega+(angle_sine.'*S).';
    if is_asymmetric
        theta_model = theta_model+(angle_cosine.'*C).';
        basis = [S;C].';
        residual = zeros(2*nmode,1);
        Jacobian = zeros(2*nmode,2*nmode);
    else
        basis = S.';
        residual = zeros(nmode,1);
        Jacobian = zeros(nmode,nmode);
    end
    theta_wrapped = mod(theta_model,2*pi);
    R = interp1(theta_ext,R_ext,theta_wrapped,'pchip');
    Z = interp1(theta_ext,Z_ext,theta_wrapped,'pchip');
    Rp = interp1(theta_ext,Rp_ext,theta_wrapped,'pchip');
    Zp = interp1(theta_ext,Zp_ext,theta_wrapped,'pchip');

    % The negative-mode real parts are zero in the EQUIL gauge.
    for row = 1:nmode
        harmonic = row+1;
        cm = cos(harmonic*omega);
        sm = sin(harmonic*omega);
        residual(row) = mean(R.*cm+Z.*sm);
        integrand = Rp.*cm+Zp.*sm;
        Jacobian(row,:) = mean(bsxfun(@times,basis,integrand),1);
    end

    if ~is_asymmetric
        return
    end

    % The second half removes the imaginary parts of the same negative
    % harmonics. The retained positive-mode imaginary parts are -A_m.
    for mode = 1:nmode
        harmonic = mode+1;
        cm = cos(harmonic*omega);
        sm = sin(harmonic*omega);
        row = nmode+mode;
        residual(row) = mean(R.*sm-Z.*cm);
        integrand = Rp.*sm-Zp.*cm;
        Jacobian(row,:) = mean(bsxfun(@times,basis,integrand),1);
    end
end


function derivative = periodic_derivative(values)
%PERIODIC_DERIVATIVE FFT derivative on a uniform 2*pi-periodic grid.
    N = numel(values);
    if mod(N,2) == 0
        modes = [0:N/2-1, 0, -N/2+1:-1].';
    else
        half = (N-1)/2;
        modes = [0:half, -half:-1].';
    end
    derivative = real(ifft(1i*modes.*fft(values)));
end
