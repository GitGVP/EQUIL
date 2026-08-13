function state = equil_variational_state(L, LX, x)
%EQUIL_VARIATIONAL_STATE Evaluate global profiles and first derivatives.

    if LX.eps_val <= 0
        error('eps_val must be positive.');
    end
    if numel(LX.Sbc) ~= L.P.Ns
        error('Sbc must have one entry for each of the %d S_m profiles.', ...
              L.P.Ns);
    end
    if numel(LX.Abc) ~= L.P.Na
        error('Abc must have one entry for each of the %d A_m profiles.', ...
              L.P.Na);
    end
    if L.P.vertical_shift && (~isscalar(LX.Zbc) || ~isfinite(LX.Zbc))
        error('Zbc must be one finite scalar when vertical_shift is enabled.');
    end

    state.r = L.r_q(:);
    state.omega = L.omega(:).';
    state.epsilon = LX.eps_val;
    state.Sbc = LX.Sbc(:);
    state.Abc = LX.Abc(:);

    nprofiles = 3+L.P.Ns+L.P.Na+double(L.P.vertical_shift);
    values = cell(nprofiles, 1);
    derivatives = cell(nprofiles, 1);
    edge_values = cell(nprofiles, 1);
    edge_derivatives = cell(nprofiles, 1);
    for profile = 1:nprofiles
        first = L.profile_starts(profile);
        count = L.profile_lengths(profile);
        coeff = x(first:first+count-1);
        values{profile} = L.B0{profile}*coeff;
        derivatives{profile} = L.B1{profile}*coeff;
        edge_values{profile} = L.B0_edge{profile}*coeff;
        edge_derivatives{profile} = L.B1_edge{profile}*coeff;
    end

    state.t2 = values{1};
    state.t2r = derivatives{1};
    state.delta = values{2};
    state.deltar = derivatives{2};
    state.P = values{3};
    state.Pr = derivatives{3};

    state.S = zeros(L.Nq, 1, L.P.Ns);
    state.Sr = zeros(L.Nq, 1, L.P.Ns);
    edge_S = zeros(1, 1, L.P.Ns);
    edge_Sr = zeros(1, 1, L.P.Ns);
    for is = 1:L.P.Ns
        state.S(:,:,is) = values{3+is}+L.Sbc0{is}*state.Sbc(is);
        state.Sr(:,:,is) = derivatives{3+is}+L.Sbc1{is}*state.Sbc(is);
        edge_S(:,:,is) = edge_values{3+is} ...
            +L.Sbc0_edge(is)*state.Sbc(is);
        edge_Sr(:,:,is) = edge_derivatives{3+is} ...
            +L.Sbc1_edge(is)*state.Sbc(is);
    end

    state.A = zeros(L.Nq, 1, L.P.Na);
    state.Ar = zeros(L.Nq, 1, L.P.Na);
    edge_A = zeros(1, 1, L.P.Na);
    edge_Ar = zeros(1, 1, L.P.Na);
    for ia = 1:L.P.Na
        profile = 3+L.P.Ns+ia;
        state.A(:,:,ia) = values{profile}+L.Abc0{ia}*state.Abc(ia);
        state.Ar(:,:,ia) = derivatives{profile}+L.Abc1{ia}*state.Abc(ia);
        edge_A(:,:,ia) = edge_values{profile} ...
            +L.Abc0_edge(ia)*state.Abc(ia);
        edge_Ar(:,:,ia) = edge_derivatives{profile} ...
            +L.Abc1_edge(ia)*state.Abc(ia);
    end

    state.Zshift = zeros(L.Nq,1);
    state.Zshiftr = zeros(L.Nq,1);
    edge_Zshift = 0;
    edge_Zshiftr = 0;
    if L.P.vertical_shift
        profile = 4+L.P.Ns+L.P.Na;
        state.Zshift = values{profile}+L.Zbc0*LX.Zbc;
        state.Zshiftr = derivatives{profile}+L.Zbc1*LX.Zbc;
        edge_Zshift = edge_values{profile}+L.Zbc0_edge*LX.Zbc;
        edge_Zshiftr = edge_derivatives{profile}+L.Zbc1_edge*LX.Zbc;
    end

    axis_pressure = equil_variational_pressure(L, LX, 0, 1, 1);
    state.sigma0 = axis_pressure.PB;
    state.a0 = 1-state.sigma0;
    if ~isfinite(state.a0) || state.a0 <= L.P.min_one_minus_sigma
        error('Invalid on-axis anisotropy: 1-sigma0 = %.4e.', state.a0);
    end

    epsilon = state.epsilon;
    r = state.r;
    omega = state.omega;
    c1 = cos(omega);
    s1 = sin(omega);
    modes = reshape(1:L.P.Ns, 1, 1, []);
    cm = cos(modes.*omega);
    sm = sin(modes.*omega);

    state.R = 1 + epsilon*r.*c1 - epsilon^2*state.delta ...
        + epsilon^2*sum(state.S.*cm, 3) + epsilon^3*state.P.*c1;
    state.Z = epsilon*r.*s1 - epsilon^2*sum(state.S.*sm, 3) ...
        + epsilon^3*state.P.*s1;
    state.Rr = epsilon*c1 - epsilon^2*state.deltar ...
        + epsilon^2*sum(state.Sr.*cm, 3) + epsilon^3*state.Pr.*c1;
    state.Zr = epsilon*s1 - epsilon^2*sum(state.Sr.*sm, 3) ...
        + epsilon^3*state.Pr.*s1;
    state.Rw = -epsilon*r.*s1 - epsilon^2*sum(modes.*state.S.*sm, 3) ...
        - epsilon^3*state.P.*s1;
    state.Zw = epsilon*r.*c1 - epsilon^2*sum(modes.*state.S.*cm, 3) ...
        + epsilon^3*state.P.*c1;
    if L.P.Na > 0
        angular_modes = reshape(L.P.A_modes-1, 1, 1, []);
        sam = sin(angular_modes.*omega);
        cam = cos(angular_modes.*omega);
        state.R = state.R+epsilon^2*sum(state.A.*sam,3);
        state.Z = state.Z+epsilon^2*sum(state.A.*cam,3);
        state.Rr = state.Rr+epsilon^2*sum(state.Ar.*sam,3);
        state.Zr = state.Zr+epsilon^2*sum(state.Ar.*cam,3);
        state.Rw = state.Rw+epsilon^2*sum( ...
            angular_modes.*state.A.*cam,3);
        state.Zw = state.Zw-epsilon^2*sum( ...
            angular_modes.*state.A.*sam,3);
    end
    if L.P.vertical_shift
        state.Z = state.Z+epsilon^2*state.Zshift;
        state.Zr = state.Zr+epsilon^2*state.Zshiftr;
    end

    state.JoverR = state.Rr.*state.Zw-state.Rw.*state.Zr;
    state.J = state.R.*state.JoverR;
    state.goo = state.Rw.^2+state.Zw.^2;
    if any(state.R(:) <= 0) || any(state.J(:) <= 0)
        error('The trial state has non-positive R or an inverted Jacobian.');
    end

    state.q = LX.qfun(r);
    if any(~isfinite(state.q)) || any(state.q <= 0)
        error('q(r) must be finite and positive.');
    end
    state.T = state.a0 + epsilon^2*state.t2;
    state.Tr = epsilon^2*state.t2r;
    state.psir = epsilon*r.*state.T./(state.q*state.a0);

    % J is formed with d/d(rhat).  The physical-normalized radial
    % Jacobian is J/epsilon.
    state.Bp2 = (epsilon*state.psir).^2.*state.goo./state.J.^2;

    state.edge = make_geometry(1, state.omega, epsilon, state.a0, ...
        edge_values{1}, edge_derivatives{1}, ...
        edge_values{2}, edge_derivatives{2}, ...
        edge_values{3}, edge_derivatives{3}, edge_S, edge_Sr, ...
        edge_A, edge_Ar, L.P.A_modes,edge_Zshift,edge_Zshiftr, ...
        LX.qfun(1));
    state.edge.sigma0 = state.sigma0;
    state.edge.a0 = state.a0;
    state.edge.epsilon = epsilon;
end

function edge = make_geometry(r, omega, epsilon, a0, t2, t2r, ...
        delta, deltar, P, Pr, S, Sr, A, Ar, A_modes, ...
        Zshift,Zshiftr,q)
    Ns = size(S, 3);
    modes = reshape(1:Ns, 1, 1, []);
    cm = cos(modes.*omega);
    sm = sin(modes.*omega);
    c1 = cos(omega);
    s1 = sin(omega);

    edge.r = r;
    edge.omega = omega;
    edge.t2 = t2;
    edge.t2r = t2r;
    edge.delta = delta;
    edge.deltar = deltar;
    edge.P = P;
    edge.Pr = Pr;
    edge.S = S;
    edge.Sr = Sr;
    edge.A = A;
    edge.Ar = Ar;
    edge.Zshift = Zshift;
    edge.Zshiftr = Zshiftr;
    edge.R = 1+epsilon*r.*c1-epsilon^2*delta ...
        +epsilon^2*sum(S.*cm,3)+epsilon^3*P.*c1;
    edge.Z = epsilon*r.*s1-epsilon^2*sum(S.*sm,3)+epsilon^3*P.*s1;
    edge.Rr = epsilon*c1-epsilon^2*deltar ...
        +epsilon^2*sum(Sr.*cm,3)+epsilon^3*Pr.*c1;
    edge.Zr = epsilon*s1-epsilon^2*sum(Sr.*sm,3)+epsilon^3*Pr.*s1;
    edge.Rw = -epsilon*r.*s1-epsilon^2*sum(modes.*S.*sm,3) ...
        -epsilon^3*P.*s1;
    edge.Zw = epsilon*r.*c1-epsilon^2*sum(modes.*S.*cm,3) ...
        +epsilon^3*P.*c1;
    Na = size(A,3);
    if Na > 0
        angular_modes = reshape(A_modes-1,1,1,[]);
        sam = sin(angular_modes.*omega);
        cam = cos(angular_modes.*omega);
        edge.R = edge.R+epsilon^2*sum(A.*sam,3);
        edge.Z = edge.Z+epsilon^2*sum(A.*cam,3);
        edge.Rr = edge.Rr+epsilon^2*sum(Ar.*sam,3);
        edge.Zr = edge.Zr+epsilon^2*sum(Ar.*cam,3);
        edge.Rw = edge.Rw+epsilon^2*sum(angular_modes.*A.*cam,3);
        edge.Zw = edge.Zw-epsilon^2*sum(angular_modes.*A.*sam,3);
    end
    edge.Z = edge.Z+epsilon^2*Zshift;
    edge.Zr = edge.Zr+epsilon^2*Zshiftr;
    edge.JoverR = edge.Rr.*edge.Zw-edge.Rw.*edge.Zr;
    edge.J = edge.R.*edge.JoverR;
    edge.goo = edge.Rw.^2+edge.Zw.^2;
    edge.q = q;
    edge.T = a0+epsilon^2*t2;
    edge.Tr = epsilon^2*t2r;
    edge.psir = epsilon*r.*edge.T/(q*a0);
    edge.Bp2 = (epsilon*edge.psir).^2.*edge.goo./edge.J.^2;
    if any(edge.R(:) <= 0) || any(edge.J(:) <= 0)
        error('The trial state has a non-positive edge R or Jacobian.');
    end
end
