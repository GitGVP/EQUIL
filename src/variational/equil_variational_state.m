function state = equil_variational_state(L, LX, x)
%EQUIL_VARIATIONAL_STATE Evaluate global profiles and first derivatives.

    if LX.eps_val <= 0
        error('eps_val must be positive.');
    end
    if numel(LX.Sbc) ~= L.P.Ns
        error('Sbc must have one entry for each of the %d S_m profiles.', ...
              L.P.Ns);
    end
    if numel(LX.Vbc) ~= L.P.Nh
        error('Vbc must have one entry for each of the %d V_n profiles.', ...
              L.P.Nh);
    end

    state.r = L.r_q(:);
    state.omega = L.omega(:).';
    state.epsilon = LX.eps_val;
    state.Sbc = LX.Sbc(:);
    state.Vbc = LX.Vbc(:);

    nprofiles = 3+L.P.Ns+L.P.Nh;
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

    state.V = zeros(L.Nq, 1, L.P.Nh);
    state.Vr = zeros(L.Nq, 1, L.P.Nh);
    edge_V = zeros(1, 1, L.P.Nh);
    edge_Vr = zeros(1, 1, L.P.Nh);
    for iv = 1:L.P.Nh
        profile = 3+L.P.Ns+iv;
        state.V(:,:,iv) = values{profile}+L.Vbc0{iv}*state.Vbc(iv);
        state.Vr(:,:,iv) = derivatives{profile}+L.Vbc1{iv}*state.Vbc(iv);
        edge_V(:,:,iv) = edge_values{profile} ...
            +L.Vbc0_edge(iv)*state.Vbc(iv);
        edge_Vr(:,:,iv) = edge_derivatives{profile} ...
            +L.Vbc1_edge(iv)*state.Vbc(iv);
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
    if L.P.Nh > 0
        vmodes = reshape(1:L.P.Nh, 1, 1, []);
        svm = sin(vmodes.*omega);
        cvm = cos(vmodes.*omega);
        state.R = state.R+epsilon^2*sum(state.V.*svm,3);
        state.Z = state.Z-epsilon^2*sum(state.V.*cvm,3);
        state.Rr = state.Rr+epsilon^2*sum(state.Vr.*svm,3);
        state.Zr = state.Zr-epsilon^2*sum(state.Vr.*cvm,3);
        state.Rw = state.Rw+epsilon^2*sum(vmodes.*state.V.*cvm,3);
        state.Zw = state.Zw+epsilon^2*sum(vmodes.*state.V.*svm,3);
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
        edge_V, edge_Vr, LX.qfun(1));
    state.edge.sigma0 = state.sigma0;
    state.edge.a0 = state.a0;
    state.edge.epsilon = epsilon;
end

function edge = make_geometry(r, omega, epsilon, a0, t2, t2r, ...
        delta, deltar, P, Pr, S, Sr, V, Vr, q)
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
    edge.V = V;
    edge.Vr = Vr;
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
    Nh = size(V,3);
    if Nh > 0
        vmodes = reshape(1:Nh,1,1,[]);
        svm = sin(vmodes.*omega);
        cvm = cos(vmodes.*omega);
        edge.R = edge.R+epsilon^2*sum(V.*svm,3);
        edge.Z = edge.Z-epsilon^2*sum(V.*cvm,3);
        edge.Rr = edge.Rr+epsilon^2*sum(Vr.*svm,3);
        edge.Zr = edge.Zr-epsilon^2*sum(Vr.*cvm,3);
        edge.Rw = edge.Rw+epsilon^2*sum(vmodes.*V.*cvm,3);
        edge.Zw = edge.Zw+epsilon^2*sum(vmodes.*V.*svm,3);
    end
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
