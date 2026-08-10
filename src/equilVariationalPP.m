function LY = equilVariationalPP(L, LX, LY, cache)
%EQUILVARIATIONALPP Plot-ready postprocessing for the variational solver.

    if nargin < 4 || isempty(cache)
        [~,cache] = equil_variational_residual(L,LX,LY.x,[]);
    end
    state = cache.state;
    epsilon = LX.eps_val;
    nprofiles = 3+L.P.Ns+L.P.Nh;

    qvalues = cell(nprofiles,1);
    qderivatives = cell(nprofiles,1);
    axis_values = cell(nprofiles,1);
    axis_derivatives = cell(nprofiles,1);
    edge_values = cell(nprofiles,1);
    edge_derivatives = cell(nprofiles,1);
    for profile = 1:nprofiles
        rows = block_rows(L,profile);
        coeff = LY.x(rows);
        qvalues{profile} = L.B0{profile}*coeff;
        qderivatives{profile} = L.B1{profile}*coeff;
        axis_values{profile} = full(L.B0_axis{profile}*coeff);
        axis_derivatives{profile} = full(L.B1_axis{profile}*coeff);
        edge_values{profile} = full(L.B0_edge{profile}*coeff);
        edge_derivatives{profile} = full(L.B1_edge{profile}*coeff);
    end

    r_plt = [0;L.r_q;1];
    omega_plt = [L.omega,2*pi];
    t2 = extend_profile(axis_values{1},qvalues{1},edge_values{1});
    t2r = extend_profile(axis_derivatives{1},qderivatives{1},edge_derivatives{1});
    delta = extend_profile(axis_values{2},qvalues{2},edge_values{2});
    deltar = extend_profile(axis_derivatives{2},qderivatives{2},edge_derivatives{2});
    gaugeP = extend_profile(axis_values{3},qvalues{3},edge_values{3});
    gaugePr = extend_profile(axis_derivatives{3},qderivatives{3},edge_derivatives{3});

    S = zeros(numel(r_plt),1,L.P.Ns);
    Sr = zeros(numel(r_plt),1,L.P.Ns);
    for is = 1:L.P.Ns
        profile = 3+is;
        axis_S = axis_values{profile}+L.Sbc0_axis(is)*LX.Sbc(is);
        axis_Sr = axis_derivatives{profile}+L.Sbc1_axis(is)*LX.Sbc(is);
        quad_S = qvalues{profile}+L.Sbc0{is}*LX.Sbc(is);
        quad_Sr = qderivatives{profile}+L.Sbc1{is}*LX.Sbc(is);
        edge_S = edge_values{profile}+L.Sbc0_edge(is)*LX.Sbc(is);
        edge_Sr = edge_derivatives{profile}+L.Sbc1_edge(is)*LX.Sbc(is);
        S(:,:,is) = extend_profile(axis_S,quad_S,edge_S);
        Sr(:,:,is) = extend_profile(axis_Sr,quad_Sr,edge_Sr);
    end

    V = zeros(numel(r_plt),1,L.P.Nh);
    Vr = zeros(numel(r_plt),1,L.P.Nh);
    for iv = 1:L.P.Nh
        profile = 3+L.P.Ns+iv;
        axis_V = axis_values{profile}+L.Vbc0_axis(iv)*LX.Vbc(iv);
        axis_Vr = axis_derivatives{profile}+L.Vbc1_axis(iv)*LX.Vbc(iv);
        quad_V = qvalues{profile}+L.Vbc0{iv}*LX.Vbc(iv);
        quad_Vr = qderivatives{profile}+L.Vbc1{iv}*LX.Vbc(iv);
        edge_V = edge_values{profile}+L.Vbc0_edge(iv)*LX.Vbc(iv);
        edge_Vr = edge_derivatives{profile}+L.Vbc1_edge(iv)*LX.Vbc(iv);
        V(:,:,iv) = extend_profile(axis_V,quad_V,edge_V);
        Vr(:,:,iv) = extend_profile(axis_Vr,quad_Vr,edge_Vr);
    end

    modes = reshape(1:L.P.Ns,1,1,[]);
    cm = cos(modes.*omega_plt);
    sm = sin(modes.*omega_plt);
    c1 = cos(omega_plt);
    s1 = sin(omega_plt);
    R = 1+epsilon*r_plt.*c1-epsilon^2*delta ...
        +epsilon^2*sum(S.*cm,3)+epsilon^3*gaugeP.*c1;
    Z = epsilon*r_plt.*s1-epsilon^2*sum(S.*sm,3) ...
        +epsilon^3*gaugeP.*s1;
    Rr = epsilon*c1-epsilon^2*deltar ...
        +epsilon^2*sum(Sr.*cm,3)+epsilon^3*gaugePr.*c1;
    Zr = epsilon*s1-epsilon^2*sum(Sr.*sm,3) ...
        +epsilon^3*gaugePr.*s1;
    Rw = -epsilon*r_plt.*s1-epsilon^2*sum(modes.*S.*sm,3) ...
        -epsilon^3*gaugeP.*s1;
    Zw = epsilon*r_plt.*c1-epsilon^2*sum(modes.*S.*cm,3) ...
        +epsilon^3*gaugeP.*c1;
    if L.P.Nh > 0
        vmodes = reshape(1:L.P.Nh,1,1,[]);
        svm = sin(vmodes.*omega_plt);
        cvm = cos(vmodes.*omega_plt);
        R = R+epsilon^2*sum(V.*svm,3);
        Z = Z-epsilon^2*sum(V.*cvm,3);
        Rr = Rr+epsilon^2*sum(Vr.*svm,3);
        Zr = Zr-epsilon^2*sum(Vr.*cvm,3);
        Rw = Rw+epsilon^2*sum(vmodes.*V.*cvm,3);
        Zw = Zw+epsilon^2*sum(vmodes.*V.*svm,3);
    end
    JoverR = Rr.*Zw-Rw.*Zr;
    J = R.*JoverR;
    goo = Rw.^2+Zw.^2;

    q = LX.qfun(r_plt);
    T = state.a0+epsilon^2*t2;
    psir = epsilon*r_plt.*T./(q*state.a0);
    Bp2 = zeros(size(R));
    Bp2(2:end,:) = (epsilon*psir(2:end)).^2.*goo(2:end,:)./J(2:end,:).^2;

    plot_state.r = r_plt(2:end);
    plot_state.R = R(2:end,:);
    plot_state.Bp2 = Bp2(2:end,:);
    plot_state.T = T(2:end);
    local = equil_variational_local_B(L,LX,plot_state,[]);
    B = [ones(1,numel(omega_plt));local.B];
    Bphi = [ones(1,numel(omega_plt));local.Bphi];
    sigma = [state.sigma0*ones(1,numel(omega_plt));local.sigma];

    axis_pressure = equil_variational_pressure(L,LX,0,1,1);
    axis_sigma_B = axis_pressure.PBB-axis_pressure.PB;
    GB = [(-axis_sigma_B+state.a0)*ones(1,numel(omega_plt));local.GB];
    Pi_parallel = [axis_pressure.Pi*ones(1,numel(omega_plt)); ...
                   local.pressure.Pi];
    Pi_perp = [axis_pressure.Pperp*ones(1,numel(omega_plt)); ...
               local.pressure.Pperp];

    psi = cumtrapz(r_plt,r_plt.*T./(q*state.a0));
    psiN = (psi-psi(1))/(psi(end)-psi(1));

    LY.r = L.r_q;
    LY.r_plt = r_plt;
    LY.omega = L.omega;
    LY.omega_plt = omega_plt;
    LY.t2 = t2;
    LY.t2p = t2r;
    LY.delta = delta;
    LY.deltap = deltar;
    LY.P = gaugeP;
    LY.Pp = gaugePr;
    LY.S = S;
    LY.Sp = Sr;
    LY.V = V;
    LY.Vp = Vr;
    LY.RR = R;
    LY.ZZ = Z;
    LY.J = J;
    LY.JoverR = JoverR;
    LY.goo = goo;
    LY.BB = B;
    LY.Bphi = Bphi;
    LY.Bp2 = Bp2;
    LY.Bphi2 = Bphi.^2;
    LY.B2 = B.^2;
    LY.sigma = sigma;
    LY.GB = GB;
    LY.firehose_margin = min(1-sigma(:));
    LY.mirror_margin = min(GB(:));
    LY.Pi_parallel = Pi_parallel;
    LY.Pi_perp = Pi_perp;
    LY.betapar = Pi_parallel/epsilon^2;
    LY.betaperp = Pi_perp/epsilon^2;
    LY.psi = psi;
    LY.psiN = psiN;
    LY.psiN_q = interp1(r_plt,psiN,L.r_q,'pchip');
    LY.T = T;
    LY.psir = psir;
    LY.q = q;
    LY.a0 = state.a0;
    LY.sigma0 = state.sigma0;
    LY.gauge_error = cache.fields.gauge;
    LY.local_B_residual = max(LY.local_B_residual,local.max_abs_G);
    scaled_J = J(2:end,:)./(epsilon^2*r_plt(2:end));
    LY.min_J = min(J(2:end,:),[],'all');
    LY.min_J_over_eps2r = min(scaled_J,[],'all');
    LY.min_edge_J_over_eps2 = min(J(end,:))/epsilon^2;

    toroidal_flux = mean(abs(Bphi).*JoverR,2);
    LY.rhotor = sqrt(max(toroidal_flux,0)/toroidal_flux(end));
    rmin = min(R,[],2);
    rmax = max(R,[],2);
    [zmin,izmin] = min(Z,[],2);
    [zmax,izmax] = max(Z,[],2);
    nr = numel(r_plt);
    rzmin = R(sub2ind(size(R),(1:nr).',izmin));
    rzmax = R(sub2ind(size(R),(1:nr).',izmax));
    LY.kappa = (zmax-zmin)./(rmax-rmin);
    LY.rgeom = (rmax+rmin)/2;
    LY.zgeom = (zmax+zmin)/2;
    LY.aminor = (rmax-rmin)/2;
    LY.deltal = (LY.rgeom-rzmin)./LY.aminor;
    LY.deltau = (LY.rgeom-rzmax)./LY.aminor;
    LY.deltatrig = (LY.deltal+LY.deltau)/2;
    LY.kappa(1) = 1;
    LY.deltal(1) = 0;
    LY.deltau(1) = 0;
    LY.deltatrig(1) = 0;
end

function value = extend_profile(axis_value,quadrature_value,edge_value)
    value = [axis_value;quadrature_value;edge_value];
end

function rows = block_rows(L,profile)
    first = L.profile_starts(profile);
    rows = first:first+L.profile_lengths(profile)-1;
end
