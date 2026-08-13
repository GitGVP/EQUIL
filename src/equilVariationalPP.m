function LY = equilVariationalPP(L, LX, LY, cache)
%EQUILVARIATIONALPP Plot-ready postprocessing for the variational solver.

    if nargin < 4 || isempty(cache)
        [~,cache] = equil_variational_residual(L,LX,LY.x,[]);
    end
    state = cache.state;
    epsilon = LX.eps_val;
    nprofiles = 3+L.P.Ns+L.P.Na+double(L.P.vertical_shift);

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

    A = zeros(numel(r_plt),1,L.P.Na);
    Ar = zeros(numel(r_plt),1,L.P.Na);
    for ia = 1:L.P.Na
        profile = 3+L.P.Ns+ia;
        axis_A = axis_values{profile}+L.Abc0_axis(ia)*LX.Abc(ia);
        axis_Ar = axis_derivatives{profile}+L.Abc1_axis(ia)*LX.Abc(ia);
        quad_A = qvalues{profile}+L.Abc0{ia}*LX.Abc(ia);
        quad_Ar = qderivatives{profile}+L.Abc1{ia}*LX.Abc(ia);
        edge_A = edge_values{profile}+L.Abc0_edge(ia)*LX.Abc(ia);
        edge_Ar = edge_derivatives{profile}+L.Abc1_edge(ia)*LX.Abc(ia);
        A(:,:,ia) = extend_profile(axis_A,quad_A,edge_A);
        Ar(:,:,ia) = extend_profile(axis_Ar,quad_Ar,edge_Ar);
    end

    Zshift = zeros(numel(r_plt),1);
    Zshiftr = zeros(numel(r_plt),1);
    if L.P.vertical_shift
        profile = 4+L.P.Ns+L.P.Na;
        axis_Zshift = axis_values{profile}+L.Zbc0_axis*LX.Zbc;
        axis_Zshiftr = axis_derivatives{profile}+L.Zbc1_axis*LX.Zbc;
        quad_Zshift = qvalues{profile}+L.Zbc0*LX.Zbc;
        quad_Zshiftr = qderivatives{profile}+L.Zbc1*LX.Zbc;
        edge_Zshift = edge_values{profile}+L.Zbc0_edge*LX.Zbc;
        edge_Zshiftr = edge_derivatives{profile}+L.Zbc1_edge*LX.Zbc;
        Zshift = extend_profile(axis_Zshift,quad_Zshift,edge_Zshift);
        Zshiftr = extend_profile( ...
            axis_Zshiftr,quad_Zshiftr,edge_Zshiftr);
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
    if L.P.Na > 0
        angular_modes = reshape(L.P.A_modes-1,1,1,[]);
        sam = sin(angular_modes.*omega_plt);
        cam = cos(angular_modes.*omega_plt);
        R = R+epsilon^2*sum(A.*sam,3);
        Z = Z+epsilon^2*sum(A.*cam,3);
        Rr = Rr+epsilon^2*sum(Ar.*sam,3);
        Zr = Zr+epsilon^2*sum(Ar.*cam,3);
        Rw = Rw+epsilon^2*sum(angular_modes.*A.*cam,3);
        Zw = Zw-epsilon^2*sum(angular_modes.*A.*sam,3);
    end
    if L.P.vertical_shift
        Z = Z+epsilon^2*Zshift;
        Zr = Zr+epsilon^2*Zshiftr;
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

    % Legacy-normalized surface and volume diagnostics.  Pressure is
    % stored internally as Pi=mu0*p/B0^2=epsilon^2*beta, whereas Wk uses
    % the historical beta normalization and is multiplied by R0^3*P0 by
    % dimensional callers.  Gauss weights are used on the solve grid.
    beta_parallel_q = cache.fields.pressure.Pi/epsilon^2;
    beta_perp_q = cache.fields.pressure.Pperp/epsilon^2;
    beta_R_q = cache.fields.pressure.PR/epsilon^2;
    J_q = state.J;
    volume_integral = @(field) ...
        (2*pi)^2*sum(L.w_r.*mean(field.*J_q,2));
    Wkpar = volume_integral(beta_parallel_q);
    Wkperp = volume_integral(beta_perp_q);
    Wk = 0.5*(Wkpar+Wkperp);
    % The legacy Wkrot convention integrates d(beta_parallel)/dR rather
    % than a kinetic energy density; retain that compatibility meaning.
    Wkrot = volume_integral(beta_R_q);
    Wp = volume_integral(state.Bp2);
    edge_current = epsilon*state.edge.psir.*state.edge.goo./state.edge.J;
    Ip = 2*pi*mean(edge_current);
    rBt = state.edge.T;
    Ftt = 2*pi*sum(L.w_r.*mean( ...
        cache.local.Bphi.*state.JoverR,2));
    Ft0 = 2*pi*sum(L.w_r.*mean(state.JoverR./state.R,2));
    Ft = Ftt-rBt*Ft0;

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
    LY.A = A;
    LY.Ap = Ar;
    LY.Zshift = Zshift;
    LY.Zshiftp = Zshiftr;
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
    % Compatibility aliases used by legacy postprocessing scripts.
    LY.BBp2 = LY.Bp2;
    LY.BBt2 = LY.Bphi2;
    LY.BB2 = LY.B2;
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
    LY.dbetapardB0 = state.sigma0/epsilon^2;
    LY.rBt = rBt;
    LY.Wk = Wk;
    LY.Wkpar = Wkpar;
    LY.Wkperp = Wkperp;
    LY.Wkrot = Wkrot;
    LY.Wp = Wp;
    LY.Ip = Ip;
    LY.Ftt = Ftt;
    LY.Ft0 = Ft0;
    LY.Ft = Ft;
    LY.bp = 4*epsilon^2*Wk/Ip^2;
    LY.bppar = 2*epsilon^2*Wkpar/Ip^2;
    LY.bpperp = 2*epsilon^2*Wkperp/Ip^2;
    LY.bprot = (4/3)*epsilon^2*Wkrot/Ip^2;
    LY.li = 2*Wp/Ip^2;
    LY.bpli2 = LY.bppar+LY.bpperp+LY.bprot+LY.li/2;
    LY.gavg = mean(state.R.*cache.local.Bphi,2);
    LY.gauge_error = cache.fields.gauge;
    LY.local_B_residual = max(LY.local_B_residual,local.max_abs_G);
    scaled_J = J(2:end,:)./(epsilon^2*r_plt(2:end));
    interior_J = J(2:end,:);
    LY.min_J = min(interior_J(:));
    LY.min_J_over_eps2r = min(scaled_J(:));
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

    [LY.theta_SFL,LY.RR_sfl,LY.ZZ_sfl, ...
        LY.dthetaSFLdomega,LY.theta_SFL_periodicity_error] = ...
        straight_field_line_geometry( ...
            R,Z,J,sigma,state.a0,epsilon,r_plt,omega_plt);
end

function value = extend_profile(axis_value,quadrature_value,edge_value)
    value = [axis_value;quadrature_value;edge_value];
end

function rows = block_rows(L,profile)
    first = L.profile_starts(profile);
    rows = first:first+L.profile_lengths(profile)-1;
end

function [theta,R_sfl,Z_sfl,dtheta,periodicity_error] = ...
        straight_field_line_geometry( ...
            R,Z,J,sigma,a0,epsilon,r,omega)
%STRAIGHT_FIELD_LINE_GEOMETRY Integrate and invert the exact SFL angle.
% dtheta/domega = a0*J/[epsilon^2*r*R^2*(1-sigma)].
    dtheta = a0*J./(epsilon^2*r.*R.^2.*(1-sigma));
    dtheta(1,:) = 1;
    theta = cumtrapz(omega,dtheta,2);
    theta(1,:) = omega;
    periodicity_error = 2*pi*mean(dtheta(:,1:end-1),2)-2*pi;
    % The discrete gauge constraint uses the solver's periodic angular
    % quadrature, which can differ slightly from cumulative trapezoidal
    % integration.  Normalize so the remapping is exactly 2*pi-periodic.
    theta(2:end,:) = theta(2:end,:).* ...
        (2*pi./theta(2:end,end));

    nr = numel(r);
    nomega = numel(omega);
    theta_uniform = linspace(0,2*pi,nomega);
    R_sfl = ones(nr,nomega);
    Z_sfl = zeros(nr,nomega);
    for ir = 2:nr
        theta_extended = [theta(ir,end)-2*pi,theta(ir,:), ...
            theta(ir,1)+2*pi];
        R_extended = [R(ir,end),R(ir,:),R(ir,1)];
        Z_extended = [Z(ir,end),Z(ir,:),Z(ir,1)];
        [theta_extended,unique_indices] = unique(theta_extended);
        R_extended = R_extended(unique_indices);
        Z_extended = Z_extended(unique_indices);
        R_sfl(ir,:) = interp1( ...
            theta_extended,R_extended,theta_uniform,'spline');
        Z_sfl(ir,:) = interp1( ...
            theta_extended,Z_extended,theta_uniform,'spline');
    end
end
