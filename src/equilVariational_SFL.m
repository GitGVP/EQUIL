function LY = equilVariational_SFL(L,LX,LY,native)
%EQUILVARIATIONAL_SFL Native straight-field-line equilibrium reconstruction.
%   All radial derivatives are evaluated from the variational B-spline
%   representation and all angular derivatives from the Fourier geometry.
%   The returned *_sfl fields include R, Z, their first SFL derivatives,
%   metric and magnetic tensors, anisotropic pressure maps, curvature, and
%   the tensor derivatives consumed by the isotropic and anisotropic Mishka
%   models.  They are
%   therefore suitable for stability interfaces without differentiating
%   sampled R(s,theta) and Z(s,theta) arrays.

    r = LY.r_plt(:);
    omega = LY.omega_plt(:).';
    epsilon = LX.eps_val;
    one_minus_sigma_cyl = LY.one_minus_sigma_cyl(:);
    sigma_cyl_r = LY.sigma_cyl_r(:);
    if isempty(sigma_cyl_r)
        error(['The cylindrical-reference anisotropy needs PrB when ', ...
            'do_SFL=true.']);
    end
    R = LY.RR;
    Z = LY.ZZ;
    Rr = native.Rr;
    Zr = native.Zr;
    Rw = native.Rw;
    Zw = native.Zw;
    J = LY.J;

    [Rrr,Zrr,Rrw,Zrw,Rww,Zww] = ...
        native_second_geometry(L,LX,LY,omega);
    H = Rr.*Zw-Rw.*Zr;
    Jr = Rr.*H+R.*(Rrr.*Zw+Rr.*Zrw-Rrw.*Zr-Rw.*Zrr);
    LY.Jr = Jr;
    Jw = Rw.*H+R.*(Rrw.*Zw+Rr.*Zww-Rww.*Zr-Rw.*Zrw);
    goo = Rw.^2+Zw.^2;
    goor = 2*(Rw.*Rrw+Zw.*Zrw);
    goow = 2*(Rw.*Rww+Zw.*Zww);

    local_state.r = r;
    local_state.R = R;
    local_state.Bp2 = LY.Bp2;
    local_state.T = LY.T;
    fixed = equil_variational_B_constraint(L,LX,local_state,LY.BB);
    pressure = fixed.pressure;
    B = LY.BB;
    Bphi_cyl = fixed.Bphi;
    sigma = fixed.sigma;

    [sigma_r_partial,sigma_R,pressure_PrB,~] = ...
        sigma_partial_derivatives(pressure,B,sigma);
    q = LX.qfun(r); q = q(:);
    qr = LX.qpfun(r); qr = qr(:);
    T = LY.T(:);
    Tr = epsilon^2*LY.t2p(:);
    psir = LY.psir(:);
    radial_flux_derivative = (T+r.*Tr)./q-r.*T.*qr./q.^2;
    psirr = epsilon*( ...
        radial_flux_derivative./one_minus_sigma_cyl ...
        +r.*T.*sigma_cyl_r./(q.*one_minus_sigma_cyl.^2));

    Bp2r = zeros(size(R));
    Bp2w = zeros(size(R));
    interior = 2:numel(r);
    C = epsilon^2;
    Bp2r(interior,:) = C*( ...
        2*psir(interior).*psirr(interior).*goo(interior,:)./J(interior,:).^2 ...
        +psir(interior).^2.*goor(interior,:)./J(interior,:).^2 ...
        -2*psir(interior).^2.*goo(interior,:).*Jr(interior,:)./J(interior,:).^3);
    Bp2w(interior,:) = C*psir(interior).^2.*( ...
        goow(interior,:)./J(interior,:).^2 ...
        -2*goo(interior,:).*Jw(interior,:)./J(interior,:).^3);

    Gr_fixed_B = -(sigma_r_partial+sigma_R.*Rr) ...
        -Tr./(R.*Bphi_cyl)+T.*Rr./(R.^2.*Bphi_cyl) ...
        -T.*Bp2r./(2*R.*Bphi_cyl.^3);
    Gw_fixed_B = -sigma_R.*Rw+T.*Rw./(R.^2.*Bphi_cyl) ...
        -T.*Bp2w./(2*R.*Bphi_cyl.^3);
    Br = -Gr_fixed_B./fixed.GB;
    LY.Br = Br;
    Bw = -Gw_fixed_B./fixed.GB;
    sigmar = sigma_r_partial+sigma_R.*Rr+fixed.sigma_B.*Br;
    sigmaw = sigma_R.*Rw+fixed.sigma_B.*Bw;

    thetaomega = zeros(size(J));
    thetaomega(interior,:) = one_minus_sigma_cyl(interior).*J(interior,:)./( ...
        epsilon^2*r(interior).*R(interior,:).^2.*(1-sigma(interior,:)));
    thetaomega(1,:) = 1;
    thetaomega_r = zeros(size(thetaomega));
    thetaomega_r(interior,:) = thetaomega(interior,:).*( ...
        Jr(interior,:)./J(interior,:)-1./r(interior) ...
        -2*Rr(interior,:)./R(interior,:) ...
        -sigma_cyl_r(interior)./one_minus_sigma_cyl(interior) ...
        +sigmar(interior,:)./(1-sigma(interior,:)));
    thetaomega_w = zeros(size(thetaomega));
    thetaomega_w(interior,:) = thetaomega(interior,:).*( ...
        Jw(interior,:)./J(interior,:)-2*Rw(interior,:)./R(interior,:) ...
        +sigmaw(interior,:)./(1-sigma(interior,:)));
    [thetaomega,thetaomega_r,thetaomega_w] = ...
        regularize_legacy_axis_cubic(L,r,thetaomega,thetaomega_r,thetaomega_w);
    if max(abs(sigma(:))) < 1e-13
        [B,~,Bw] = regularize_legacy_axis_cubic(L,r,B,Br,Bw);
    end

    theta_raw = cumtrapz(omega,thetaomega,2);
    thetar_raw = cumtrapz(omega,thetaomega_r,2);
    period = theta_raw(:,end);
    periodicity_error = period-2*pi;
    scale = ones(size(r));
    scale(interior) = 2*pi./period(interior);
    period_r_over_period = zeros(size(r));
    period_r_over_period(interior) = ...
        thetar_raw(interior,end)./period(interior);
    theta = scale.*theta_raw;
    thetaomega = scale.*thetaomega;
    thetar = scale.*(thetar_raw-theta_raw.*period_r_over_period);
    thetaomega_r = scale.*( ...
        thetaomega_r-(thetaomega./scale).*period_r_over_period);
    thetaomega_w = scale.*thetaomega_w;
    theta(1,:) = omega;
    thetaomega(1,:) = 1;
    thetar(1,:) = 0;
    thetaomega_r(1,:) = 0;
    thetaomega_w(1,:) = 0;

    ratio = thetar./thetaomega;
    ratio_w = (thetaomega_r.*thetaomega-thetar.*thetaomega_w) ...
        ./thetaomega.^2;
    Rs = Rr-ratio.*Rw;
    Zs = Zr-ratio.*Zw;
    Rt = Rw./thetaomega;
    Zt = Zw./thetaomega;
    Rst = (Rrw-ratio_w.*Rw-ratio.*Rww)./thetaomega;
    Zst = (Zrw-ratio_w.*Zw-ratio.*Zww)./thetaomega;
    Rtt = Rww./thetaomega.^2-Rw.*thetaomega_w./thetaomega.^3;
    Ztt = Zww./thetaomega.^2-Zw.*thetaomega_w./thetaomega.^3;

    sigmas = sigmar-ratio.*sigmaw;
    sigmat = sigmaw./thetaomega;
    Bs = Br-ratio.*Bw;
    Bt = Bw./thetaomega;
    g11 = Rs.^2+Zs.^2;
    g12 = Rs.*Rt+Zs.*Zt;
    g22 = Rt.^2+Zt.^2;
    dg12dt = Rst.*Rt+Rs.*Rtt+Zst.*Zt+Zs.*Ztt;
    dg22dt = 2*(Rt.*Rtt+Zt.*Ztt);

    % Ja=J_SFL/s is the regular Jacobian used by VENUS at the axis.
    Ja = epsilon^2*R.^2.*(1-sigma)./one_minus_sigma_cyl;
    dJads = Ja.*(2*Rs./R-sigmas./(1-sigma) ...
        +sigma_cyl_r./one_minus_sigma_cyl);
    dJadt = Ja.*(2*Rt./R-sigmat./(1-sigma));
    dR2ds = 2*R.*Rs;
    dR2dt = 2*R.*Rt;
    B2 = B.^2;
    dB2ds = 2*B.*Bs;
    dB2dt = 2*B.*Bt;

    J_sfl = r.*Ja;
    Btheta = zeros(size(R));
    % psir is normalized with the minor-radius derivative, whereas the SFL
    % geometry uses s=r/a. The extra epsilon converts B^theta to the same
    % R0-normalized coordinate basis as (R_s,R_theta) and gives
    % B^phi/B^theta=q.
    Btheta(interior,:) = epsilon*psir(interior)./J_sfl(interior,:);
    Btheta(1,:) = (epsilon^2*T(1)/( ...
        q(1)*one_minus_sigma_cyl(1)))./Ja(1,:);
    Bphi = T./(R.^2.*(1-sigma));
    Btheta_t = -Btheta.*dJadt./Ja;
    Bphi_t = Bphi.*(-2*Rt./R+sigmat./(1-sigma));

    btheta = Btheta./B;
    bphi = Bphi./B;
    btheta_t = (Btheta_t.*B-Btheta.*Bt)./B.^2;
    bphi_t = (Bphi_t.*B-Bphi.*Bt)./B.^2;
    bR = btheta.*Rt;
    btor = R.*bphi;
    dbRdt = btheta_t.*Rt+btheta.*Rtt;
    dbZdt = btheta_t.*Zt+btheta.*Ztt;
    dbtordt = Rt.*bphi+R.*bphi_t;
    kappaR = btheta.*dbRdt-bphi.*btor;
    kappaZ = btheta.*dbZdt;
    kappator = btheta.*dbtordt+bphi.*bR;
    kappas = zeros(size(R));
    kappat = zeros(size(R));
    kappas(interior,:) = R(interior,:)./J_sfl(interior,:).*( ...
        kappaR(interior,:).*Zt(interior,:) ...
        -kappaZ(interior,:).*Rt(interior,:));
    kappat(interior,:) = R(interior,:)./J_sfl(interior,:).*( ...
        -kappaR(interior,:).*Zs(interior,:) ...
        +kappaZ(interior,:).*Rs(interior,:));
    kappaphi = kappator./R;

    % Analytically kappa=(b dot grad)b is perpendicular to B.  Independent
    % evaluation of the differentiated geometry leaves a small parallel
    % roundoff/interpolation component which grows in the force balance as
    % (P_parallel-P_perp) increases.  Remove that unphysical component
    % algebraically instead of treating it as an equilibrium residual.
    Bmetric2 = g22.*Btheta.^2+R.^2.*Bphi.^2;
    kappa_cov_theta = g12.*kappas+g22.*kappat;
    kappa_cov_phi = R.^2.*kappaphi;
    kappa_dot_B_raw = Btheta.*kappa_cov_theta+Bphi.*kappa_cov_phi;
    curvature_scale = max(1,max(abs([ ...
        Btheta(:).*kappa_cov_theta(:); ...
        Bphi(:).*kappa_cov_phi(:)])));
    LY.sfl_curvature_parallel_error_raw = ...
        max(abs(kappa_dot_B_raw(:)))/curvature_scale;
    kappa_parallel = kappa_dot_B_raw./Bmetric2;
    kappat = kappat-Btheta.*kappa_parallel;
    kappaphi = kappaphi-Bphi.*kappa_parallel;

    Dparallel = pressure.Pr;
    Dperp     = pressure.Pr - B.*pressure_PrB;

    % Use the solved anisotropic force balance to construct the equilibrium
    % current consumed by the stability operator.  Reconstructing j=curl(B)
    % inside the symbolic model makes the marginal internal-kink energy the
    % difference of independently interpolated second-geometry terms.  The
    % representation below is equivalent for a converged equilibrium, has
    % the exact isotropic MISHKA limit, and uses Ampere's law for j^theta to
    % fix the otherwise undetermined field-aligned current.
    Pperp_R = pressure.PR-B.*pressure.PRB;
    Pperp_B = -B.*pressure.PBB;
    Pperp_s = Dperp+Pperp_R.*Rs+Pperp_B.*Bs;
    Pperp_t = Pperp_R.*Rt+Pperp_B.*Bt;
    Bdot_grad_Pperp = Btheta.*Pperp_t;
    BD1 = g12.*Btheta;
    BD2 = g22.*Btheta;
    BD3 = R.^2.*Bphi;
    pressure_delta = pressure.Pi-pressure.Pperp;
    force_s = Pperp_s-BD1.*Bdot_grad_Pperp./B2 ...
        +pressure_delta.*(g11.*kappas+g12.*kappat);
    force_t = Pperp_t-BD2.*Bdot_grad_Pperp./B2 ...
        +pressure_delta.*(g12.*kappas+g22.*kappat);
    force_phi = -BD3.*Bdot_grad_Pperp./B2 ...
        +pressure_delta.*R.^2.*kappaphi;

    % j cross B can represent only the perpendicular force.  The exact EOS
    % and curvature identities make this projection a no-op analytically;
    % numerically it is the least-change removal of the independently
    % reconstructed parallel component.
    force_dot_B_raw = Btheta.*force_t+Bphi.*force_phi;
    force_scale = max(1,max(abs([force_s(:);force_t(:);force_phi(:)])));
    LY.sfl_force_parallel_error_raw = ...
        max(abs(force_dot_B_raw(:)))/force_scale;
    force_parallel = force_dot_B_raw./Bmetric2;
    force_s = force_s-BD1.*force_parallel;
    force_t = force_t-BD2.*force_parallel;
    force_phi = force_phi-BD3.*force_parallel;

    Bphi_s = Tr./(R.^2.*(1-sigma)) ...
        +Bphi.*(-2*Rs./R+sigmas./(1-sigma));
    current_s = zeros(size(R));
    current_theta = zeros(size(R));
    current_phi = zeros(size(R));
    current_s(interior,:) = force_phi(interior,:)./( ...
        J_sfl(interior,:).*Btheta(interior,:));
    current_theta(interior,:) = -( ...
        dR2ds(interior,:).*Bphi(interior,:) ...
        +R(interior,:).^2.*Bphi_s(interior,:))./J_sfl(interior,:);
    current_phi(interior,:) = ( ...
        current_theta(interior,:).*Bphi(interior,:) ...
        -force_s(interior,:)./J_sfl(interior,:))./Btheta(interior,:);
    current_s = extrapolate_axis(current_s,r);
    current_theta = extrapolate_axis(current_theta,r);
    current_phi = extrapolate_axis(current_phi,r);

    force_balance_t_error = force_t+J_sfl.*current_s.*Bphi;
    LY.sfl_force_closure_error = ...
        max(abs(force_balance_t_error(:)))/force_scale;

    source.R = R;
    source.Z = Z;
    source.Rs = Rs;
    source.Rtheta = Rt;
    source.Zs = Zs;
    source.Ztheta = Zt;
    source.B = B;
    source.B2 = B2;
    source.dBds = Bs;
    source.dBdtheta = Bt;
    source.dB2ds = dB2ds;
    source.dB2dtheta = dB2dt;
    source.sigma = sigma;
    source.sigmas = sigmas;
    source.sigmatheta = sigmat;
    source.Pi_parallel = pressure.Pi;
    source.Pi_perp = pressure.Pperp;
    source.Dparallel = Dparallel;
    source.Dperp = Dperp;
    source.kappas = kappas;
    source.kappatheta = kappat;
    source.kappaphi = kappaphi;
    source.current_s = current_s;
    source.current_theta = current_theta;
    source.current_phi = current_phi;
    source.g11 = g11;
    source.g12 = g12;
    source.g22 = g22;
    source.dg12dtheta = dg12dt;
    source.dg22dtheta = dg22dt;
    source.Ja = Ja;
    source.dJads = dJads;
    source.dJadtheta = dJadt;
    source.Bs = zeros(size(R));
    source.Btheta = Btheta;
    source.Bphi = Bphi;
    source.dR2ds = dR2ds;
    source.dR2dtheta = dR2dt;

    theta_uniform = linspace(0,2*pi,numel(omega));
    sfl = remap_fields(source,theta,theta_uniform);
    % Preserve the exact algebraic identities after interpolation.  Direct
    % interpolation of both sides would otherwise introduce a small and
    % entirely artificial metric inconsistency.
    sfl.g11 = sfl.Rs.^2+sfl.Zs.^2;
    sfl.g12 = sfl.Rs.*sfl.Rtheta+sfl.Zs.*sfl.Ztheta;
    sfl.g22 = sfl.Rtheta.^2+sfl.Ztheta.^2;
    sfl.Ja = epsilon^2*sfl.R.^2.*(1-sfl.sigma) ...
        ./one_minus_sigma_cyl;
    % Btheta and Bphi were remapped independently above, whereas Ja is
    % reconstructed nonlinearly from the remapped R and sigma.  Rebuild the
    % contravariant field from its flux profile so that interpolation does
    % not destroy Ja*Btheta=g or Bphi=q*Btheta.
    g_native = epsilon^2*T./(q.*one_minus_sigma_cyl);
    sfl.Btheta = g_native./sfl.Ja;
    sfl.Bphi = q.*sfl.Btheta;
    sfl.dJads = sfl.Ja.*(2*sfl.Rs./sfl.R ...
        -sfl.sigmas./(1-sfl.sigma) ...
        +sigma_cyl_r./one_minus_sigma_cyl);
    sfl.dJadtheta = sfl.Ja.*(2*sfl.Rtheta./sfl.R ...
        -sfl.sigmatheta./(1-sfl.sigma));
    sfl.dR2ds = 2*sfl.R.*sfl.Rs;
    sfl.dR2dtheta = 2*sfl.R.*sfl.Rtheta;
    sfl.B2 = sfl.B.^2;
    sfl.dB2ds = 2*sfl.B.*sfl.dBds;
    sfl.dB2dtheta = 2*sfl.B.*sfl.dBdtheta;
    names = fieldnames(sfl);
    for k = 1:numel(names)
        LY.([names{k},'_sfl']) = sfl.(names{k});
    end
    LY.RR_sfl = sfl.R;
    LY.ZZ_sfl = sfl.Z;
    LY.theta_SFL = theta;
    LY.dthetaSFLdomega = thetaomega;
    LY.d2thetaSFLdrdomega = thetaomega_r;
    LY.dthetaSFLdr = thetar;
    LY.theta_SFL_periodicity_error = periodicity_error;
end

function value = extrapolate_axis(value,r)
% Cubic one-sided continuation of a regular contravariant current field.
    anchors = 2:min(5,numel(r));
    if numel(anchors) < 3
        error('At least three positive-radius samples are needed at the axis.');
    end
    V = [ones(numel(anchors),1),r(anchors), ...
        r(anchors).^2,r(anchors).^3];
    coefficient = V\value(anchors,:);
    value(1,:) = coefficient(1,:);
end

function [Rrr,Zrr,Rrw,Zrw,Rww,Zww] = ...
        native_second_geometry(L,LX,LY,omega)
    nprofiles = 3+L.P.Ns+L.P.Na+double(L.P.vertical_shift);
    second = cell(nprofiles,1);
    for profile = 1:nprofiles
        rows = block_rows(L,profile);
        coefficient = LY.x(rows);
        second{profile} = [ ...
            full(L.B2_axis{profile}*coefficient); ...
            L.B2{profile}*coefficient; ...
            full(L.B2_edge{profile}*coefficient)];
    end

    delta2 = second{2};
    gauge2 = second{3};
    S2 = zeros(numel(LY.r_plt),1,L.P.Ns);
    for is = 1:L.P.Ns
        profile = 3+is;
        S2(:,:,is) = second{profile}+[ ...
            L.Sbc2_axis(is);L.Sbc2{is};L.Sbc2_edge(is)]*LX.Sbc(is);
    end
    A2 = zeros(numel(LY.r_plt),1,L.P.Na);
    for ia = 1:L.P.Na
        profile = 3+L.P.Ns+ia;
        A2(:,:,ia) = second{profile}+[ ...
            L.Abc2_axis(ia);L.Abc2{ia};L.Abc2_edge(ia)]*LX.Abc(ia);
    end
    Zshift2 = zeros(numel(LY.r_plt),1);
    if L.P.vertical_shift
        profile = nprofiles;
        Zshift2 = second{profile}+[ ...
            L.Zbc2_axis;L.Zbc2;L.Zbc2_edge]*LX.Zbc;
    end

    epsilon = LX.eps_val;
    c1 = cos(omega);
    s1 = sin(omega);
    modes = reshape(1:L.P.Ns,1,1,[]);
    cm = cos(modes.*omega);
    sm = sin(modes.*omega);
    Rrr = -epsilon^2*delta2+epsilon^2*sum(S2.*cm,3) ...
        +epsilon^3*gauge2.*c1;
    Zrr = -epsilon^2*sum(S2.*sm,3)+epsilon^3*gauge2.*s1 ...
        +epsilon^2*Zshift2;
    Rrw = -epsilon*s1-epsilon^2*sum(modes.*LY.Sp.*sm,3) ...
        -epsilon^3*LY.Pp.*s1;
    Zrw = epsilon*c1-epsilon^2*sum(modes.*LY.Sp.*cm,3) ...
        +epsilon^3*LY.Pp.*c1;
    Rww = -epsilon*LY.r_plt.*c1 ...
        -epsilon^2*sum(modes.^2.*LY.S.*cm,3)-epsilon^3*LY.P.*c1;
    Zww = -epsilon*LY.r_plt.*s1 ...
        +epsilon^2*sum(modes.^2.*LY.S.*sm,3)-epsilon^3*LY.P.*s1;
    if L.P.Na > 0
        angular_modes = reshape(L.P.A_modes-1,1,1,[]);
        sam = sin(angular_modes.*omega);
        cam = cos(angular_modes.*omega);
        Rrr = Rrr+epsilon^2*sum(A2.*sam,3);
        Zrr = Zrr+epsilon^2*sum(A2.*cam,3);
        Rrw = Rrw+epsilon^2*sum(angular_modes.*LY.Ap.*cam,3);
        Zrw = Zrw-epsilon^2*sum(angular_modes.*LY.Ap.*sam,3);
        Rww = Rww-epsilon^2*sum(angular_modes.^2.*LY.A.*sam,3);
        Zww = Zww-epsilon^2*sum(angular_modes.^2.*LY.A.*cam,3);
    end
end

function [sigma_r,sigma_R,PrB,PRB] = ...
        sigma_partial_derivatives(pressure,B,sigma)
    if isfield(pressure,'PrB') && isfield(pressure,'PRB')
        PrB = pressure.PrB;
        PRB = pressure.PRB;
        sigma_r = PrB./B;
        sigma_R = PRB./B;
    elseif max(abs(sigma(:))) < 1e-13
        PrB = zeros(size(B));
        PRB = zeros(size(B));
        sigma_r = zeros(size(B));
        sigma_R = zeros(size(B));
    else
        error(['An anisotropic pressure_map used with do_SFL=true must ', ...
            'provide PrB and PRB.']);
    end
end

function mapped = remap_fields(source,theta,theta_uniform)
    names = fieldnames(source);
    nr = size(theta,1);
    nw = size(theta,2);
    values = zeros(nw,numel(names));
    for k = 1:numel(names)
        mapped.(names{k}) = zeros(nr,nw);
    end
    for ir = 1:nr
        for k = 1:numel(names)
            field = source.(names{k});
            values(:,k) = field(ir,:).';
        end
        theta_extended = [theta(ir,end)-2*pi,theta(ir,:), ...
            theta(ir,1)+2*pi];
        values_extended = [values(end,:);values;values(1,:)];
        [theta_extended,unique_indices] = unique(theta_extended);
        values_extended = values_extended(unique_indices,:);
        interpolated = interp1(theta_extended,values_extended, ...
            theta_uniform,'spline');
        for k = 1:numel(names)
            mapped.(names{k})(ir,:) = interpolated(:,k).';
        end
    end
end

function [value,radial,angular] = ...
        regularize_legacy_axis_cubic(L,r,value,radial,angular)
% The legacy-compatible space does not impose polar Taylor regularity.
% High-degree profiles can consequently have large, cancelling second
% derivatives in the first span even when R and Z are well resolved.  A
% cubic continuation of theta_omega from the resolved outer part of that
% span supplies the regular axis limit without differentiating sampled SFL
% geometry.  Axis-regular spaces use the exact expressions unchanged.
    if ~strcmp(L.P.radial_discretization,'legacy')
        return
    end
    first_span = find(r < 1/L.P.m);
    positive = first_span(r(first_span) > 0);
    if numel(positive) < 3
        error('Legacy SFL axis regularization requires three first-span points.');
    end
    anchors = [1;positive(end-2:end)];
    radius = r(positive(end));
    x = r(anchors)/radius;
    V = [ones(4,1),x,x.^2,x.^3];
    coefficient = V\value(anchors,:);
    angular_coefficient = V\angular(anchors,:);
    rows = first_span;
    xr = r(rows)/radius;
    E = [ones(numel(rows),1),xr,xr.^2,xr.^3];
    E1 = [zeros(numel(rows),1),ones(numel(rows),1), ...
        2*xr,3*xr.^2]/radius;
    value(rows,:) = E*coefficient;
    radial(rows,:) = E1*coefficient;
    angular(rows,:) = E*angular_coefficient;
    if any(value(rows,:) <= 0)
        error('The regularized legacy theta_omega is not monotone.');
    end
end

function rows = block_rows(L,profile)
    first = L.profile_starts(profile);
    rows = first:first+L.profile_lengths(profile)-1;
end
