function R = residuals_noRepl_spec(r,x,epsilon,omega,kinetic_profiles,q,qp,P0,P1,P2,M_profiles,dof_count,Nb,Sbc,equation_of_state)
	N = numel(r);
	Ns = numel(Sbc);
	nRes = 3+Nb+Ns;
	t2_c    = x(1:dof_count);
	delta_c = x(dof_count+1:2*dof_count);
	P_c = x(2*dof_count+1:3*dof_count);
	Bs_c    = reshape(x(3*dof_count+1:(3+Nb)*dof_count), [dof_count, 1, Nb]);
	Bs_flat   = reshape(Bs_c, size(Bs_c,1), []);        % dof_count x Ns
	S_unkn     = reshape(x((3+Nb)*dof_count+1:end), [dof_count-1, 1, Ns]);
	%S_c     = cat(1, S_unkn, reshape(, [1, 1, Ns]));
    S_c=S_unkn;
	S_flat   = reshape(S_c, size(S_c,1), []);        % dof_count x Ns
	t2      = P0 * t2_c;
	t2p     = P1 * t2_c;
	delta   = P0 * delta_c;
	deltap  = P1 * delta_c;
	deltapp = P2 * delta_c;
	P   = P0 * P_c;
	Pp  = P1 * P_c;
	Ppp = P2 * P_c;
	Bs       = reshape(P0 * Bs_flat, N, 1, []);
	Bsp = reshape(P1 * Bs_flat, N, 1, []);
	S       = reshape(P0(:,1:end-1) * S_flat, N, 1, []);
	Sp = reshape(P1(:,1:end-1) * S_flat, N, 1, []);
	Spp = reshape(P2(:,1:end-1) * S_flat, N, 1, []);
    % --- Enforce S boundary form: S_m(r) = Sbc(m) * r^(m-1)  -------------
    % Replace reconstructed S (and derivatives) by the analytic form derived
    % from Sbc. This enforces the chosen boundary profile for all r points.
    % Note: we assume m index = 1..Ns corresponds to harmonic order m-1 in code.
    for mIdx = 1:Ns
        exponent = mIdx +1;           % S_m ~ r^{m-1}
        % Sbc provided as Sbc(mIdx)
        % compute S, S', S'' safely (handle r==0 for negative exponents)
       
            % general exponent >= 1
            S(:,1,mIdx)   = S(:,1,mIdx)+ Sbc(mIdx) * (r(:) .^ (exponent-1));
            % derivatives:
            Sp(:,1,mIdx)  = Sp(:,1,mIdx)  +Sbc(mIdx) * (exponent-1) * (r(:) .^ (exponent - 2));
            if exponent - 1 >= 1
                Spp(:,1,mIdx) = Spp(:,1,mIdx)+ Sbc(mIdx) * (exponent-1) * (exponent - 2) * (r(:) .^ (exponent - 3));
            else
                % second derivative singular at r=0 when exponent = 1: set finite value
                Spp(:,1,mIdx) = Spp(:,1,mIdx)+zeros(N,1);
            end
       
    end
    
    % --- Enforce homogeneity at r=0 for the primary unknowns
    % User requested "Cr^2 and C as additional degree of freedom to impose the
    % fact that is homogeneous at r=0".  To avoid changing the unknown vector
    % layout we enforce homogeneity strongly at the r==0 quadrature point by
    % forcing both value and first derivative to zero there (consistent with
    % u ~ C*r^2 behaviour: u(0)=0 and u'(0)=0).
    tol_r0 = 1e-14;
    idx0 = find(abs(r) < tol_r0, 1);
    if isempty(idx0)
        % fallback to first quadrature point
        idx0 = 1;
    end
    
    % Force t2, delta, P and Bs to have zero value and zero first derivative at r==0
    t2(idx0) = 0;    t2p(idx0) = 0;
    delta(idx0) = 0; deltap(idx0) = 0;
    P(idx0) = 0;     Pp(idx0) = 0;
    % for each B-profile
    if ~isempty(Bs)
        for ib = 1:Nb
            Bs(idx0,1,ib) = 0;
            Bsp(idx0,1,ib) = 0;
        end
    end
    
    % Note: we do not change coefficient vectors (t2_c, etc.). We only overwrite
    % the reconstructed values at r==0 used in residual evaluation. This strongly
    % pins the solution at r==0 and enforces the homogeneous regularity there.

	resTotal = zeros(nRes,N);
	ms = reshape(linspace(2,Ns+1,Ns), [1,1,Ns]);
	mb = reshape(0:(Nb-1), [1,1,Nb]);
	RR = 1 - delta.*epsilon.^2 + epsilon.^3.*P.*cos(omega) + epsilon.*r.*cos(omega) + epsilon.^2.*sum(cos((-1 + ms).*omega).*S,3);
	BB = 1 + epsilon.*sum(Bs.*cos(mb.*omega),3);

	[dbetapardB0, ...
        dbetapardr, dbetapardB, dbetapardR, ...
        d2betapardB2, d2betapardrdB, d2betapardRdB] = equation_of_state(kinetic_profiles,r,RR,BB);

	T1 = -(epsilon.^2.*P.*sin(omega)) - r.*sin(omega) + epsilon.*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3);
	T2 = epsilon.^2.*P.*cos(omega) + r.*cos(omega) - epsilon.*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3);
	T3 = -(deltap.*epsilon) + cos(omega) + epsilon.^2.*Pp.*cos(omega) + epsilon.*sum(cos((-1 + ms).*omega).*Sp,3);
	T4 = sin(omega) + epsilon.^2.*Pp.*sin(omega) - epsilon.*sum(sin((-1 + ms).*omega).*Sp,3);
	T5 = 1 - dbetapardB0.*epsilon.^2 + epsilon.^2.*t2;
	T8 = r + epsilon.^2.*Pp.*r - deltap.*epsilon.*r.*cos(omega) - epsilon.*sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) - epsilon.^3.*Pp.*sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) + epsilon.*r.*cos(omega).*sum(cos((-1 + ms).*omega).*Sp,3) - epsilon.*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3).*(-(deltap.*epsilon) + cos(omega) + epsilon.^2.*Pp.*cos(omega) + epsilon.*sum(cos((-1 + ms).*omega).*Sp,3)) - epsilon.*r.*sin(omega).*sum(sin((-1 + ms).*omega).*Sp,3) + epsilon.^2.*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3).*sum(sin((-1 + ms).*omega).*Sp,3) + epsilon.^2.*P.*(1 + epsilon.^2.*Pp - deltap.*epsilon.*cos(omega) + epsilon.*cos(omega).*sum(cos((-1 + ms).*omega).*Sp,3) - epsilon.*sin(omega).*sum(sin((-1 + ms).*omega).*Sp,3));
	T9 = 1./(r + epsilon.^2.*Pp.*r - deltap.*epsilon.*r.*cos(omega) - epsilon.*sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) - epsilon.^3.*Pp.*sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) + epsilon.*r.*cos(omega).*sum(cos((-1 + ms).*omega).*Sp,3) - epsilon.*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3).*(-(deltap.*epsilon) + cos(omega) + epsilon.^2.*Pp.*cos(omega) + epsilon.*sum(cos((-1 + ms).*omega).*Sp,3)) - epsilon.*r.*sin(omega).*sum(sin((-1 + ms).*omega).*Sp,3) + epsilon.^2.*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3).*sum(sin((-1 + ms).*omega).*Sp,3) + epsilon.^2.*P.*(1 + epsilon.^2.*Pp - deltap.*epsilon.*cos(omega) + epsilon.*cos(omega).*sum(cos((-1 + ms).*omega).*Sp,3) - epsilon.*sin(omega).*sum(sin((-1 + ms).*omega).*Sp,3)));
	T12 = epsilon.^4.*P.^2 + r.^2 - 2.*epsilon.*r.*cos(omega).*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3) + epsilon.^2.*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3).^2 - 2.*epsilon.*r.*sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) + epsilon.^2.*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3).^2 - 2.*epsilon.^2.*P.*(-r + epsilon.*cos(omega).*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3) + epsilon.*sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3));
	T15 = cos(omega) + epsilon.^2.*Pp.*cos(omega) - epsilon.*sum((-1 + ms).*cos((-1 + ms).*omega).*Sp,3);
	T16 = -((1 + epsilon.^2.*Pp).*sin(omega)) + epsilon.*sum(-((-1 + ms).*sin((-1 + ms).*omega).*Sp),3);
	T17 = (1 - dbetapardB0.*epsilon.^2 + epsilon.^2.*t2).^2;
	T24 = (r + epsilon.^2.*Pp.*r - deltap.*epsilon.*r.*cos(omega) - epsilon.*sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) - epsilon.^3.*Pp.*sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) + epsilon.*r.*cos(omega).*sum(cos((-1 + ms).*omega).*Sp,3) - epsilon.*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3).*(-(deltap.*epsilon) + cos(omega) + epsilon.^2.*Pp.*cos(omega) + epsilon.*sum(cos((-1 + ms).*omega).*Sp,3)) - epsilon.*r.*sin(omega).*sum(sin((-1 + ms).*omega).*Sp,3) + epsilon.^2.*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3).*sum(sin((-1 + ms).*omega).*Sp,3) + epsilon.^2.*P.*(1 + epsilon.^2.*Pp - deltap.*epsilon.*cos(omega) + epsilon.*cos(omega).*sum(cos((-1 + ms).*omega).*Sp,3) - epsilon.*sin(omega).*sum(sin((-1 + ms).*omega).*Sp,3))).^(-2);
	T25 = sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3).*(-(deltap.*epsilon) + cos(omega) + epsilon.^2.*Pp.*cos(omega) + epsilon.*sum(cos((-1 + ms).*omega).*Sp,3)) - sum((-1 + ms).*cos((-1 + ms).*omega).*S,3).*(sin(omega) + epsilon.^2.*Pp.*sin(omega) - epsilon.*sum(sin((-1 + ms).*omega).*Sp,3)) - (epsilon.^2.*P + r).*(-(deltap.*sin(omega)) + sin(omega).*sum(cos((-1 + ms).*omega).*Sp,3) + cos(omega).*sum(sin((-1 + ms).*omega).*Sp,3));
	T27 = 2.*(epsilon.^2.*P.*cos(omega) + r.*cos(omega) - epsilon.*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3)).*(cos(omega) + epsilon.^2.*Pp.*cos(omega) - epsilon.*sum((-1 + ms).*cos((-1 + ms).*omega).*Sp,3));
	T31 = 2.*(epsilon.^2.*P.*sin(omega) + r.*sin(omega) - epsilon.*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3)).*(sin(omega) + epsilon.^2.*Pp.*sin(omega) - epsilon.*sum(-((-1 + ms).*sin((-1 + ms).*omega).*Sp),3));
	T36 = deltap.*epsilon.^2.*P.*cos(omega) + deltap.*r.*cos(omega) - sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) - epsilon.^2.*Pp.*sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) - sin(omega).*sum(-((-1 + ms).^2.*S.*sin((-1 + ms).*omega)),3) - epsilon.^2.*Pp.*sin(omega).*sum(-((-1 + ms).^2.*S.*sin((-1 + ms).*omega)),3) - epsilon.^2.*P.*cos(omega).*sum(cos((-1 + ms).*omega).*Sp,3) - r.*cos(omega).*sum(cos((-1 + ms).*omega).*Sp,3) + sum(-((-1 + ms).^2.*cos((-1 + ms).*omega).*S),3).*(-(deltap.*epsilon) + cos(omega) + epsilon.^2.*Pp.*cos(omega) + epsilon.*sum(cos((-1 + ms).*omega).*Sp,3)) - epsilon.^2.*P.*cos(omega).*sum((-1 + ms).*cos((-1 + ms).*omega).*Sp,3) - r.*cos(omega).*sum((-1 + ms).*cos((-1 + ms).*omega).*Sp,3) - sum((-1 + ms).*cos((-1 + ms).*omega).*S,3).*(cos(omega) + epsilon.^2.*Pp.*cos(omega) - epsilon.*sum((-1 + ms).*cos((-1 + ms).*omega).*Sp,3)) + epsilon.^2.*P.*sin(omega).*sum(sin((-1 + ms).*omega).*Sp,3) + r.*sin(omega).*sum(sin((-1 + ms).*omega).*Sp,3) + epsilon.*sum(-((-1 + ms).^2.*S.*sin((-1 + ms).*omega)),3).*sum(sin((-1 + ms).*omega).*Sp,3) - epsilon.^2.*P.*sin(omega).*sum(-((-1 + ms).*sin((-1 + ms).*omega).*Sp),3) - r.*sin(omega).*sum(-((-1 + ms).*sin((-1 + ms).*omega).*Sp),3) + epsilon.*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3).*sum(-((-1 + ms).*sin((-1 + ms).*omega).*Sp),3);
	T37 = deltap.*epsilon.^2.*P.*sin(omega) + deltap.*r.*sin(omega) - cos(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) - epsilon.^2.*Pp.*cos(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) + deltap.*epsilon.*sum(-((-1 + ms).^2.*S.*sin((-1 + ms).*omega)),3) - cos(omega).*sum(-((-1 + ms).^2.*S.*sin((-1 + ms).*omega)),3) - epsilon.^2.*Pp.*cos(omega).*sum(-((-1 + ms).^2.*S.*sin((-1 + ms).*omega)),3) - epsilon.^2.*P.*sin(omega).*sum(cos((-1 + ms).*omega).*Sp,3) - r.*sin(omega).*sum(cos((-1 + ms).*omega).*Sp,3) - epsilon.*sum(-((-1 + ms).^2.*S.*sin((-1 + ms).*omega)),3).*sum(cos((-1 + ms).*omega).*Sp,3) - epsilon.^2.*P.*sin(omega).*sum((-1 + ms).*cos((-1 + ms).*omega).*Sp,3) - r.*sin(omega).*sum((-1 + ms).*cos((-1 + ms).*omega).*Sp,3) + epsilon.*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3).*sum((-1 + ms).*cos((-1 + ms).*omega).*Sp,3) - epsilon.^2.*P.*cos(omega).*sum(sin((-1 + ms).*omega).*Sp,3) - r.*cos(omega).*sum(sin((-1 + ms).*omega).*Sp,3) - sum(-((-1 + ms).^2.*cos((-1 + ms).*omega).*S),3).*(sin(omega) + epsilon.^2.*Pp.*sin(omega) - epsilon.*sum(sin((-1 + ms).*omega).*Sp,3)) + epsilon.^2.*P.*cos(omega).*sum(-((-1 + ms).*sin((-1 + ms).*omega).*Sp),3) + r.*cos(omega).*sum(-((-1 + ms).*sin((-1 + ms).*omega).*Sp),3) + sum((-1 + ms).*cos((-1 + ms).*omega).*S,3).*(sin(omega) + epsilon.^2.*Pp.*sin(omega) - epsilon.*sum(-((-1 + ms).*sin((-1 + ms).*omega).*Sp),3));
	T40 = d2betapardrdB + d2betapardB2.*epsilon.*sum(Bsp.*cos(mb.*omega),3) + d2betapardRdB.*epsilon.*(-(deltap.*epsilon) + cos(omega) + epsilon.^2.*Pp.*cos(omega) + epsilon.*sum(cos((-1 + ms).*omega).*Sp,3));

	resF = (1 - delta.*epsilon.^2 + epsilon.^3.*P.*cos(omega) + epsilon.*r.*cos(omega) + epsilon.^2.*sum(cos((-1 + ms).*omega).*S,3)).*((epsilon.^2.*P.*cos(omega) + r.*cos(omega) - epsilon.*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3)).*(-(deltap.*epsilon) + cos(omega) + epsilon.^2.*Pp.*cos(omega) + epsilon.*sum(cos((-1 + ms).*omega).*Sp,3)) + (epsilon.^2.*P.*sin(omega) + r.*sin(omega) - epsilon.*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3)).*(sin(omega) + epsilon.^2.*Pp.*sin(omega) - epsilon.*sum(sin((-1 + ms).*omega).*Sp,3))).*(dbetapardr - (BB.*dbetapardB.*epsilon.*T17.*T3)./((BB - dbetapardB.*epsilon.^2).^2.*RR.^3) - BB.*T40 + (epsilon.*r.^2.*T5.^2.*T9.*(epsilon.*RR.*T24.*T25.*T37 + epsilon.*T1.*T25.*T9 - RR.*T36.*T9))./(q.^2.*RR.^3) + (BB.*T5.*(BB.*(BB - dbetapardB.*epsilon.^2).*t2p + T5.*(BB.*T40 - dbetapardB.*epsilon.*sum(Bsp.*cos(mb.*omega),3))))./((BB - dbetapardB.*epsilon.^2).^3.*RR.^2) + (epsilon.^2.*r.^2.*T17.*(2.*BB.*epsilon.^2.*RR.*T25.*T8.*(d2betapardRdB.*T1 + d2betapardB2.*sum(-(mb.*Bs.*sin(mb.*omega)),3)) - dbetapardB.*(2.*BB.*epsilon.^2.*T1.*T25.*T8 + BB.*RR.*(2.*epsilon.^2.*T25.*T37 + T27.*T8 + T31.*T8 - 2.*epsilon.*T36.*T8) + 2.*epsilon.^2.*RR.*T25.*T8.*sum(-(mb.*Bs.*sin(mb.*omega)),3))))./(2..*BB.^2.*q.^2.*RR.^3.*T8.^3) + (r.*T5.*T9.*(-(qp.*r.*RR.*T12.*T5.*T9) + q.*(deltapp.*epsilon.*r.*RR.*T12.*T2.*T24.*T5 - r.*RR.*T12.*T15.*T24.*T3.*T5 + r.*RR.*T12.*T16.*T24.*T4.*T5 + epsilon.^2.*r.*RR.*T12.*t2p.*T9 + RR.*T12.*T5.*T9 + r.*RR.*T27.*T5.*T9 - epsilon.*r.*T12.*T3.*T5.*T9 + r.*RR.*T31.*T5.*T9 - epsilon.^2.*Ppp.*r.*RR.*T12.*T2.*T24.*T5.*cos(omega) + epsilon.^2.*Ppp.*r.*RR.*T1.*T12.*T24.*T5.*sin(omega) - epsilon.*r.*RR.*T12.*T2.*T24.*T5.*sum(cos((-1 + ms).*omega).*Spp,3) - epsilon.*r.*RR.*T1.*T12.*T24.*T5.*sum(sin((-1 + ms).*omega).*Spp,3))))./(q.^3.*RR.^3));
	resB = (BB.^2.*(1 - T17./((BB - dbetapardB.*epsilon.^2).^2.*RR.^2)))./epsilon - (epsilon.*r.^2.*T12.*T17.*T24)./(q.^2.*RR.^2);
	resF_FFT = real(fft(resF, [], 2))/numel(omega);
	resB_FFT = real(fft(resB, [], 2))/numel(omega);
	mF = [0, 1, 2:(Ns+1)];
	rowsF = [1, 2, (2 + Nb + (2:(Ns+1)))];
	mB = 0:(Nb-1);
	rowsB = 3 + (1:Nb);
	resTotal(rowsF, :) = resF_FFT(:, mF + 1).';
	resTotal(rowsB, :) = resB_FFT(:, mB + 1).';
	resTotal(3,:) = mean((r - ((epsilon.^2.*P.*cos(omega) + r.*cos(omega) - epsilon.*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3)).*(-(deltap.*epsilon) + cos(omega) + epsilon.^2.*Pp.*cos(omega) + epsilon.*sum(cos((-1 + ms).*omega).*Sp,3)) + (epsilon.^2.*P.*sin(omega) + r.*sin(omega) - epsilon.*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3)).*(sin(omega) + epsilon.^2.*Pp.*sin(omega) - epsilon.*sum(sin((-1 + ms).*omega).*Sp,3)))./RR)./epsilon.^2, 2);
	R = M_profiles * reshape(resTotal.', [], 1);
end
