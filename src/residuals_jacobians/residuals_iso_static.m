function R = residuals_iso_static(r,x,epsilon,omega,kinetic_profiles,q,qp,P0,P1,P2,M_profiles,dof_count,Nb,Sbc,equation_of_state)
	N = numel(r);
	Ns = numel(Sbc);
	nRes = 3+Nb+Ns;
	[t2, t2p, delta, deltap, deltapp, P, Pp, Ppp, Bs, ~, S, Sp, Spp] = ...
          equil_unpack_state(x, dof_count, Nb, Sbc, P0, P1, P2);
	resTotal = zeros(nRes,N);
	ms = reshape(linspace(2,Ns+1,Ns), [1,1,Ns]);
	mb = reshape(0:(Nb-1), [1,1,Nb]);
	RR = 1 - delta.*epsilon.^2 + epsilon.^3.*P.*cos(omega) + epsilon.*r.*cos(omega) + epsilon.^2.*sum(cos((-1 + ms).*omega).*S,3);
    dbetapardr = kinetic_profiles.betap(r);
    RBphi = 1 + epsilon.^2.*t2;
    temp2 = epsilon.^2.*P.*sin(omega) + r.*sin(omega) - epsilon.*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3);
    T1 = -(epsilon.^2.*P.*sin(omega)) - r.*sin(omega) + epsilon.*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3);
	T2 = epsilon.^2.*P.*cos(omega) + r.*cos(omega) - epsilon.*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3);
	T3 = -(deltap.*epsilon) + cos(omega) + epsilon.^2.*Pp.*cos(omega) + epsilon.*sum(cos((-1 + ms).*omega).*Sp,3);
	temp1 = (r + epsilon.^2.*Pp.*r - deltap.*epsilon.*r.*cos(omega) - epsilon.*sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) - epsilon.^3.*Pp.*sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) + epsilon.*r.*cos(omega).*sum(cos((-1 + ms).*omega).*Sp,3) - epsilon.*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3).*(T3) - epsilon.*r.*sin(omega).*sum(sin((-1 + ms).*omega).*Sp,3) + epsilon.^2.*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3).*sum(sin((-1 + ms).*omega).*Sp,3) + epsilon.^2.*P.*(1 + epsilon.^2.*Pp - deltap.*epsilon.*cos(omega) + epsilon.*cos(omega).*sum(cos((-1 + ms).*omega).*Sp,3) - epsilon.*sin(omega).*sum(sin((-1 + ms).*omega).*Sp,3)));
    T4 = sin(omega) + epsilon.^2.*Pp.*sin(omega) - epsilon.*sum(sin((-1 + ms).*omega).*Sp,3);
	T8 = 1./temp1;
	T11 = epsilon.^4.*P.^2 + r.^2 - 2.*epsilon.*r.*cos(omega).*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3) + epsilon.^2.*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3).^2 - 2.*epsilon.*r.*sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) + epsilon.^2.*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3).^2 - 2.*epsilon.^2.*P.*(-r + epsilon.*cos(omega).*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3) + epsilon.*sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3));
	T14 = cos(omega) + epsilon.^2.*Pp.*cos(omega) - epsilon.*sum((-1 + ms).*cos((-1 + ms).*omega).*Sp,3);
	T15 = -((1 + epsilon.^2.*Pp).*sin(omega)) + epsilon.*sum(-((-1 + ms).*sin((-1 + ms).*omega).*Sp),3);
	T16 = temp1.^(-2);
	T18 = -(epsilon.^2.*P.*cos(omega)) - r.*cos(omega) + epsilon.*sum(-((-1 + ms).^2.*cos((-1 + ms).*omega).*S),3);
	T20 = -(epsilon.^2.*P.*sin(omega)) - r.*sin(omega) - epsilon.*sum(-((-1 + ms).^2.*S.*sin((-1 + ms).*omega)),3);
	T23 = sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3).*(T3) - sum((-1 + ms).*cos((-1 + ms).*omega).*S,3).*(T4) - (epsilon.^2.*P + r).*(-(deltap.*sin(omega)) + sin(omega).*sum(cos((-1 + ms).*omega).*Sp,3) + cos(omega).*sum(sin((-1 + ms).*omega).*Sp,3));

	resF = (RR).*((T2).*(T3) + (temp2).*(T4)).*(dbetapardr + ((RBphi).*t2p)./RR.^2 - ((r + epsilon.^2.*r.*t2).^2.*T8.*(-(epsilon.^2.*T1.*T23.*T8) + RR.*(epsilon.*T1.*T14.*T16.*T23 - epsilon.*T15.*T16.*T2.*T23 - epsilon.*T16.*T20.*T23.*T3 + epsilon.*T16.*T18.*T23.*T4 + T1.*T15.*T8 + T14.*T2.*T8 + T18.*T3.*T8 + T20.*T4.*T8)))./(q.^2.*RR.^3) + (r.*(RBphi).*T8.*(q.*RR.*T11.*(RBphi).*T8 - qp.*r.*RR.*T11.*(RBphi).*T8 + 2.*q.*r.*RR.*(RBphi).*(T1.*T15 + T14.*T2).*T8 + epsilon.^2.*q.*r.*RR.*T11.*t2p.*T8 - epsilon.*q.*r.*T11.*(RBphi).*T3.*T8 - q.*r.*RR.*T11.*T16.*(RBphi).*(-(deltapp.*epsilon.*T2) + T14.*T3 - T15.*T4 + epsilon.^2.*Ppp.*T2.*cos(omega) - epsilon.^2.*Ppp.*T1.*sin(omega) + epsilon.*T2.*sum(cos((-1 + ms).*omega).*Spp,3) + epsilon.*T1.*sum(sin((-1 + ms).*omega).*Spp,3))))./(q.^3.*RR.^3));
	% dummy residual for Bs
	resB = epsilon.*sum(Bs.*cos(mb.*omega),3);
    resF_FFT = real(fft(resF, [], 2))/numel(omega);
	resB_FFT = real(fft(resB, [], 2))/numel(omega);
	mF = [0, 1, 2:(Ns+1)];
	rowsF = [1, 2, (2 + Nb + (2:(Ns+1)))];
	mB = 0:(Nb-1);
	rowsB = 3 + (1:Nb);
	resTotal(rowsF, :) = resF_FFT(:, mF + 1).';
	resTotal(rowsB, :) = resB_FFT(:, mB + 1).';
	resTotal(3,:) = mean((r - ((T2).*(T3) + (temp2).*(T4))./RR)./epsilon.^2, 2);
	R = M_profiles * reshape(resTotal.', [], 1);
end
