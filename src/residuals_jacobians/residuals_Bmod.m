function R = residuals_Bmod(r,x,epsilon,omega,kinetic_profiles,q,qp,P0,P1,P2,M_profiles,dof_count,Nb,Sbc,equation_of_state)
	N = numel(r);
	Ns = numel(Sbc);
	nRes = 3+Nb+Ns;
	t2_c    = x(1:dof_count);
	delta_c = x(dof_count+1:2*dof_count);
	P_c = x(2*dof_count+1:3*dof_count);
	Bs_c    = reshape(x(3*dof_count+1:(3+Nb)*dof_count), [dof_count, 1, Nb]);
	Bs_flat   = reshape(Bs_c, size(Bs_c,1), []);        % dof_count x Ns
	S_unkn     = reshape(x((3+Nb)*dof_count+1:end), [dof_count-1, 1, Ns]);
	S_c     = cat(1, S_unkn, reshape(Sbc, [1, 1, Ns]));
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
	S       = reshape(P0 * S_flat, N, 1, []);
	Sp = reshape(P1 * S_flat, N, 1, []);
	Spp = reshape(P2 * S_flat, N, 1, []);
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
	T9 = 1./(r + epsilon.^2.*Pp.*r - deltap.*epsilon.*r.*cos(omega) - epsilon.*sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) - epsilon.^3.*Pp.*sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) + epsilon.*r.*cos(omega).*sum(cos((-1 + ms).*omega).*Sp,3) - epsilon.*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3).*(-(deltap.*epsilon) + cos(omega) + epsilon.^2.*Pp.*cos(omega) + epsilon.*sum(cos((-1 + ms).*omega).*Sp,3)) - epsilon.*r.*sin(omega).*sum(sin((-1 + ms).*omega).*Sp,3) + epsilon.^2.*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3).*sum(sin((-1 + ms).*omega).*Sp,3) + epsilon.^2.*P.*(1 + epsilon.^2.*Pp - deltap.*epsilon.*cos(omega) + epsilon.*cos(omega).*sum(cos((-1 + ms).*omega).*Sp,3) - epsilon.*sin(omega).*sum(sin((-1 + ms).*omega).*Sp,3)));
	T12 = epsilon.^4.*P.^2 + r.^2 - 2.*epsilon.*r.*cos(omega).*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3) + epsilon.^2.*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3).^2 - 2.*epsilon.*r.*sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) + epsilon.^2.*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3).^2 - 2.*epsilon.^2.*P.*(-r + epsilon.*cos(omega).*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3) + epsilon.*sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3));
	T13 = (1 - dbetapardB0.*epsilon.^2 + epsilon.^2.*t2).^2;
	T16 = cos(omega) + epsilon.^2.*Pp.*cos(omega) - epsilon.*sum((-1 + ms).*cos((-1 + ms).*omega).*Sp,3);
	T17 = -((1 + epsilon.^2.*Pp).*sin(omega)) + epsilon.*sum(-((-1 + ms).*sin((-1 + ms).*omega).*Sp),3);
	T18 = (r + epsilon.^2.*Pp.*r - deltap.*epsilon.*r.*cos(omega) - epsilon.*sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) - epsilon.^3.*Pp.*sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) + epsilon.*r.*cos(omega).*sum(cos((-1 + ms).*omega).*Sp,3) - epsilon.*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3).*(-(deltap.*epsilon) + cos(omega) + epsilon.^2.*Pp.*cos(omega) + epsilon.*sum(cos((-1 + ms).*omega).*Sp,3)) - epsilon.*r.*sin(omega).*sum(sin((-1 + ms).*omega).*Sp,3) + epsilon.^2.*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3).*sum(sin((-1 + ms).*omega).*Sp,3) + epsilon.^2.*P.*(1 + epsilon.^2.*Pp - deltap.*epsilon.*cos(omega) + epsilon.*cos(omega).*sum(cos((-1 + ms).*omega).*Sp,3) - epsilon.*sin(omega).*sum(sin((-1 + ms).*omega).*Sp,3))).^(-2);
	T22 = sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3).*(-(deltap.*epsilon) + cos(omega) + epsilon.^2.*Pp.*cos(omega) + epsilon.*sum(cos((-1 + ms).*omega).*Sp,3)) - sum((-1 + ms).*cos((-1 + ms).*omega).*S,3).*(sin(omega) + epsilon.^2.*Pp.*sin(omega) - epsilon.*sum(sin((-1 + ms).*omega).*Sp,3)) - (epsilon.^2.*P + r).*(-(deltap.*sin(omega)) + sin(omega).*sum(cos((-1 + ms).*omega).*Sp,3) + cos(omega).*sum(sin((-1 + ms).*omega).*Sp,3));
	T27 = d2betapardB2.*sum(-(mb.*Bs.*sin(mb.*omega)),3) - d2betapardRdB.*(epsilon.^2.*P.*sin(omega) + r.*sin(omega) - epsilon.*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3));
	T28 = d2betapardrdB + d2betapardB2.*epsilon.*sum(Bsp.*cos(mb.*omega),3) + d2betapardRdB.*epsilon.*(-(deltap.*epsilon) + cos(omega) + epsilon.^2.*Pp.*cos(omega) + epsilon.*sum(cos((-1 + ms).*omega).*Sp,3));
	T30 = -(epsilon.^2.*P.*cos(omega)) - r.*cos(omega) + epsilon.*sum(-((-1 + ms).^2.*cos((-1 + ms).*omega).*S),3);
	T32 = -(epsilon.^2.*P.*sin(omega)) - r.*sin(omega) - epsilon.*sum(-((-1 + ms).^2.*S.*sin((-1 + ms).*omega)),3);
	T33 = (1 - delta.*epsilon.^2 + epsilon.^3.*P.*cos(omega) + epsilon.*r.*cos(omega) + epsilon.^2.*sum(cos((-1 + ms).*omega).*S,3)).^(-2);

	resF = (1 - delta.*epsilon.^2 + epsilon.^3.*P.*cos(omega) + epsilon.*r.*cos(omega) + epsilon.^2.*sum(cos((-1 + ms).*omega).*S,3)).*((epsilon.^2.*P.*cos(omega) + r.*cos(omega) - epsilon.*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3)).*(-(deltap.*epsilon) + cos(omega) + epsilon.^2.*Pp.*cos(omega) + epsilon.*sum(cos((-1 + ms).*omega).*Sp,3)) + (epsilon.^2.*P.*sin(omega) + r.*sin(omega) - epsilon.*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3)).*(sin(omega) + epsilon.^2.*Pp.*sin(omega) - epsilon.*sum(sin((-1 + ms).*omega).*Sp,3))).*(-(dbetapardR.*epsilon.*T3) - (r.^2.*T5.^2.*T9.*(-(epsilon.^2.*T1.*T22.*T9) + RR.*(epsilon.*T1.*T16.*T18.*T22 - epsilon.*T17.*T18.*T2.*T22 - epsilon.*T18.*T22.*T3.*T32 + epsilon.*T18.*T22.*T30.*T4 + T1.*T17.*T9 + T16.*T2.*T9 + T3.*T30.*T9 + T32.*T4.*T9)))./(q.^2.*RR.^3) + (BB.*(dbetapardr - BB.*T28 + dbetapardR.*epsilon.*T3 + dbetapardB.*epsilon.*sum(Bsp.*cos(mb.*omega),3)))./(BB - dbetapardB.*epsilon.^2) + (BB.*T5.*(BB.*(BB - dbetapardB.*epsilon.^2).*t2p + T5.*(BB.*T28 - dbetapardB.*epsilon.*sum(Bsp.*cos(mb.*omega),3))))./((BB - dbetapardB.*epsilon.^2).^3.*RR.^2) + (BB.*epsilon.^2.*T22.*(-(dbetapardR.*(BB - dbetapardB.*epsilon.^2).^2.*RR.^2.*T1) + (BB.^2.*RR.^2 - 2.*BB.*dbetapardB.*epsilon.^2.*RR.^2 + dbetapardB.^2.*epsilon.^4.*RR.^2 - T13).*(BB.*T27 - dbetapardB.*sum(-(mb.*Bs.*sin(mb.*omega)),3))))./((BB - dbetapardB.*epsilon.^2).^3.*RR.^2.*T12) + (r.*T5.*T9.*(-(qp.*r.*RR.*T12.*T5.*T9) + q.*(deltapp.*epsilon.*r.*RR.*T12.*T18.*T2.*T5 - r.*RR.*T12.*T16.*T18.*T3.*T5 + r.*RR.*T12.*T17.*T18.*T4.*T5 + epsilon.^2.*r.*RR.*T12.*t2p.*T9 + RR.*T12.*T5.*T9 + 2.*r.*RR.*T1.*T17.*T5.*T9 + 2.*r.*RR.*T16.*T2.*T5.*T9 - epsilon.*r.*T12.*T3.*T5.*T9 - epsilon.^2.*Ppp.*r.*RR.*T12.*T18.*T2.*T5.*cos(omega) + epsilon.^2.*Ppp.*r.*RR.*T1.*T12.*T18.*T5.*sin(omega) - epsilon.*r.*RR.*T12.*T18.*T2.*T5.*sum(cos((-1 + ms).*omega).*Spp,3) - epsilon.*r.*RR.*T1.*T12.*T18.*T5.*sum(sin((-1 + ms).*omega).*Spp,3))))./(q.^3.*RR.^3));
	resB = (BB - sqrt(T13.*T33.*((epsilon.^2.*r.^2.*T12.*T18)./q.^2 + (-1 + (dbetapardB.*epsilon.^2)./(1 + epsilon.*sum(Bs.*cos(mb.*omega),3))).^(-2))))./epsilon;
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
