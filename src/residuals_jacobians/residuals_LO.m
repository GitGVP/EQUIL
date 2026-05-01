function R = residuals_LO(r,x,epsilon,omega,kinetic_profiles,q,qp,P0,P1,P2,M_profiles,dof_count,Nb,Sbc,equation_of_state)
    N = numel(r);
    Ns = numel(Sbc);
    nRes = 3+Nb+Ns;
    [t2, t2p, delta, deltap, deltapp, P, Pp, Ppp, Bs, Bsp, S, Sp, Spp] = ...
          equil_unpack_state(x, dof_count, Nb, Sbc, P0, P1, P2);
    resTotal = zeros(nRes,N);
	ms = reshape(linspace(2,Ns+1,Ns), [1,1,Ns]);
	mb = reshape(0:(Nb-1), [1,1,Nb]);
	RR = 1 - delta.*epsilon.^2 + epsilon.^3.*P.*cos(omega) + epsilon.*r.*cos(omega) + epsilon.^2.*sum(cos((-1 + ms).*omega).*S,3);
	BB = 1 + epsilon.*sum(Bs.*cos(mb.*omega),3);

	[dbetapardB0, ...
        dbetapardr, dbetapardB, dbetapardR, ...
        d2betapardB2, d2betapardrdB, d2betapardRdB] = equation_of_state(kinetic_profiles,r,RR,BB);
    d2betapardrdR=0;
	resF = (qp.*r.^2.*(cos(omega).*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3) + sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) + r.*(-1 - deltap.*cos(omega) + r.*cos(omega) + cos(omega).*sum(cos((-1 + ms).*omega).*Sp,3) - sin(omega).*sum(sin((-1 + ms).*omega).*Sp,3))) - q.^3.*(-(dbetapardr.*r) + dbetapardr.*deltap.*r.*cos(omega) + 2.*d2betapardrdB.*r.^2.*cos(omega) - d2betapardrdR.*r.^2.*cos(omega) - dbetapardr.*r.^2.*cos(omega) + d2betapardrdB.*r.*sum(Bs.*cos(mb.*omega),3) + dbetapardB.*r.*(cos(omega) + sum(Bsp.*cos(mb.*omega),3)) + dbetapardr.*cos(omega).*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3) + dbetapardr.*sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) - dbetapardr.*r.*cos(omega).*sum(cos((-1 + ms).*omega).*Sp,3) + dbetapardr.*r.*sin(omega).*sum(sin((-1 + ms).*omega).*Sp,3) + t2p.*(cos(omega).*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3) + sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) + r.*(-1 + deltap.*cos(omega) + r.*cos(omega) - cos(omega).*sum(cos((-1 + ms).*omega).*Sp,3) + sin(omega).*sum(sin((-1 + ms).*omega).*Sp,3)))) - q.*r.*(cos(omega).*sum(-((-1 + ms).^2.*cos((-1 + ms).*omega).*S),3) - sin(omega).*sum(-((-1 + ms).^2.*S.*sin((-1 + ms).*omega)),3) + r.*(-2 - deltap.*cos(omega) + 3.*r.*cos(omega) - deltapp.*r.*cos(omega) + cos(omega).*sum(cos((-1 + ms).*omega).*Sp,3) - sin(omega).*sum(sin((-1 + ms).*omega).*Sp,3) + r.*cos(omega).*sum(cos((-1 + ms).*omega).*Spp,3) - r.*sin(omega).*sum(sin((-1 + ms).*omega).*Spp,3))))./q.^3;
	resB = 2.*(r.*cos(omega) + sum(Bs.*cos(mb.*omega),3));
	resF_FFT = real(fft(resF, [], 2))/numel(omega);
	resB_FFT = real(fft(resB, [], 2))/numel(omega);
	mF = [0, 1, 2:(Ns+1)];
	rowsF = [1, 2, (2 + Nb + (2:(Ns+1)))];
	mB = 0:(Nb-1);
	rowsB = 3 + (1:Nb);
	resTotal(rowsF, :) = resF_FFT(:, mF + 1).';
	resTotal(rowsB, :) = resB_FFT(:, mB + 1).';
	resTotal(3,:) = mean(-(P.*cos(omega).^2) - Pp.*r.*cos(omega).^2 - P.*sin(omega).^2 - Pp.*r.*sin(omega).^2 - r.*(delta + r.^2.*cos(omega).^2 - sum(cos((-1 + ms).*omega).*S,3)) + sum((-1 + ms).*cos((-1 + ms).*omega).*S,3).*(-deltap + sum(cos((-1 + ms).*omega).*Sp,3)) - sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3).*sum(sin((-1 + ms).*omega).*Sp,3) - r.*cos(omega).*(cos(omega).*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3) + sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) + r.*(deltap.*cos(omega) - cos(omega).*sum(cos((-1 + ms).*omega).*Sp,3) + sin(omega).*sum(sin((-1 + ms).*omega).*Sp,3))), 2);
	R = M_profiles * reshape(resTotal.', [], 1);
end