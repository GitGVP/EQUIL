function R = residuals_LO(r,x,epsilon,omega,kinetic_profiles,q,qp,P0,P1,P2,M_profiles,dof_count,Nb,Sbc,equation_of_state)
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

[~, ...
        dbetapardr, dbetapardB, ~, ...
        ~, d2betapardrdB, ~,d2betapardrdR,~,...
		~,~,~,~,~] = equation_of_state(kinetic_profiles,r,ones(N,numel(omega)),ones(N,numel(omega)));

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
