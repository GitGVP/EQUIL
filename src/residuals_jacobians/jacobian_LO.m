function J = jacobian_LO(r,x,epsilon,omega,kinetic_profiles,q,qp,P0,P1,P2,dof_count,Nb,Sbc,P_templates,M_extended,equation_of_state)
    Nq = numel(r);
    Ns = numel(Sbc);
    nRes = 3+Nb+Ns;
    [t2, t2p, delta, deltap, deltapp, P, Pp, Ppp, Bs, Bsp, S, Sp, Spp] = ...
      equil_unpack_state(x, dof_count, Nb, Sbc, P0, P1, P2);
	ms = reshape(linspace(2,Ns+1,Ns), [1,1,Ns]);
	mb = reshape(0:(Nb-1), [1,1,Nb]);
	BB = 1 + epsilon.*sum(Bs.*cos(mb.*omega),3);
	RR = 1 - delta.*epsilon.^2 + epsilon.^3.*P.*cos(omega) + epsilon.*r.*cos(omega) + epsilon.^2.*sum(cos((-1 + ms).*omega).*S,3);
[dbetapardB0, ...
        dbetapardr, dbetapardB, ~, ... % dbetapardR not used
        d2betapardB2, d2betapardrdB, d2betapardRdB,d2betapardrdR,~,... %d2betapardR2 not used
		d3betapardrdB2,d3betapardB3,d3betapardBdR2,d3betapardB2dR,d3betapardrdRdB] = equation_of_state(kinetic_profiles,r,RR,BB);
    d2betapardrdR=0;
    Ntheta = numel(omega);
    mF   = [0,1,2:(Ns+1)];
    rowsF = [1,2,(2 + Nb + (2:(Ns+1)))];
    mB   = 0:(Nb-1);
    rowsB = 3 + (1:Nb); 
    nmF = numel(mF);
    nmB = numel(mB);
	jacTotal = zeros(3,nRes,nRes,Nq); % (which deriv, which equation, which profile)

	% --- resF modes ---
	j21F = -(cos(omega).*sum((-1 + ms).*cos((-1 + ms).*omega).*S,3)) - sin(omega).*sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3) - r.*(-1 + deltap.*cos(omega) + r.*cos(omega) - cos(omega).*sum(cos((-1 + ms).*omega).*Sp,3) + sin(omega).*sum(sin((-1 + ms).*omega).*Sp,3));
	j21F_FFT = real(fft(j21F, [], 2)) / Ntheta;
	jacTotal(2,1, rowsF, :) = reshape(j21F_FFT(:, mF+1).', [1,1,nmF, Nq]);

	j22F = -((r.*(-(q.*r) + qp.*r.^2 + q.^3.*(dbetapardr + t2p)).*cos(omega))./q.^3);
	j22F_FFT = real(fft(j22F, [], 2)) / Ntheta;
	jacTotal(2,2, rowsF, :) = reshape(j22F_FFT(:, mF+1).', [1,1,nmF, Nq]);

	j32F = (r.^3.*cos(omega))./q.^2;
	j32F_FFT = real(fft(j32F, [], 2)) / Ntheta;
	jacTotal(3,2, rowsF, :) = reshape(j32F_FFT(:, mF+1).', [1,1,nmF, Nq]);

	% derivatives wrt B_l (bundle)
	for lp = 0:Nb-1
		j13lpF = -(d2betapardrdB.*r.*cos(lp.*omega));
		j13lpF_FFT = real(fft(j13lpF, [], 2)) / Ntheta;
		jacTotal(1, 3+lp+1, rowsF, :) = reshape(j13lpF_FFT(:, mF+1).', [1,1,nmF, Nq]);
		j23lpF = -(dbetapardB.*r.*cos(lp.*omega));
		j23lpF_FFT = real(fft(j23lpF, [], 2)) / Ntheta;
		jacTotal(2, 3+lp+1, rowsF, :) = reshape(j23lpF_FFT(:, mF+1).', [1,1,nmF, Nq]);
	end

	% derivatives wrt S_k (bundle), S modes indexed k=2..Ns+1
	for kp = 2:(Ns+1)
		j13kpF = (qp.*r.^2.*(cos(omega).*(-1 + kp).*cos((-1 + kp).*omega).*1 + sin(omega).*-((-1 + kp).*1.*sin((-1 + kp).*omega))))./q.^3 - (dbetapardr + t2p).*(cos(omega).*(-1 + kp).*cos((-1 + kp).*omega).*1 + sin(omega).*-((-1 + kp).*1.*sin((-1 + kp).*omega))) + (-(r.*cos(omega).*-((-1 + kp).^2.*cos((-1 + kp).*omega).*1)) + r.*sin(omega).*-((-1 + kp).^2.*1.*sin((-1 + kp).*omega)))./q.^2;
		j13kpF_FFT = real(fft(j13kpF, [], 2)) / Ntheta;
		jacTotal(1, 3 + Nb + (kp - 1), rowsF, :) = reshape(j13kpF_FFT(:, mF+1).', [1,1,nmF, Nq]);

		j23kpF = (r.*(-(q.*r) + qp.*r.^2 + q.^3.*(dbetapardr + t2p)).*(cos(omega).*cos((-1 + kp).*omega).*1 - sin(omega).*1.*sin((-1 + kp).*omega)))./q.^3;
		j23kpF_FFT = real(fft(j23kpF, [], 2)) / Ntheta;
		jacTotal(2, 3 + Nb + (kp - 1), rowsF, :) = reshape(j23kpF_FFT(:, mF+1).', [1,1,nmF, Nq]);

		j33kpF = -((r.^3.*(cos(omega).*cos((-1 + kp).*omega).*1 - sin(omega).*1.*sin((-1 + kp).*omega)))./q.^2);
		j33kpF_FFT = real(fft(j33kpF, [], 2)) / Ntheta;
		jacTotal(3, 3 + Nb + (kp - 1), rowsF, :) = reshape(j33kpF_FFT(:, mF+1).', [1,1,nmF, Nq]);

	end

	% --- resP ---
	gamma = 3;
	jacTotal(1,2, gamma, :) = -r;
	jacTotal(1,3, gamma, :) = -1;
	jacTotal(2,2, gamma, :) = mean(-(r.^2.*cos(omega).^2) - sum((-1 + ms).*cos((-1 + ms).*omega).*S,3), 2);
	jacTotal(2,3, gamma, :) = -r;
	for kp = 2:(Ns+1)
		jacTotal(1, 3 + Nb + (kp - 1), gamma, :) = mean(r.*cos((-1 + kp).*omega).*1 - (-1 + kp).*cos((-1 + kp).*omega).*1.*(deltap + r.*cos(omega).^2 - sum(cos((-1 + ms).*omega).*Sp,3)) - -((-1 + kp).*1.*sin((-1 + kp).*omega)).*(r.*cos(omega).*sin(omega) + sum(sin((-1 + ms).*omega).*Sp,3)), 2);
		jacTotal(2, 3 + Nb + (kp - 1), gamma, :) = mean(cos((-1 + kp).*omega).*1.*(r.^2.*cos(omega).^2 + sum((-1 + ms).*cos((-1 + ms).*omega).*S,3)) - 1.*sin((-1 + kp).*omega).*(r.^2.*cos(omega).*sin(omega) + sum(-((-1 + ms).*S.*sin((-1 + ms).*omega)),3)), 2);
	end

	% --- resB modes ---
	% derivatives wrt B_l (bundle)
	for lp = 0:Nb-1
		j13lpB = 2.*cos(lp.*omega).* ones(Nq,1);
		j13lpB_FFT = real(fft(j13lpB, [], 2)) / Ntheta;
		jacTotal(1, 3+lp+1, rowsB, :) = reshape(j13lpB_FFT(:, mB+1).', [1,1,nmB, Nq]);
	end

	%% Construct J
    totalDofs = dof_count * (2 + Nb) + (dof_count-1) *(Ns+1);
    nCols = (3+ Nb + Ns)^2 * Nq;
    % precomputed: M_extended (totalDofs x nCols), P_templates{d} with fields i,j,v_template
    J = spalloc(totalDofs, totalDofs, 1 * nnz(M_extended));  % or a tighter nnz estimate

    for d = 1:3
        % Extract jacobian block and flatten in the SAME order P_templates was built.
        jac_d = squeeze(jacTotal(d,:,:,:));    % size: [nProfiles, nProfiles, Nq] (alpha,gamma,q)
        % permute to [q, gamma, alpha] then linearize -> ordering matches P_templates{d}.i
        jacTotal_flat = reshape(permute(jac_d, [3, 1, 2]), [], 1);   % length == nCols

        % Row-scale P template by jac values (this gives diag(jac)*P_extended)
        vP_scaled = P_templates{d}.v_template .* jacTotal_flat(P_templates{d}.i);

        % Create sparse scaled-P and accumulate J via M_extended * P_scaled
        P_scaled = sparse(P_templates{d}.i, P_templates{d}.j, vP_scaled, nCols, totalDofs);
        J = J + M_extended * P_scaled;
    end
end
