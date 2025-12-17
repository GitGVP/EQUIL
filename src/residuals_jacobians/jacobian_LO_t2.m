function J = jacobian_LO_t2(r,x,epsilon,omega,kinetic_profiles,q,qp,P0,P1,P2,dof_count,Nb,Sbc,P_templates,M_extended,equation_of_state)
	Nq = numel(r);
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
	Bs       = reshape(P0 * Bs_flat, Nq, 1, []);
	Bsp = reshape(P1 * Bs_flat, Nq, 1, []);
	S       = reshape(P0 * S_flat, Nq, 1, []);
	Sp = reshape(P1 * S_flat, Nq, 1, []);
	Spp = reshape(P2 * S_flat, Nq, 1, []);
	ms = reshape(linspace(2,Ns+1,Ns), [1,1,Ns]);
	mb = reshape(0:(Nb-1), [1,1,Nb]);
    Ntheta = numel(omega);
    mF   = [0,1,2:(Ns+1)];
    rowsF = [1,2,(2 + Nb + (2:(Ns+1)))];
    mB   = 0:(Nb-1);
    rowsB = 3 + (1:Nb); 
    nmF = numel(mF);
    nmB = numel(mB);
	jacTotal = zeros(3,nRes,nRes,Nq); % (which deriv, which equation, which profile)

	% --- resF modes ---
	j12F = cos(omega).*ones(Nq,1);
	j12F_FFT = real(fft(j12F, [], 2)) / Ntheta;
	jacTotal(1,2, rowsF, :) = reshape(j12F_FFT(:, mF+1).', [1,1,nmF, Nq]);

	j21F = q.^3.*ones(size(omega));
	j21F_FFT = real(fft(j21F, [], 2)) / Ntheta;
	jacTotal(2,1, rowsF, :) = reshape(j21F_FFT(:, mF+1).', [1,1,nmF, Nq]);

	% derivatives wrt S_k (bundle), S modes indexed k=2..Ns+1
	for kp = 2:(Ns+1)
		j33kpF = cos(kp.*omega).*ones(Nq,1);
		j33kpF_FFT = real(fft(j33kpF, [], 2)) / Ntheta;
		jacTotal(3, 3 + Nb + (kp - 1), rowsF, :) = reshape(j33kpF_FFT(:, mF+1).', [1,1,nmF, Nq]);
	end

	% --- resP ---
	gamma = 3;
	jacTotal(1,3, gamma, :) = 1;

	% --- resB modes ---
	% derivatives wrt B_l (bundle)
	for lp = 0:Nb-1
		j13lpB = epsilon.*cos(lp.*omega).*ones(Nq,1);
		j13lpB_FFT = real(fft(j13lpB, [], 2)) / Ntheta;
		jacTotal(1, 3+lp+1, rowsB, :) = reshape(j13lpB_FFT(:, mB+1).', [1,1,nmB, Nq]);
	end

	%% Construct J
    totalDofs = dof_count * (3 + Nb) + (dof_count-1) *Ns;
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
