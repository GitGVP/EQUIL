function LY = equilY(L, LX)
  r_fine = linspace(0,1,0.5*1e4); 
  r_fine = r_fine(2:end);
  if true %isequal(L.P.residuals_fun, @residuals_rotation) || isequal(L.P.residuals_fun, @residuals_vec)
      [dbetapardB0, dbetapardr, dbetapardB, dbetapardR, d2betapardB2, ...
    d2betapardrdB, d2betapardRdB, d2betapardrdR, d2betapardR2, d3betapardrdB2, ...
    d3betapardB3, d3betapardBdR2, d3betapardB2dR, d3betapardrdRdB, d3betapardrdR2, betapar, ...
    betaperp] = L.P.equation_of_state(LX.kinetic_profiles,r_fine.',ones(numel(r_fine), numel(L.omega)),ones(numel(r_fine), numel(L.omega)));
      dbetapardr = mean(dbetapardr,2).';
      d2betapardrdB = mean(d2betapardrdB,2).';
      d3betapardrdB2 = mean(d3betapardrdB2,2).';
      d3betapardrdRdB = mean(d3betapardrdRdB,2).';
      d2betapardrdR = mean(d2betapardrdR,2).';
      d3betapardrdR2 = mean(d3betapardrdR2,2).';
  else
      [dbetapardr, ~, ~, d2betapardrdB] = L.P.equation_of_state(LX.betafun(r_fine.'),LX.betapfun(r_fine.'),...
        LX.Ahfun(r_fine.'),LX.Ahpfun(r_fine.'),ones(numel(r_fine), numel(L.omega)),zeros(numel(r_fine), numel(L.omega)));
      dbetapardr = mean(dbetapardr,2).';
      d2betapardrdB = mean(d2betapardrdB,2).';
  end
  % made anisotropic
  LY.t2p_fine = - 2 * r_fine ./ LX.qfun(r_fine) .^ 2 + r_fine .^ 2 .* LX.qpfun(r_fine) ./ LX.qfun(r_fine) .^ 3 - dbetapardr;
  LY.t2_fine =  cumtrapz(r_fine, LY.t2p_fine);
  LY.t2_ana = interp1(r_fine, LY.t2_fine, L.r_q, 'spline');
  LY.t2p_ana = interp1(r_fine, LY.t2p_fine, L.r_q, 'spline');

  LY.deltap_fine = LX.qfun(r_fine).^2 ./ r_fine .^ 3 .* cumtrapz(r_fine, r_fine.^3 ./ LX.qfun(r_fine).^2 + r_fine.^2 .* (-2 * dbetapardr - d2betapardrdR + d2betapardrdB));
  LY.delta_fine = cumtrapz(r_fine, LY.deltap_fine);
  LY.delta_ana = interp1(r_fine, LY.delta_fine, L.r_q, 'spline');
  LY.deltap_ana = interp1(r_fine, LY.deltap_fine, L.r_q, 'spline');


  [LY.S2_fine, LY.S2p_fine] = solve_S_eq(r_fine, LX.S2bc,LX.qfun,LX.qpfun, 2);
  [LY.S3_fine, LY.S3p_fine] = solve_S_eq(r_fine, LX.S3bc,LX.qfun,LX.qpfun, 3);
  LY.P_fine = -r_fine.^3/8 + LY.S2_fine .^2 ./ (2* r_fine) - r_fine .* LY.delta_fine / 2;  % todo: add other shaping
  LY.Pp_fine = -3*r_fine.^2/8 + 2 * LY.S2_fine .* LY.S2p_fine  ./ (2* r_fine) - LY.S2_fine.^2 ./ (2* r_fine.^2) -  LY.delta_fine / 2 - r_fine .* LY.deltap_fine / 2;
  LY.S2_ana = interp1(r_fine, LY.S2_fine, L.r_q, 'spline');
  LY.S3_ana = interp1(r_fine, LY.S3_fine, L.r_q, 'spline');
  LY.P_ana = interp1(r_fine, LY.P_fine, L.r_q, 'spline');
  LY.Pp_ana = interp1(r_fine, LY.Pp_fine, L.r_q, 'spline');

  %NLO shift
  delta1_RHS = (LY.S2_fine.*(-3.*(-d2betapardrdB + d2betapardrdR + 2.*dbetapardr).*LX.qfun(r_fine).^3.*r_fine + 4.*LX.qpfun(r_fine).*r_fine.*LY.S3_fine - 2.*LX.qfun(r_fine).*(-3.*LY.deltap_fine.*r_fine + 14.*LY.S3_fine + 5.*r_fine.*LY.S3p_fine)) + r_fine.*LY.S2p_fine.*(-3.*(-d2betapardrdB + d2betapardrdR + 2.*dbetapardr).*LX.qfun(r_fine).^3.*r_fine + 6.*LX.qpfun(r_fine).*r_fine.^2.*(LY.deltap_fine - LY.S3p_fine) + 2.*LX.qfun(r_fine).*(-10.*LY.S3_fine + r_fine.*(-3.*LY.deltap_fine + r_fine + LY.S3p_fine))))./(2..*LX.qfun(r_fine).^3);
  LY.delta1p_ana = LX.qfun(r_fine).^2 ./ r_fine .^ 3 .* cumtrapz(r_fine, delta1_RHS);
  LY.delta1_ana = cumtrapz(r_fine, LY.delta1p_ana);
  
  %NLO elongation
  S2_1_RHS = (LX.qfun(r_fine).^2./r_fine.^3) .* (r_fine.*(6.*LY.deltap_fine.*LX.qpfun(r_fine).*r_fine.^2.*(LY.deltap_fine - 2.*LY.S3p_fine) + 2.*LX.qfun(r_fine).*(-6.*LY.deltap_fine.^2.*r_fine + r_fine.^3 + 2.*LY.deltap_fine.*(r_fine.^2 - 8.*LY.S3_fine) - 2.*r_fine.*(-2.*LY.deltap_fine + r_fine).*LY.S3p_fine) + LX.qfun(r_fine).^3.*(6.*d2betapardrdB.*LY.deltap_fine.*r_fine - 2.*d2betapardrdB.*r_fine.^2 + d3betapardrdR2.*r_fine.^2 + d3betapardrdB2.*r_fine.^2 - 2.*d3betapardrdRdB.*r_fine.^2 - 8.*d2betapardrdB.*LY.S3_fine - 6.*d2betapardrdB.*r_fine.*LY.S3p_fine + d2betapardrdR.*(-6.*LY.deltap_fine.*r_fine + 4.*r_fine.^2 + 8.*LY.S3_fine + 6.*r_fine.*LY.S3p_fine) + 2.*dbetapardr.*(8.*LY.S3_fine + r_fine.*(-6.*LY.deltap_fine + r_fine + 6.*LY.S3p_fine)))))./(4..*LX.qfun(r_fine).^3);
  S2_1_RHS_fun = @(r) interp1(r_fine, S2_1_RHS, r, 'spline');
  [LY.S2_1_fine, LY.S2_1p_fine] = solve_S_eq(r_fine, 0,LX.qfun,LX.qpfun, 2, S2_1_RHS_fun);
  
  LY.r_fine = r_fine;
  LY.B1_ana =  - L.r_q;
  
  profiles_q = zeros(L.Nq * (3 + L.P.Nb + L.P.Ns), 1);

  % Calculate segment sizes
  base_offset = 0;

  % Assign base variables (t2, delta, P)
  profiles_q(base_offset + (1:L.Nq)) = LY.t2_ana;
  base_offset = base_offset + L.Nq;
  profiles_q(base_offset + (1:L.Nq)) = LY.delta_ana;
  base_offset = base_offset + L.Nq;
  profiles_q(base_offset + (1:L.Nq)) = LY.P_ana;
  base_offset = base_offset + L.Nq;

  % Assign B1
  profiles_q(base_offset + (1:L.Nq)) = 0;%LY.B1_ana;
  base_offset = base_offset + L.Nq;

  % Remaining B coefficients (B2 through B_Nb) are already zeros from preallocation
  base_offset = base_offset + (L.P.Nb - 1) * L.Nq;

  % Assign S coefficients
  profiles_q(base_offset + (1:L.Nq)) = LY.S2_ana;
  base_offset = base_offset + L.Nq;

  if L.P.Ns >= 2
    profiles_q(base_offset + (1:L.Nq)) = LY.S3_ana;
    base_offset = base_offset + L.Nq;
    
    if L.P.Ns > 2
        for i = 3:L.P.Ns
            profiles_q(base_offset + (1:L.Nq)) = LX.Sbc(i) * L.r_q.^(i-1);
            base_offset = base_offset + L.Nq;
        end
    end
  end

  if not(L.P.hot_restart)
    % testing 0 for all profiles besides shaping
    profiles_q(1:(3+L.P.Nb)*L.Nq) = 0;
    x = profiles_to_x(profiles_q, L.A_global, L.M_profiles, L.P0_end, L.P.Nb, LX.Sbc, L.profile_lengths, L.Nq);
  else
      x = LX.x;
  end
  res_norms = [];LY.isconverged = false;
  % Newton iterations
  for k = 1:L.P.nk
    %Compute residuals and Jacobian at quadrature points
    if true%isequal(L.P.residuals_fun, @residuals_rotation) || isequal(L.P.residuals_fun, @residuals_vec)
        residuals = L.P.residuals_fun(L.r_q,x,LX.eps_val,L.omega,...
            LX.kinetic_profiles,LX.q_vec,LX.qp_vec, L.P0,L.P1,L.P2,L.M_profiles,L.dof_count,L.P.Nb,LX.Sbc,L.P.equation_of_state);
        J = L.P.jacobian_fun(L.r_q,x,LX.eps_val,L.omega,...
            LX.kinetic_profiles,LX.q_vec,LX.qp_vec, L.P0, L.P1, L.P2,L.dof_count,L.P.Nb,LX.Sbc,L.P_templates, L.M_extended,L.P.equation_of_state);
        update = J \ residuals;
    else
        residuals = L.P.residuals_fun(L.r_q,x,LX.eps_val,L.omega,...
            LX.beta_vec,LX.betap_vec,LX.q_vec,LX.qp_vec,LX.Ah_vec,LX.Ahp_vec, L.P0,L.P1,L.P2,L.M_profiles,L.dof_count,L.P.Nb,LX.Sbc,L.P.equation_of_state);
        J = L.P.jacobian_fun(L.r_q,x,LX.eps_val,L.omega,...
            LX.beta_vec,LX.betap_vec,LX.q_vec,LX.qp_vec,LX.Ah_vec,LX.Ahp_vec, L.P0, L.P1, L.P2,L.dof_count,L.P.Nb,LX.Sbc,L.P_templates, L.M_extended,L.P.equation_of_state);

        if isequal(L.P.jacobian_fun, @jacobian_noEL) || isequal(L.P.jacobian_fun, @jacobian_noEL_2)
            update = J \ residuals;
        else
            [LL,UU,PP,QQ] = lu(J,'vector');  % P,Q permutations for stability
            precond = @(v) apply_preconditioner(v,LL,UU,PP,QQ);
            if isequal(L.P.residuals_fun, @residuals_noEL)
              dv = 1e-8; % seems to work well.
              % Complex arithmetic doesn't work with FFT so whatever.

              Jv = @(v) (L.P.residuals_fun(L.r_q,x+dv*v,LX.eps_val,L.omega,...
              LX.beta_vec,LX.betap_vec,LX.q_vec,LX.qp_vec,LX.Ah_vec,LX.Ahp_vec, L.P0,L.P1,L.P2,L.M_profiles,L.dof_count,L.P.Nb,LX.Sbc,L.P.equation_of_state) - ...
              L.P.residuals_fun(L.r_q,x-dv*v,LX.eps_val,L.omega,...
              LX.beta_vec,LX.betap_vec,LX.q_vec,LX.qp_vec,LX.Ah_vec,LX.Ahp_vec, L.P0,L.P1,L.P2,L.M_profiles,L.dof_count,L.P.Nb,LX.Sbc,L.P.equation_of_state)) / (2*dv);
            else
              dv = 1e-10;
              Jv = @(v) imag(L.P.residuals_fun(L.r_q,x+1i*dv*v,LX.eps_val,L.omega,...
              LX.beta_vec,LX.betap_vec,LX.q_vec,LX.qp_vec,LX.Ah_vec,LX.Ahp_vec, L.P0,L.P1,L.P2,L.M_profiles,L.dof_count,L.P.Nb,LX.Sbc,L.P.equation_of_state)) / dv;
            end
            gmres_tol = 1e-2;
            [update, ~] = gmres(Jv, residuals, [],gmres_tol, 100, precond);
        end
    end
    %update = gmres(Jv, residuals, [],gmres_tol, 100);
    % try with damping
    x = x - L.P.damping*update;
    res_norm = norm(residuals);
    res_norms(k) = res_norm;
    if L.P.debug > 1
      fprintf('Iter %d, eps %.1e, |res| = %.4e, |Δx| = %.3e\n', ...
              k, LX.eps_val, res_norm, norm(update));
    end
    % Break if converged
    if res_norm < 5e-14
        LY.isconverged = true;
      break;
    end
  end
  LY.x = x;
  LY.res_norms = res_norms;
  LY = equilPP(L,LX,LY);
end

function y = apply_preconditioner(v,L,U,P,Q)
    % Apply LU with permutations as preconditioner solve
    Q_inv = zeros(size(Q));
    Q_inv(Q) = 1:length(Q);
    
    % Apply permutation P on input
    w = v(P);      
    
    % Solve LU system
    z = U \ (L \ w);
    
    % Apply inverse permutation of Q
    y = zeros(size(z));
    y(Q) = z;
end



function x = profiles_to_x(profiles_q, A_global, M_profiles, P0_end, Nb, Sbc, profile_lengths, Nq)
% Robust conversion + solve.
% profiles_q : (nProfiles*Nq) x 1  (ordering: profile1 all q, profile2 all q, ...)
% profile_lengths : numeric vector of dofs per profile (will be coerced)

% 1) coerce profile_lengths to numeric
if ~isnumeric(profile_lengths)
    if iscell(profile_lengths)
        try
            profile_lengths = cell2mat(profile_lengths);
        catch
            error('profile_lengths is a cell that cannot be converted to numeric.');
        end
    else
        try
            profile_lengths = double(profile_lengths);
        catch
            error('profile_lengths must be numeric. Got %s.', class(profile_lengths));
        end
    end
end
profile_lengths = profile_lengths(:).';    % row vector

% 2) sizes and basic checks
nProfiles = numel(profile_lengths);
totalDofs = sum(profile_lengths);
if size(M_profiles,1) ~= totalDofs
    error('M_profiles has %d rows but totalDofs = %d.', size(M_profiles,1), totalDofs);
end
if numel(profiles_q) ~= nProfiles * Nq
    error('profiles_q length (%d) != nProfiles*Nq (%d).', numel(profiles_q), nProfiles * Nq);
end
if numel(P0_end) ~= Nq
    error('P0_end length (%d) != Nq (%d).', numel(P0_end), Nq);
end
if size(A_global,1) ~= totalDofs || size(A_global,2) ~= totalDofs
    error('A_global must be square totalDofs x totalDofs. Found %dx%d.', size(A_global,1), size(A_global,2));
end

% 3) build BC_flat in same per-profile block order used by M_profiles
BC_flat = zeros(nProfiles * Nq, 1);
Ns = numel(Sbc);
if Ns > 0
    first_S_idx = nProfiles - Ns + 1;
    for kk = 1:Ns
        profile_idx = first_S_idx + (kk-1);
        slot = (profile_idx-1)*Nq + (1:Nq);
        BC_flat(slot) = Sbc(kk) * P0_end(:);
    end
end

% 4) compute rhs and solve
rhs = M_profiles * (profiles_q - BC_flat);

% ensure A_global is sparse for efficient solve
if ~issparse(A_global)
    A_global = sparse(A_global);
end

x = A_global \ rhs;
end
