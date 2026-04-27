function x = equil_x0(L, LX, profiles_q)
    if nargin < 3
        profiles_q = zeros(L.Nq * (3 + L.P.Nb + L.P.Ns), 1);

        %% Initial guess TODO: (rotating, anisotropic cases)
        % The order is such that we solve the uncoupled equations first:
        % 1) P, 2) t2, 3) delta
        % test initial guess: linear P
        profiles_q((2*L.Nq+1:3*L.Nq)) = (-1 + sqrt(1+LX.Sbc(1)^2 * LX.eps_val^2)) .* L.r_q / LX.eps_val^2; % positive branch
        % test initial guess, quadratic t2
        beta0 = LX.kinetic_profiles.beta(0);
        % needs to not be singular!
        q0 = LX.qfun(0);
        % This uses P0 above, but is independent of the branch selected
        % above
        profiles_q((1:L.Nq)) = (beta0 - (1 + 2 *LX.Sbc(1)^2 * LX.eps_val^2)/q0^2) * L.r_q.^2;

        % test initial guess: quadratic delta
        %profiles_q((L.Nq+1:2*L.Nq)) = (q0 + 4 * beta0 * q0^3 - 7 * L.r_q * qp0) ./ (8*q0 - 10 * L.r_q *qp0) .* L.r_q.^2;
        profiles_q((L.Nq+1:2*L.Nq)) = (LX.eps_val*(-1 + 4*beta0*q0^2)*LX.Sbc(1) + (1 + 4*beta0*q0^2)*sqrt(1 + ...
        LX.eps_val^2*LX.Sbc(1)^2))/(8.*(-2*LX.eps_val*LX.Sbc(1) - 3*LX.eps_val^3*LX.Sbc(1)^3 + ...
        sqrt(1 + LX.eps_val^2*LX.Sbc(1)^2) + 3*LX.eps_val^2*LX.Sbc(1)^2*sqrt(1 + ...
        LX.eps_val^2*LX.Sbc(1)^2))) .* L.r_q.^2;


        base_offset = (L.P.Nb + 3) * L.Nq;
        for i = 1:L.P.Ns
            profiles_q(base_offset + (1:L.Nq)) = LX.Sbc(i) * L.r_q.^(i);
            base_offset = base_offset + L.Nq;
        end
    end
    base_offset = (L.P.Nb + 3) * L.Nq;
    for i = 1:L.P.Ns
        profiles_q(base_offset + (1:L.Nq)) = profiles_q(base_offset + (1:L.Nq)) - LX.Sbc(i) * L.P0_end;
        base_offset = base_offset + L.Nq;
    end
    rhs = L.M_profiles * profiles_q;
    x   = L.lu.Q * (L.lu.U \ (L.lu.L \ (L.lu.P * rhs)));
end