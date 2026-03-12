function x = equil_x0(L, LX, profiles_q)
    if nargin < 3
        profiles_q = zeros(L.Nq * (3 + L.P.Nb + L.P.Ns), 1);
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