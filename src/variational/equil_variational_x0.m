function x = equil_variational_x0(L, LX)
%EQUIL_VARIATIONAL_X0 Smooth large-aspect-ratio initial profiles.

    r = L.r_q;
    epsilon = LX.eps_val;
    q0 = LX.qfun(0);
    beta0 = LX.kinetic_profiles.beta(0);
    if L.P.Ns > 0
        S2 = LX.Sbc(1);
    else
        S2 = 0;
    end

    t2_target = (beta0-(1+2*S2^2*epsilon^2)/q0^2)*r.^2;
    root = sqrt(1+epsilon^2*S2^2);
    delta_den = 8*(-2*epsilon*S2-3*epsilon^3*S2^3 ...
        + root+3*epsilon^2*S2^2*root);
    delta_num = epsilon*(-1+4*beta0*q0^2)*S2 ...
        +(1+4*beta0*q0^2)*root;
    delta_target = (delta_num/delta_den)*r.^2;
    P_target = ((root-1)/epsilon^2)*r;

    x = zeros(L.total_dofs, 1);
    x(block_rows(L,1)) = L.B0{1}\t2_target;
    x(block_rows(L,2)) = L.B0{2}\delta_target;
    x(block_rows(L,3)) = L.B0{3}\P_target;
    for is = 1:L.P.Ns
        profile = 3+is;
        target = LX.Sbc(is)*r.^is;
        target = target-L.Sbc0{is}*LX.Sbc(is);
        x(block_rows(L,profile)) = L.B0{profile}\target;
    end
    for iv = 1:L.P.Nh
        profile = 3+L.P.Ns+iv;
        target = LX.Vbc(iv)*r.^iv;
        target = target-L.Vbc0{iv}*LX.Vbc(iv);
        x(block_rows(L,profile)) = L.B0{profile}\target;
    end
end

function rows = block_rows(L, profile)
    first = L.profile_starts(profile);
    rows = first:first+L.profile_lengths(profile)-1;
end
