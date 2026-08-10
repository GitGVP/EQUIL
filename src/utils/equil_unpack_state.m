function [t2, t2p, delta, deltap, deltapp, P, Pp, Ppp, Bs, Bsp, S, Sp, Spp] = ...
          equil_unpack_state(x, dof_count, Nb, Sbc, P0, P1, P2)

Ns = numel(Sbc);
idx = 1;

if iscell(P0)
    % Axis-regular basis: each profile owns its projection matrices.
    % Block lengths are inferred from those matrices, allowing each
    % profile to keep only the powers it can actually use.
    n_t2 = size(P0{1},2);
    n_delta = size(P0{2},2);
    n_P = size(P0{3},2);
    t2_c = x(idx:idx+n_t2-1);          idx = idx + n_t2;
    delta_c = x(idx:idx+n_delta-1);    idx = idx + n_delta;
    P_c = x(idx:idx+n_P-1);            idx = idx + n_P;

    t2      = P0{1} * t2_c;
    t2p     = P1{1} * t2_c;
    delta   = P0{2} * delta_c;
    deltap  = P1{2} * delta_c;
    deltapp = P2{2} * delta_c;
    P       = P0{3} * P_c;
    Pp      = P1{3} * P_c;
    Ppp     = P2{3} * P_c;

    Nq = size(P0{1},1);
    Bs  = zeros(Nq,1,Nb);
    Bsp = zeros(Nq,1,Nb);
    for ib = 1:Nb
        alpha = 3 + ib;
        n_B = size(P0{alpha},2);
        B_c = x(idx:idx+n_B-1);        idx = idx + n_B;
        Bs(:,:,ib)  = P0{alpha} * B_c;
        Bsp(:,:,ib) = P1{alpha} * B_c;
    end

    S   = zeros(Nq,1,Ns);
    Sp  = zeros(Nq,1,Ns);
    Spp = zeros(Nq,1,Ns);
    for is = 1:Ns
        alpha = 3 + Nb + is;
        n_S_unknown = size(P0{alpha},2)-1;
        S_c = [x(idx:idx+n_S_unknown-1); Sbc(is)];
        idx = idx + n_S_unknown;
        S(:,:,is)   = P0{alpha} * S_c;
        Sp(:,:,is)  = P1{alpha} * S_c;
        Spp(:,:,is) = P2{alpha} * S_c;
    end
    if idx-1 ~= numel(x)
        error('Axis-regular state layout consumed %d of %d entries.', ...
              idx-1,numel(x));
    end
else
    % Original open-knot B-spline basis.
    t2_c = x(idx:idx+dof_count-1);             idx = idx + dof_count;
    delta_c = x(idx:idx+dof_count-2);          idx = idx + (dof_count-1);
    delta_cf = [0; delta_c];
    P_c = x(idx:idx+dof_count-1);              idx = idx + dof_count;
    Bs_c = reshape(x(idx:idx+Nb*dof_count-1), ...
                   [dof_count, 1, Nb]);        idx = idx + Nb*dof_count;
    S_unkn = reshape(x(idx:end), [dof_count-1, 1, Ns]);

    S_c = cat(1, S_unkn, reshape(Sbc, [1, 1, Ns]));
    Bs_flat = reshape(Bs_c, size(Bs_c,1), []);
    S_flat  = reshape(S_c,  size(S_c,1),  []);

    t2      = P0 * t2_c;
    t2p     = P1 * t2_c;
    delta   = P0 * delta_cf;
    deltap  = P1 * delta_cf;
    deltapp = P2 * delta_cf;
    P       = P0 * P_c;
    Pp      = P1 * P_c;
    Ppp     = P2 * P_c;
    Bs      = reshape(P0 * Bs_flat, size(P0,1), 1, []);
    Bsp     = reshape(P1 * Bs_flat, size(P1,1), 1, []);
    S       = reshape(P0 * S_flat, size(P0,1), 1, []);
    Sp      = reshape(P1 * S_flat, size(P1,1), 1, []);
    Spp     = reshape(P2 * S_flat, size(P2,1), 1, []);
end

end
