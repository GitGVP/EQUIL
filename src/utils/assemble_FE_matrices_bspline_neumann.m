function [r_q, ...
          P0_full, P1_full, P2_full,...
          M_profiles, M_extended, P_templates, P_extended, ...
          A_global, profile_lengths, profile_starts, P0_end, ...
          ML2, ML2_q] = ...
          assemble_FE_matrices_bspline_neumann(r_nodes, Nb, Ns, nq, p)

m = (numel(r_nodes)-1)/2;
m_s = m;
nb = m_s + p;
active_dofs = nb - 1;       % dof_count: after removing c_1 (Dirichlet at r=0)

[gauss_pts, gauss_wts] = lgwt(nq, 0, 1);
gauss_pts = flip(gauss_pts);

Nq = nq * m_s;
dof_count = active_dofs;

est_nnz = Nq * (p+1);
P0_i = zeros(est_nnz,1); P0_j = zeros(est_nnz,1); P0_v = zeros(est_nnz,1);
P1_i = zeros(est_nnz,1); P1_j = zeros(est_nnz,1); P1_v = zeros(est_nnz,1);
P2_i = zeros(est_nnz,1); P2_j = zeros(est_nnz,1); P2_v = zeros(est_nnz,1);
M_i  = zeros(est_nnz,1); M_j  = zeros(est_nnz,1); M_v  = zeros(est_nnz,1);

idx = 1;
r_q = zeros(Nq,1);
loc_offset = 0:p;

for e = 1:m_s
    nodes = e + loc_offset;
    alpha_clustering = 1;
    internal_knots = ((1:(m_s-1))/m_s).^alpha_clustering;
    T = [zeros(1,p+1), internal_knots, ones(1,p+1)];

    r0 = T(e+p);
    r1 = T(e+p+1);
    elem_len = r1 - r0;

    r_eval_e = r0 + elem_len * gauss_pts(:)';
    [Nall_e, N1all_e, N2all_e] = bspline_eval_all(T, p, r_eval_e);

    for q = 1:nq
        q_global = (e-1)*nq + q;
        r_q(q_global) = r_eval_e(q);
        for i = 1:(p+1)
            node_global = nodes(i);
            if node_global == 1
                continue
            end
            col = node_global - 1;
            phi_iq   = Nall_e(node_global, q);
            dphi_iq  = N1all_e(node_global, q);
            d2phi_iq = N2all_e(node_global, q);
            P0_i(idx) = q_global; P0_j(idx) = col; P0_v(idx) = phi_iq;
            P1_i(idx) = q_global; P1_j(idx) = col; P1_v(idx) = dphi_iq;
            P2_i(idx) = q_global; P2_j(idx) = col; P2_v(idx) = d2phi_iq;
            M_i(idx)  = col;      M_j(idx)  = q_global; M_v(idx) = elem_len * phi_iq * gauss_wts(q);
            idx = idx + 1;
        end
    end
end

valid = 1:idx-1;
P0_full = sparse(P0_i(valid), P0_j(valid), P0_v(valid), Nq, dof_count);
P1_full = sparse(P1_i(valid), P1_j(valid), P1_v(valid), Nq, dof_count);
P2_full = sparse(P2_i(valid), P2_j(valid), P2_v(valid), Nq, dof_count);
M_full  = sparse(M_i(valid),  M_j(valid),  M_v(valid),  dof_count, Nq);

% S-block: right Dirichlet BC removes last DOF
P0_S = P0_full(:, 1:(dof_count-1));
P1_S = P1_full(:, 1:(dof_count-1));
P2_S = P2_full(:, 1:(dof_count-1));
M_S  = M_full(1:(dof_count-1), :);

% Delta-block: both c_1=0 and c_2=0 enforced, removes first two DOFs
% (c_1 already removed by col=node_global-1 shift; c_2 corresponds to col=1)
P0_D = P0_full(:, 2:end);      % drop col 1  (was c_2, now first active col)
P1_D = P1_full(:, 2:end);
P2_D = P2_full(:, 2:end);
M_D  = M_full(2:end, :);       % drop row 1

% profile bookkeeping
% layout: [t2 | delta | P | Bs(1..Nb) | S(1..Ns)]
% t2, P, Bs: dof_count DOFs each  (3+Nb profiles)
% delta:     dof_count-2 DOFs     (1 profile, index 2)
% S:         dof_count-1 DOFs     (Ns profiles)
nProfiles = 3 + Nb + Ns;
profile_lengths = [dof_count, dof_count-1, dof_count, ...   % t2, delta, P
                   repmat(dof_count, 1, Nb), ...             % Bs
                   repmat(dof_count-1, 1, Ns)];              % S
totalDofs = sum(profile_lengths);
profile_starts = cumsum([1, profile_lengths(1:end-1)]);
profile_ends   = profile_starts + profile_lengths - 1;

% --- helper: return the right (ii,jj,vv,block_nnz) for profile alpha ---
[ii_M_full, jj_M_full, vv_M_full] = find(M_full);
[ii_M_S,    jj_M_S,    vv_M_S   ] = find(M_S);
[ii_M_D,    jj_M_D,    vv_M_D   ] = find(M_D);
nnz_M_full = numel(vv_M_full);
nnz_M_S    = numel(vv_M_S);
nnz_M_D    = numel(vv_M_D);

get_M_block = @(alpha) get_block(alpha, Nb, Ns, ...
    ii_M_full, jj_M_full, vv_M_full, nnz_M_full, ...
    ii_M_D,    jj_M_D,    vv_M_D,    nnz_M_D, ...
    ii_M_S,    jj_M_S,    vv_M_S,    nnz_M_S);

% M_profiles
total_Mp_nnz = nnz_M_full + nnz_M_D + (1+Nb)*nnz_M_full + Ns*nnz_M_S;
% recount properly:
total_Mp_nnz = 0;
for alpha = 1:nProfiles
    [~,~,~,bnnz] = get_M_block(alpha);
    total_Mp_nnz = total_Mp_nnz + bnnz;
end
Mp_i = zeros(total_Mp_nnz,1);
Mp_j = zeros(total_Mp_nnz,1);
Mp_v = zeros(total_Mp_nnz,1);
ptrM = 1;
for alpha = 1:nProfiles
    [ii,jj,vv,bnnz] = get_M_block(alpha);
    rows_alpha = profile_starts(alpha):profile_ends(alpha);
    col_block_start = (alpha-1)*Nq;
    slice = ptrM:(ptrM+bnnz-1);
    Mp_i(slice) = rows_alpha(ii)';
    Mp_j(slice) = col_block_start + jj;
    Mp_v(slice) = vv;
    ptrM = ptrM + bnnz;
end
M_profiles = sparse(Mp_i, Mp_j, Mp_v, totalDofs, nProfiles*Nq);

% M_extended
% again recount:
total_M_nnz = 0;
for alpha = 1:nProfiles
    [~,~,~,bnnz] = get_M_block(alpha);
    total_M_nnz = total_M_nnz + nProfiles*bnnz;
end
M_tri_i = zeros(total_M_nnz,1);
M_tri_j = zeros(total_M_nnz,1);
M_tri_v = zeros(total_M_nnz,1);
ptr = 1;
for alpha = 1:nProfiles
    [ii,jj,vv,bnnz] = get_M_block(alpha);
    rows_alpha = profile_starts(alpha):profile_ends(alpha);
    for gamma = 1:nProfiles
        col_block_start = ((alpha-1)*nProfiles + (gamma-1))*Nq;
        sliceIdx = ptr:(ptr+bnnz-1);
        M_tri_i(sliceIdx) = rows_alpha(ii)';
        M_tri_j(sliceIdx) = col_block_start + jj;
        M_tri_v(sliceIdx) = vv;
        ptr = ptr + bnnz;
    end
end
M_extended = sparse(M_tri_i, M_tri_j, M_tri_v, totalDofs, nProfiles^2*Nq);

% P_extended and P_templates
P_blocks_full = {P0_full, P1_full, P2_full};
P_blocks_S    = {P0_S,    P1_S,    P2_S};
P_blocks_D    = {P0_D,    P1_D,    P2_D};

P_templates = cell(1,3);
P_extended  = cell(1,3);

for d = 1:3
    [ii_P,   jj_P,   vv_P  ] = find(P_blocks_full{d});
    [ii_PS,  jj_PS,  vv_PS ] = find(P_blocks_S{d});
    [ii_PD,  jj_PD,  vv_PD ] = find(P_blocks_D{d});
    nnz_P_full = numel(vv_P);
    nnz_P_S    = numel(vv_PS);
    nnz_P_D    = numel(vv_PD);

    total_Pd_nnz = 0;
    for gamma = 1:nProfiles
        if gamma == 2                          % delta profile
            total_Pd_nnz = total_Pd_nnz + nProfiles * nnz_P_D;
        elseif gamma <= 3+Nb
            total_Pd_nnz = total_Pd_nnz + nProfiles * nnz_P_full;
        else
            total_Pd_nnz = total_Pd_nnz + nProfiles * nnz_P_S;
        end
    end

    P_tri_i = zeros(total_Pd_nnz,1);
    P_tri_j = zeros(total_Pd_nnz,1);
    P_tri_v = zeros(total_Pd_nnz,1);
    ptrP = 1;

    for gamma = 1:nProfiles
        cols_gamma = profile_starts(gamma):profile_ends(gamma);
        if gamma == 2
            ii = ii_PD; jj = jj_PD; vv = vv_PD; bnnz = nnz_P_D;
        elseif gamma <= 3+Nb
            ii = ii_P;  jj = jj_P;  vv = vv_P;  bnnz = nnz_P_full;
        else
            ii = ii_PS; jj = jj_PS; vv = vv_PS; bnnz = nnz_P_S;
        end

        for alpha = 1:nProfiles
            row_block_start = ((alpha-1)*nProfiles + (gamma-1))*Nq;
            sliceIdx = ptrP:(ptrP+bnnz-1);
            P_tri_i(sliceIdx) = row_block_start + ii;
            P_tri_j(sliceIdx) = cols_gamma(jj)';
            P_tri_v(sliceIdx) = vv;
            ptrP = ptrP + bnnz;
        end
    end

    if ptrP-1 ~= total_Pd_nnz
        error('Mismatch in P%d triplet count.', d);
    end

    P_extended{d}  = sparse(P_tri_i, P_tri_j, P_tri_v, nProfiles^2*Nq, totalDofs);
    P_templates{d}.i          = P_tri_i;
    P_templates{d}.j          = P_tri_j;
    P_templates{d}.v_template = P_tri_v;
end

% A_global (block diagonal)
A_full = M_full * P0_full;
A_S    = M_S    * P0_S;
A_D    = M_D    * P0_D;
[ii_A_full, jj_A_full, vv_A_full] = find(A_full);
[ii_A_S,    jj_A_S,    vv_A_S   ] = find(A_S);
[ii_A_D,    jj_A_D,    vv_A_D   ] = find(A_D);
nnz_A_full = numel(vv_A_full);
nnz_A_S    = numel(vv_A_S);
nnz_A_D    = numel(vv_A_D);

total_A_nnz = 0;
for alpha = 1:nProfiles
    if alpha == 2
        total_A_nnz = total_A_nnz + nnz_A_D;
    elseif alpha <= 3+Nb
        total_A_nnz = total_A_nnz + nnz_A_full;
    else
        total_A_nnz = total_A_nnz + nnz_A_S;
    end
end
Ai = zeros(total_A_nnz,1); Aj = zeros(total_A_nnz,1); Av = zeros(total_A_nnz,1);
ptrA = 1;
for alpha = 1:nProfiles
    block_rows = profile_starts(alpha):profile_ends(alpha);
    if alpha == 2
        ii = ii_A_D; jj = jj_A_D; vv = vv_A_D; bnnz = nnz_A_D;
    elseif alpha <= 3+Nb
        ii = ii_A_full; jj = jj_A_full; vv = vv_A_full; bnnz = nnz_A_full;
    else
        ii = ii_A_S; jj = jj_A_S; vv = vv_A_S; bnnz = nnz_A_S;
    end
    slice = ptrA:(ptrA+bnnz-1);
    Ai(slice) = block_rows(ii)';
    Aj(slice) = block_rows(jj)';
    Av(slice) = vv;
    ptrA = ptrA + bnnz;
end
A_global = sparse(Ai, Aj, Av, totalDofs, totalDofs);

P0_end = P0_full(:, end);

if nargout > 12
    X    = A_full \ M_full;
    ML2  = M_full' * X;
    w_q  = sum(M_full,1)';
    ML2_q = spdiags(w_q, 0, Nq, Nq);
end
end

% -------------------------------------------------------------------------
function [ii, jj, vv, bnnz] = get_block(alpha, Nb, Ns, ...
        ii_full, jj_full, vv_full, nnz_full, ...
        ii_D,    jj_D,    vv_D,    nnz_D, ...
        ii_S,    jj_S,    vv_S,    nnz_S)
    if alpha == 2
        ii = ii_D; jj = jj_D; vv = vv_D; bnnz = nnz_D;
    elseif alpha <= 3+Nb
        ii = ii_full; jj = jj_full; vv = vv_full; bnnz = nnz_full;
    else
        ii = ii_S; jj = jj_S; vv = vv_S; bnnz = nnz_S;
    end
end