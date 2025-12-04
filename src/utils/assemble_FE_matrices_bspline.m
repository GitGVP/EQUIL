function [r_q, ...
          P0_full, P1_full, P2_full,...
          M_profiles, M_extended, P_templates, P_extended, ...
          A_global, profile_lengths, profile_starts, P0_end] = ...
          assemble_FE_matrices_bspline(r_nodes, Nb, Ns, nq, p)

% r_nodes vector, Nb number of "extra full" profiles, Ns number with right-BC
m = (numel(r_nodes)-1)/2;


m_s = m;                    % number of spline spans must equal FE elements
nb = m_s + p;               % control points = spans + degree
active_dofs = nb - 1;       % DOFs after eliminating DOF at r=0


% quadrature (we expect caller set nq; here just prepare)
[gauss_pts, gauss_wts] = lgwt(nq, 0, 1);
gauss_pts = flip(gauss_pts);

Nq = nq * m_s;
dof_count = active_dofs;

% estimate nonzeros: per quad point there are (p+1) contributions
est_nnz = Nq * (p+1);

% triplets for P and M (we only store columns for active DOFs: col = node_global-1)
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

    r0 = T(e+p);         % start of element e
    r1 = T(e+p+1);       % end of element e
    elem_len = r1 - r0;

    r_eval_e = r0 + elem_len * gauss_pts(:)';  
    [Nall_e, N1all_e, N2all_e] = bspline_eval_all(T, p, r_eval_e);
    %[Nall_e, N1all_e, N2all_e] = bspline_eval_all([zeros(1,p+1),(1:(m_s-1))/m_s,ones(1,p+1)], p, r_eval_e);
    for q = 1:nq
        q_global = (e-1)*nq + q;
        r_q(q_global) = r_eval_e(q);
        for i = 1:(p+1)
            node_global = nodes(i);
            if node_global == 1
                continue                       % skip DOF1; active col = node_global-1
            end
            col = node_global - 1;
            phi_iq   = Nall_e(node_global, q);
            dphi_iq  = N1all_e(node_global, q);
            d2phi_iq = N2all_e(node_global, q);
            P0_i(idx) = q_global; P0_j(idx) = col; P0_v(idx) = phi_iq;
            P1_i(idx) = q_global; P1_j(idx) = col; P1_v(idx) = dphi_iq;
            P2_i(idx) = q_global; P2_j(idx) = col; P2_v(idx) = d2phi_iq;
            M_i(idx)  = col;     M_j(idx)  = q_global; M_v(idx)  = elem_len * phi_iq * gauss_wts(q);
            idx = idx + 1;
        end
    end
end

valid = 1:idx-1;
P0_full = sparse(P0_i(valid), P0_j(valid), P0_v(valid), Nq, dof_count);
P1_full = sparse(P1_i(valid), P1_j(valid), P1_v(valid), Nq, dof_count);
P2_full = sparse(P2_i(valid), P2_j(valid), P2_v(valid), Nq, dof_count);
M_full  = sparse(M_i(valid),  M_j(valid),  M_v(valid),  dof_count, Nq);

% reduced S-blocks (profiles with an extra right Dirichlet BC have one fewer DOF)
P0_S = P0_full(:, 1:(dof_count-1));
P1_S = P1_full(:, 1:(dof_count-1));
P2_S = P2_full(:, 1:(dof_count-1));
M_S  = M_full(1:(dof_count-1), :);

% profile bookkeeping
nProfiles = 3 + Nb + Ns;
profile_lengths = [repmat(dof_count, 1, 3 + Nb), repmat(dof_count - 1, 1, Ns)];
totalDofs = sum(profile_lengths);
profile_starts = cumsum([1, profile_lengths(1:end-1)]);
profile_ends = profile_starts + profile_lengths - 1;

% build M_profiles (maps flattened per-profile residuals -> stacked DOFs)
[ii_M_full, jj_M_full, vv_M_full] = find(M_full);
[ii_M_S, jj_M_S, vv_M_S] = find(M_S);
nnz_M_full = numel(vv_M_full);
nnz_M_S = numel(vv_M_S);

% M_profiles: columns are [res(1,:)', res(2,:)', ...] so dimension totalDofs x (nProfiles*Nq)
total_Mp_triplets = (3 + Nb) * nnz_M_full + Ns * nnz_M_S;
Mp_i = zeros(total_Mp_triplets,1);
Mp_j = zeros(total_Mp_triplets,1);
Mp_v = zeros(total_Mp_triplets,1);
ptrM = 1;
for alpha = 1:nProfiles
    rows_alpha = profile_starts(alpha):profile_ends(alpha);
    if alpha <= 3 + Nb
        ii = ii_M_full; jj = jj_M_full; vv = vv_M_full; block_nnz = nnz_M_full;
    else
        ii = ii_M_S; jj = jj_M_S; vv = vv_M_S; block_nnz = nnz_M_S;
    end
    col_block_start = (alpha-1) * Nq;
    slice = ptrM:(ptrM + block_nnz - 1);
    Mp_i(slice) = rows_alpha(ii)';
    Mp_j(slice) = col_block_start + jj;
    Mp_v(slice) = vv;
    ptrM = ptrM + block_nnz;
end
M_profiles = sparse(Mp_i, Mp_j, Mp_v, totalDofs, nProfiles * Nq);

% Build M_extended and P_extended templates by block replication (same pattern as you had)
nCols = nProfiles^2 * Nq;
total_M_triplets = nProfiles * ( (3 + Nb) * nnz_M_full + Ns * nnz_M_S );
M_tri_i = zeros(total_M_triplets,1);
M_tri_j = zeros(total_M_triplets,1);
M_tri_v = zeros(total_M_triplets,1);
ptr = 1;
for alpha = 1:nProfiles
    rows_alpha = profile_starts(alpha):profile_ends(alpha);
    if alpha <= 3 + Nb
        ii = ii_M_full; jj = jj_M_full; vv = vv_M_full; block_nnz = nnz_M_full;
    else
        ii = ii_M_S; jj = jj_M_S; vv = vv_M_S; block_nnz = nnz_M_S;
    end
    for gamma = 1:nProfiles
        col_block_start = ((alpha-1)*nProfiles + (gamma-1)) * Nq;
        sliceIdx = ptr:(ptr + block_nnz - 1);
        M_tri_i(sliceIdx) = rows_alpha(ii)';
        M_tri_j(sliceIdx) = col_block_start + jj;
        M_tri_v(sliceIdx) = vv;
        ptr = ptr + block_nnz;
    end
end
M_extended = sparse(M_tri_i, M_tri_j, M_tri_v, totalDofs, nCols);

% P_extended templates and P_templates for quick scaling in Newton
P_templates = cell(1,3);
P_extended = cell(1,3);
P_blocks_full = {P0_full, P1_full, P2_full};
P_blocks_S    = {P0_S,    P1_S,    P2_S};

for d = 1:3
    [ii_P, jj_P, vv_P]   = find(P_blocks_full{d});
    [ii_PS, jj_PS, vv_PS] = find(P_blocks_S{d});
    nnz_P_full = numel(vv_P);
    nnz_P_S    = numel(vv_PS);
    total_Pd_triplets = nProfiles * ( (3 + Nb) * nnz_P_full + Ns * nnz_P_S );

    P_tri_i = zeros(total_Pd_triplets,1);
    P_tri_j = zeros(total_Pd_triplets,1);
    P_tri_v = zeros(total_Pd_triplets,1);
    ptrP = 1;

    for gamma = 1:nProfiles
        cols_gamma = profile_starts(gamma):profile_ends(gamma);
        if gamma <= 3 + Nb
            ii = ii_P; jj = jj_P; vv = vv_P; block_nnz = nnz_P_full;
        else
            ii = ii_PS; jj = jj_PS; vv = vv_PS; block_nnz = nnz_P_S;
        end

        for alpha = 1:nProfiles
            row_block_start = ((alpha-1)*nProfiles + (gamma-1)) * Nq;
            sliceIdx = ptrP:(ptrP + block_nnz - 1);
            P_tri_i(sliceIdx) = row_block_start + ii;
            P_tri_j(sliceIdx) = cols_gamma(jj)';
            P_tri_v(sliceIdx) = vv;
            ptrP = ptrP + block_nnz;
        end
    end

    if ptrP-1 ~= total_Pd_triplets
        error('Mismatch in P%d triplet count.', d);
    end

    P_extended{d} = sparse(P_tri_i, P_tri_j, P_tri_v, nCols, totalDofs);
    P_templates{d}.i = P_tri_i;
    P_templates{d}.j = P_tri_j;
    P_templates{d}.v_template = P_tri_v;
end

% store M template
M_template.i = M_tri_i;
M_template.j = M_tri_j;
M_template.v = M_tri_v;

% build A_global as block-diagonal of local A blocks (A_full and A_S)
A_full = M_full * P0_full;         % dof_count x dof_count
A_S    = M_S    * P0_S;            % (dof_count-1) x (dof_count-1)
[ii_A_full, jj_A_full, vv_A_full] = find(A_full);
[ii_A_S, jj_A_S, vv_A_S] = find(A_S);
nnz_A_full = numel(vv_A_full);
nnz_A_S = numel(vv_A_S);

n_full_profiles = 3 + Nb;
n_S_profiles = Ns;
total_A_nnz = n_full_profiles * nnz_A_full + n_S_profiles * nnz_A_S;
Ai = zeros(total_A_nnz,1); Aj = zeros(total_A_nnz,1); Av = zeros(total_A_nnz,1);
ptrA = 1;
for alpha = 1:nProfiles
    block_rows = profile_starts(alpha):profile_ends(alpha);
    if alpha <= n_full_profiles
        slice = ptrA:(ptrA + nnz_A_full - 1);
        Ai(slice) = block_rows(ii_A_full)';
        Aj(slice) = block_rows(jj_A_full)';
        Av(slice) = vv_A_full;
        ptrA = ptrA + nnz_A_full;
    else
        slice = ptrA:(ptrA + nnz_A_S - 1);
        Ai(slice) = block_rows(ii_A_S)';
        Aj(slice) = block_rows(jj_A_S)';
        Av(slice) = vv_A_S;
        ptrA = ptrA + nnz_A_S;
    end
end
A_global = sparse(Ai, Aj, Av, totalDofs, totalDofs);

% return last-column P0 for BC assembly convenience
P0_end = P0_full(:, end);

end
