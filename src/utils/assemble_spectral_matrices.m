function [r_q, ...
          P0_full, P1_full, P2_full,...
          M_profiles, M_extended, P_templates, P_extended, ...
          A_global, profile_lengths, profile_starts, P0_end] = ...
          assemble_spectral_matrices(r_nodes, Nb, Ns, nq, p)

% r_nodes vector, Nb number of "extra full" profiles, Ns number with right-BC
m = (numel(r_nodes)-1)/2;


m_s = m;                    % number of spline spans must equal FE elements
%nb = m_s + p;               % control points = spans + degree
%active_dofs = nb - 1;       % DOFs after eliminating DOF at r=0

% total number of global nodes for SEM with p-degree (p+1 per element, shared endpoints)
Nnodes = m_s * p + 1;
active_dofs = Nnodes - 1;   % remove DOF at r=0 (global node 1)


% quadrature (we expect caller set nq; here just prepare)
%[gauss_pts, gauss_wts] = lgwt(nq, 0, 1);
%gauss_pts = flip(gauss_pts);
% default quadrature points: use LGL nodes. If user supplied nq, accept it,
% but for full SEM accuracy typically nq == p+1 (GLL).
%if nargin < 4 || isempty(nq)
    nq = p + 1;
%end

% Get LGL nodes & weights on [-1,1] for quadrature; map to [0,1]
[xi_lgl, w_lgl] = legendre_gauss_lobatto(nq-1); % returns nq nodes and weights
gauss_pts = (xi_lgl + 1) / 2;      % in [0,1]
gauss_wts = w_lgl * 0.5;           % weights on [0,1]


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
%loc_offset = 0:p;
% create element partition using the same clustering parameter as original FE code
alpha_clustering = 1;
internal_knots = ((1:(m_s-1))/m_s).^alpha_clustering;
T = [zeros(1,p+1), internal_knots, ones(1,p+1)];  % breakpoints in [0,1]

% Reference nodal (LGL) positions for local polynomial basis (p+1 nodes)
[xi_ref, ~] = legendre_gauss_lobatto(p);    % p+1 nodes in [-1,1]
s_nodes_ref = (xi_ref + 1)/2;               % mapped to [0,1]

% Precompute barycentric weights for the nodal Lagrange basis on s_nodes_ref
w_bary = compute_barycentric_weights(s_nodes_ref);

% Build differentiation matrix D at nodal LGL points (p+1 x p+1).
% D_ij gives derivative of basis j evaluated at node i (in s coordinate).
D = compute_D_matrix(s_nodes_ref, w_bary);
D2 = D * D; % second derivative at nodes

% Loop elements and assemble
for e = 1:m_s
    % element interval in [0,1] param
    r0 = T(e+p);         % start of element e
    r1 = T(e+p+1);       % end of element e
    elem_len = r1 - r0;

    % physical evaluation points in this element
    r_eval_e = r0 + elem_len * gauss_pts(:)';   % 1 x nq

    % local reference coordinates s_eval in [0,1]
    s_eval = (r_eval_e - r0) / elem_len;       % 1 x nq

    % local nodal positions in s coordinate (p+1 nodes)
    s_nodes = s_nodes_ref;                     % column vector length p+1

    % Evaluate basis and derivatives at s_eval using barycentric formulas
    % This returns (p+1) x nq matrices: Nall_e, N1all_e_s, N2all_e_s
    [Nall_e, N1all_e_s, N2all_e_s] = barycentric_basis_and_derivatives(s_nodes, w_bary, s_eval, D, D2);

    % apply chain rule to map s derivatives -> r derivatives
    N1all_e = N1all_e_s / elem_len;
    N2all_e = N2all_e_s / (elem_len^2);

    % Global node indices for element (shared endpoints):
    nodes = ( (e-1)*p + (1:(p+1)) );

    % loop quadrature points and fill triplets
    for q = 1:nq
        q_global = (e-1)*nq + q;
        r_q(q_global) = r_eval_e(q);

        for i_local = 1:(p+1)
            node_global = nodes(i_local);
            if node_global == 1
                continue                       % skip DOF at r=0 (eliminated)
            end
            col = node_global - 1;            % active column index (shifted by -1)

            phi_iq   = Nall_e(i_local, q);
            dphi_iq  = N1all_e(i_local, q);
            d2phi_iq = N2all_e(i_local, q);

            P0_i(idx) = q_global; P0_j(idx) = col; P0_v(idx) = phi_iq;
            P1_i(idx) = q_global; P1_j(idx) = col; P1_v(idx) = dphi_iq;
            P2_i(idx) = q_global; P2_j(idx) = col; P2_v(idx) = d2phi_iq;

            % M: dof_count x Nq  (row = dof, col = q_global)
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

% profile bookkeeping (unchanged)
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

% Build M_extended and P_extended templates by block replication (same pattern)
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

% P_extended templates and P_templates for quick scaling in Newton (unchanged)
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

%% -------------------------
%% Helper subfunctions below
%% -------------------------
function w = compute_barycentric_weights(xs)
% compute barycentric weights w_j = 1 / prod_{k!=j} (xs(j) - xs(k))
n = numel(xs);
w = ones(n,1);
for j = 1:n
    diff = xs(j) - xs([1:j-1, j+1:n]);
    w(j) = 1 / prod(diff);
end
end

function D = compute_D_matrix(xs, w)
n = numel(xs);
D = zeros(n,n);
for i = 1:n
    for j = 1:n
        if i ~= j
            D(i,j) = w(j) / (w(i) * (xs(i) - xs(j)));
        end
    end
    D(i,i) = -sum(D(i, [1:i-1, i+1:n]));
end
end

function [Lvals, Lprime, L2] = barycentric_basis_and_derivatives(nodes, w, x_eval, D, D2)
nodes = nodes(:);
N = numel(nodes);
x_eval = x_eval(:);
M = numel(x_eval);

Lvals  = zeros(N, M);
Lprime = zeros(N, M);
L2     = zeros(N, M);

tol = 1e-14;

for ix = 1:M
    x = x_eval(ix);
    diffs = x - nodes;
    hits = find(abs(diffs) < tol);
    if ~isempty(hits)
        k = hits(1);
        Lvals(:,ix) = 0; Lvals(k,ix) = 1;
        Lprime(:,ix) = D(k, :)';
        L2(:,ix) = (D2(k, :))';
        continue
    end

    c = w ./ diffs;
    c1 = - w ./ (diffs.^2);
    c2 = 2 * w ./ (diffs.^3);

    v = sum(c);
    v1 = sum(c1);
    v2 = sum(c2);

    Lvals(:,ix) = c / v;
    num1 = c1 * v - c * v1;
    Lprime(:,ix) = num1 / (v^2);

    num = num1;
    nump = c2 * v - c * v2;
    L2(:,ix) = (nump * v - 2 * num * v1) / (v^3);
end
end

function [x, w] = legendre_gauss_lobatto(N)
if N == 0
    x = 0; w = 2; return
end
m = N;
k = (0:m)';
x = -cos(pi * k / m);
for iter = 1:100
    [Pvals, Pder] = legendreP_and_derivative(N, x);
    eps_fd = 1e-8;
    % compute derivative of P_N'(x) by central finite difference of the derivative:
    % call legendreP_and_derivative to obtain Pder at x+eps and x-eps (second output)
    [~, Pder_plus] = legendreP_and_derivative(N, x + eps_fd);
    [~, Pder_minus] = legendreP_and_derivative(N, x - eps_fd);
    Psec = (Pder_plus - Pder_minus) / (2*eps_fd);
    dx = zeros(size(x));
    interior = (abs(x) < 0.9999999999);
    % safe-guard: avoid division by zero in Psec; only update interior nodes where Psec is nonzero
    safe = interior & (abs(Psec) > eps);
    dx(safe) = - Pder(safe) ./ Psec(safe);
    x = x + dx;
    x(1) = -1; x(end) = 1;
    if max(abs(dx)) < 1e-14, break; end
end

% compute weights: w_i = 2 / (N*(N+1) * [P_N(x_i)]^2 )
[Pvals, ~] = legendreP_and_derivative(N, x);
w = 2 ./ (N * (N + 1) * (Pvals.^2));
end

function [Pvals, Pder] = legendreP_and_derivative(n, x)
x = x(:);
N = numel(x);
if n == 0
    Pvals = ones(N,1);
    Pder = zeros(N,1);
    return
end
Pnm2 = ones(N,1);
Pnm1 = x;
if n == 1
    Pvals = Pnm1;
    % derivative of P1(x) = 1
    Pder = ones(N,1);
    return
end
for k = 2:n
    Pn = ( (2*k-1) .* x .* Pnm1 - (k-1) * Pnm2 ) / k;
    Pnm2 = Pnm1;
    Pnm1 = Pn;
end
Pvals = Pn;
% derivative via recurrence: P_n'(x) = n/(x^2-1) ( x P_n(x) - P_{n-1}(x) )
Pder = ( n * ( x .* Pvals - Pnm2 ) ) ./ (x.^2 - 1);
end